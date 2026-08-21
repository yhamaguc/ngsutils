#!/usr/bin/env python3

"""
Build the id2name / name2id maps from a GENCODE GTF, as release assets

Three files are written: id2name_<release>.pkl.gz, name2id_<release>.pkl.gz, and a
manifest. The two archives are uploaded as GitHub release assets and are not
committed; the manifest is small, is committed, and is what tells id2name and
name2id which release to fetch and what its bytes must hash to.

The input is read one row at a time and only four attributes are pulled out of it,
which keeps the build inside about 0.3 GB. Reading the same file through
ngsutils.gtf.read_gtf, the way id2name and name2id do when an annotation file is
given to them, needs about 9.5 GB for the same answer, so it is not the default;
the cross-check flag below runs it and asserts the two agree. Measured 2026-08-21
by scripts/bench_build_routes.py.

Identifiers are truncated to their first 15 characters, which drops the version
suffix, because that is what both commands truncate the identifiers they look up
to. In a full GENCODE annotation the PAR_Y copies collapse onto their chrX
counterpart under that truncation, carrying the same name; the primary assembly
files have no PAR_Y entries at all. A key two rows disagree on is a hard error.

Where a key was reached by more than one row the first value in sorted order wins,
so the map does not depend on the row order of the input. Reused names are counted
in the report at the end.

Usage:
  build_id_name_maps.py --gtf <PATH> [--release <STR>] [--tag <STR>] [--repository <SLUG>] [--output-dir <PATH>] [--manifest <PATH>] [--cross-check]

Options:
  --gtf <PATH>         : Annotation file in GTF format, gzip allowed (required)
  --release <STR>      : Release label used in the file names [default: auto]
  --tag <STR>          : Git tag the assets will be uploaded under [default: auto]
  --repository <SLUG>  : The owner/name the assets are fetched from [default: auto]
  --output-dir <PATH>  : Directory the archives are written into [default: auto]
  --manifest <PATH>    : Where the committed manifest is written [default: auto]
  --cross-check        : Also build through ngsutils.gtf and assert both agree

"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import pickle
import re
import subprocess
import sys
import time

from docopt import docopt

import ngsutils

# NOTE: polars and ngsutils.gtf are imported inside cross_check() rather than here.
#   A build does not need them, and this script is meant to run on a bare checkout.


# %%
# Constants

# NOTE: id2name and name2id truncate every identifier they look up to this many
#   characters, so the keys of the maps they read have to be truncated the same way.
ID_PREFIX_LENGTH = 15

# The row types the cross-check reads. Every other row type repeats the attributes of
# the transcript row above it, so it can add no pair the row by row build did not see --
# which is what makes the cross-check sound while reading a fraction of the rows.
CROSS_CHECK_FEATURES = ["gene", "transcript"]

GENE_COLUMNS = ("gene_id", "gene_name")
TRANSCRIPT_COLUMNS = ("transcript_id", "transcript_name")
ATTRIBUTE_COLUMNS = list(GENE_COLUMNS) + list(TRANSCRIPT_COLUMNS)

# gencode.v50.primary_assembly.annotation.gtf.gz -> gencode.v50
RELEASE_FROM_FILENAME = re.compile(r"gencode\.v(\d+)", re.IGNORECASE)
# ##description: ... version 50 (Ensembl 116)
RELEASE_FROM_HEADER = re.compile(r"^##description:.*\bversion (\d+)\b", re.IGNORECASE)
# git@github.com:owner/name.git or https://github.com/owner/name.git
REPOSITORY_FROM_REMOTE = re.compile(r"github\.com[:/]([^/]+/[^/]+?)(?:\.git)?$")

ATTRIBUTE_PATTERNS = {
    column: re.compile(rf'{column} "([^"]*)"') for column in ATTRIBUTE_COLUMNS
}

# How many examples of a collapsed or duplicated key the report prints.
EXAMPLES_SHOWN = 5

# NOTE: a release asset must be under 2 GiB and there is no limit on the total size of
#   a release nor on its bandwidth (docs.github.com, "About releases", read 2026-08-21).
#   That is the whole reason the archives are assets rather than committed files, so warn
#   long before the one limit that still applies.
RELEASE_ASSET_LIMIT_BYTES = 2 * 1024 ** 3
SIZE_WARNING_BYTES = 512 * 1024 ** 2

GZIP_LEVEL = 6


# %%
# Helper functions

def log(message):
    print(message, file=sys.stderr, flush=True)


def open_gtf(gtf_path):
    if gtf_path.endswith(".gz"):
        return gzip.open(gtf_path, "rt")
    return open(gtf_path)


def read_gtf_header(gtf_path):
    header_lines = []
    with open_gtf(gtf_path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            header_lines.append(line.rstrip("\n"))
    return header_lines


def repository_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve_release(gtf_path, requested):
    """The release label, from the GTF header and the file name, which must agree."""
    if requested != "auto":
        return requested

    from_header = None
    for line in read_gtf_header(gtf_path):
        match = RELEASE_FROM_HEADER.match(line)
        if match:
            from_header = match.group(1)
            break

    match = RELEASE_FROM_FILENAME.search(os.path.basename(gtf_path))
    from_filename = match.group(1) if match else None

    if from_header and from_filename and from_header != from_filename:
        raise SystemExit(
            f"The GTF header says version {from_header} and the file name says "
            f"v{from_filename}. Pass the release label explicitly.")

    version = from_header or from_filename
    if not version:
        raise SystemExit(
            "No GENCODE version in the header or the file name. Pass the release "
            "label explicitly.")

    return f"gencode.v{version}"


def resolve_repository(requested):
    """The owner/name the assets will be fetched from, read off the git remote."""
    if requested != "auto":
        return requested

    try:
        completed = subprocess.run(
            ["git", "-C", repository_root(), "remote", "get-url", "origin"],
            capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        raise SystemExit("No git remote to read the repository from. Pass it explicitly.")

    match = REPOSITORY_FROM_REMOTE.search(completed.stdout.strip())
    if not match:
        raise SystemExit(
            f"Cannot read owner/name out of {completed.stdout.strip()!r}. "
            "Pass it explicitly.")
    return match.group(1)


def resolve_output_dir(requested):
    if requested != "auto":
        return requested
    # NOTE: dist/ is gitignored. The archives must not be committed -- that is the
    #   point of shipping them as release assets.
    return os.path.join(repository_root(), "dist", "maps")


def resolve_manifest_path(requested):
    """Where the committed manifest goes: the package's data directory."""
    if requested != "auto":
        return requested

    repository_data = os.path.join(repository_root(), "ngsutils", "data")
    if os.path.isdir(repository_data):
        return os.path.join(repository_data, "maps.json")

    return os.path.join(
        os.path.dirname(os.path.abspath(ngsutils.__file__)), "data", "maps.json")


def git_revision():
    try:
        completed = subprocess.run(
            ["git", "-C", repository_root(), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


def file_md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_pairs(gtf_path):
    """The (identifier, name) pairs of every row, as two sorted lists plus the row count.

    NOTE: this is the step that owns the shape of the maps -- the truncation length,
      which attributes are read, and which rows are read. Everything downstream,
      cross_check included, is checked against what this returns.
    """
    gene_pairs = set()
    transcript_pairs = set()
    rows = 0

    with open_gtf(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            rows += 1
            attributes = line.split("\t", 8)[-1]
            values = {}
            for column, pattern in ATTRIBUTE_PATTERNS.items():
                match = pattern.search(attributes)
                values[column] = match.group(1) if match else ""
            gene_pairs.add(
                (values["gene_id"][:ID_PREFIX_LENGTH], values["gene_name"]))
            transcript_pairs.add(
                (values["transcript_id"][:ID_PREFIX_LENGTH], values["transcript_name"]))

    return sorted(gene_pairs), sorted(transcript_pairs), rows


def collapse(pairs, key_index):
    """A dict from one side of the pairs to the other, plus the keys with more than one value.

    key_index 0 keys on the identifier, 1 on the name. The pairs arrive sorted and the
    first value wins, so the map does not depend on the row order of the input.
    """
    mapping = {}
    collided = {}
    for pair in pairs:
        key, value = pair[key_index], pair[1 - key_index]
        if key in mapping:
            collided.setdefault(key, [mapping[key]]).append(value)
        else:
            mapping[key] = value
    return mapping, collided


def report_collisions(label, collided, hard):
    """Print the keys more than one row reached. `hard` makes disagreement an error."""
    if not collided:
        log(f"  {label}: no key was reached by more than one row")
        return 0

    examples = sorted(collided.items())
    log(f"  {label}: {len(collided)} keys reached by more than one row, "
        f"{sum(len(v) for v in collided.values())} rows involved")
    for key, values in examples[:EXAMPLES_SHOWN]:
        log(f"    {key} -> {values}")
    if len(examples) > EXAMPLES_SHOWN:
        log(f"    ... and {len(examples) - EXAMPLES_SHOWN} more")

    if hard:
        raise SystemExit(
            f"{label}: rows disagree on the value of {len(collided)} keys. The "
            "15-character truncation is not safe for this annotation.")
    return len(collided)


def cross_check(gtf_path, gene_pairs, transcript_pairs):
    """Assert the pairs read row by row are the pairs ngsutils.gtf reads.

    The expensive route, and the one id2name and name2id take when an annotation file
    is given to them. Reading gene and transcript rows only is enough: read_pairs read
    every row, so a pair that appears on some other row type is already in its result,
    and equality here means the other row types added nothing.
    """
    log("Cross-check: reading the same file through ngsutils.gtf...")
    import ngsutils.gtf as gtfparse
    import polars as pl

    gtf = gtfparse.read_gtf(
        gtf_path, features=CROSS_CHECK_FEATURES, usecols=ATTRIBUTE_COLUMNS)

    missing = [c for c in ATTRIBUTE_COLUMNS if c not in gtf.columns]
    if missing:
        raise SystemExit(f"The GTF has no {missing} attribute; it is not GENCODE-like.")

    def pairs_of(id_column, name_column):
        frame = gtf.select([
            pl.col(id_column).str.slice(0, ID_PREFIX_LENGTH).alias("id"),
            pl.col(name_column).alias("name"),
        ]).unique()
        return set(frame.iter_rows())

    for label, expected, found in (
            ("gene", set(gene_pairs), pairs_of(*GENE_COLUMNS)),
            ("transcript", set(transcript_pairs), pairs_of(*TRANSCRIPT_COLUMNS))):
        only_streamed = expected - found
        only_read_gtf = found - expected
        log(f"  {label}: {len(found)} pairs through ngsutils.gtf, "
            f"{len(expected)} read row by row")
        if only_streamed or only_read_gtf:
            for pair in sorted(only_streamed)[:EXAMPLES_SHOWN]:
                log(f"    only in the row by row read: {pair}")
            for pair in sorted(only_read_gtf)[:EXAMPLES_SHOWN]:
                log(f"    only through ngsutils.gtf: {pair}")
            raise SystemExit(
                f"{label}: the two routes disagree on {len(only_streamed)} and "
                f"{len(only_read_gtf)} pairs. One of them is wrong.")

    log("  the two routes agree on every pair")


def write_archive(path, mapping):
    """Write mapping as a gzipped pickle and return what the manifest has to record."""
    payload = pickle.dumps(mapping)
    temporary_path = f"{path}.tmp"
    # NOTE: gzip.open would stamp the current time into the header, which changes the
    #   archive's SHA-256 on every rebuild of the same input. The manifest records that
    #   digest, so the archive has to be reproducible: mtime=0 and no stored file name.
    with open(temporary_path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=GZIP_LEVEL,
                           fileobj=raw, mtime=0) as f:
            f.write(payload)
    os.replace(temporary_path, path)

    archive_bytes = os.path.getsize(path)
    log(f"  wrote {path} ({archive_bytes} bytes on the wire, "
        f"{len(payload)} unpacked, {len(mapping)} keys)")

    if archive_bytes > SIZE_WARNING_BYTES:
        log(f"  WARNING: {archive_bytes / 1024 ** 2:.0f} MiB. A release asset must stay "
            f"under {RELEASE_ASSET_LIMIT_BYTES // 1024 ** 3} GiB.")

    return {
        "file": os.path.basename(path),
        "archive_bytes": archive_bytes,
        "archive_sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
        "pickle_bytes": len(payload),
        "pickle_sha256": hashlib.sha256(payload).hexdigest(),
        "keys": len(mapping),
    }


# %%
# Main logic

def main():
    options = docopt(__doc__)

    gtf_path = os.path.abspath(options["--gtf"])
    if not os.path.exists(gtf_path):
        raise SystemExit(f"No such GTF: {gtf_path}")

    release = resolve_release(gtf_path, options["--release"])
    tag = options["--tag"] if options["--tag"] != "auto" else f"maps-{release}"
    repository = resolve_repository(options["--repository"])
    output_dir = resolve_output_dir(options["--output-dir"])
    manifest_path = resolve_manifest_path(options["--manifest"])
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)

    started = time.time()
    log(f"GTF:        {gtf_path}")
    for line in read_gtf_header(gtf_path):
        log(f"            {line}")
    log(f"Release:    {release}")
    log(f"Assets:     {output_dir}  (uploaded to {repository} under tag {tag})")
    log(f"Manifest:   {manifest_path}")

    log("Reading every row, four attributes per row...")
    gene_pairs, transcript_pairs, rows = read_pairs(gtf_path)
    log(f"  {rows} rows, {len(gene_pairs)} distinct gene pairs, "
        f"{len(transcript_pairs)} distinct transcript pairs")

    if options["--cross-check"]:
        cross_check(gtf_path, gene_pairs, transcript_pairs)
    else:
        log("Cross-check skipped: nothing has compared this build against ngsutils.gtf.")

    log("Collapsing to maps...")
    # NOTE: identifier -> name must be single-valued, so disagreement is fatal. Name ->
    #   identifier is not: GENCODE reuses a handful of gene names across identifiers.
    genes_by_id, gene_id_collisions = collapse(gene_pairs, 0)
    transcripts_by_id, transcript_id_collisions = collapse(transcript_pairs, 0)
    genes_by_name, gene_name_collisions = collapse(gene_pairs, 1)
    transcripts_by_name, transcript_name_collisions = collapse(transcript_pairs, 1)

    report_collisions("gene_id truncated", gene_id_collisions, hard=True)
    report_collisions("transcript_id truncated", transcript_id_collisions, hard=True)
    gene_names_reused = report_collisions("gene_name", gene_name_collisions, hard=False)
    transcript_names_reused = report_collisions(
        "transcript_name", transcript_name_collisions, hard=False)

    maps = {
        "id2name": genes_by_id | transcripts_by_id,
        "name2id": genes_by_name | transcripts_by_name,
    }

    log("Writing...")
    assets = {}
    for name, mapping in sorted(maps.items()):
        path = os.path.join(output_dir, f"{name}_{release}.pkl.gz")
        assets[name] = write_archive(path, mapping)

    manifest = {
        "release": release,
        "repository": repository,
        "tag": tag,
        "assets": assets,
        "provenance": {
            "gtf_path": gtf_path,
            "gtf_md5": file_md5(gtf_path),
            "gtf_header": read_gtf_header(gtf_path),
            "built_by": os.path.basename(__file__),
            "git_revision": git_revision(),
            "id_prefix_length": ID_PREFIX_LENGTH,
            "build_route": "row by row",
            "rows_read": rows,
            "cross_checked_against_ngsutils_gtf": bool(options["--cross-check"]),
            "gene_pairs": len(gene_pairs),
            "transcript_pairs": len(transcript_pairs),
            "gene_names_reused": gene_names_reused,
            "transcript_names_reused": transcript_names_reused,
            "python": sys.version.split()[0],
        },
    }
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")
    log(f"  wrote {manifest_path}")

    log(f"Done in {time.time() - started:.0f} s.")
    log("Next: upload the archives as release assets, then commit the manifest.")
    log(f"  gh release create {tag} --repo {repository} --title {tag} "
        f"--notes 'id2name / name2id maps built from {release}' \\")
    log(f"      {' '.join(os.path.join(output_dir, a['file']) for a in assets.values())}")


if __name__ == "__main__":
    main()
