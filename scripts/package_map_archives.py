#!/usr/bin/env python3

"""
Wrap already-built lookup maps into release assets

bin/build_id_name_maps.py builds a map from a GTF and packages it in one step. This
does only the packaging half, for maps that already exist as plain pickles -- which
in practice means the releases that were committed into ngsutils/data before the
maps moved to release assets, recovered from git:

  git show <revision>:ngsutils/data/id2name_gencode.v45.pkl > id2name_gencode.v45.pkl

The archives are written the same way the builder writes them, gzip with the header
timestamp zeroed, so the same input always produces the same bytes and the SHA-256
recorded here can be checked later. A metadata file is written beside them in the
manifest's shape, so a release of an older map describes itself.

NOTE: what this cannot do is say where the map came from. A pickle carries no
  provenance, so the source of a recovered map is the git blob it was read out of
  and nothing more. Pass that revision so the metadata can record it.

Usage:
  package_map_archives.py --release <STR> --id2name <PATH> --name2id <PATH> [--recovered-from <REV>] [--tag <STR>] [--repository <SLUG>] [--output-dir <PATH>]

Options:
  --release <STR>          : Release label, for instance gencode.v45 (required)
  --id2name <PATH>         : The identifier to name pickle (required)
  --name2id <PATH>         : The name to identifier pickle (required)
  --recovered-from <REV>   : Git revision the pickles were read out of [default: unknown]
  --tag <STR>              : Git tag the assets go under [default: auto]
  --repository <SLUG>      : The owner/name the assets are fetched from [default: auto]
  --output-dir <PATH>      : Directory the archives are written into [default: auto]

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

from docopt import docopt

REPOSITORY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Kept the same as the builder's, so an archive from either is the same bytes.
GZIP_LEVEL = 6

REPOSITORY_FROM_REMOTE = re.compile(r"github\.com[:/]([^/]+/[^/]+?)(?:\.git)?$")


# %%
# Helper functions

def log(message):
    print(message, file=sys.stderr, flush=True)


def resolve_repository(requested):
    if requested != "auto":
        return requested
    completed = subprocess.run(
        ["git", "-C", REPOSITORY_ROOT, "remote", "get-url", "origin"],
        capture_output=True, text=True, check=True)
    match = REPOSITORY_FROM_REMOTE.search(completed.stdout.strip())
    if not match:
        raise SystemExit("Cannot read owner/name off the git remote. Pass it explicitly.")
    return match.group(1)


def write_archive(path, payload):
    """The builder's archive format: gzip, level 6, header timestamp zeroed."""
    temporary_path = f"{path}.tmp"
    with open(temporary_path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=GZIP_LEVEL,
                           fileobj=raw, mtime=0) as f:
            f.write(payload)
    os.replace(temporary_path, path)
    return os.path.getsize(path)


def package(kind, pickle_path, release, output_dir):
    with open(pickle_path, "rb") as f:
        payload = f.read()

    mapping = pickle.loads(payload)
    if not isinstance(mapping, dict) or not mapping:
        raise SystemExit(f"{pickle_path} does not hold a non-empty dict")

    archive_path = os.path.join(output_dir, f"{kind}_{release}.pkl.gz")
    archive_bytes = write_archive(archive_path, payload)
    log(f"  wrote {archive_path} ({archive_bytes} bytes on the wire, "
        f"{len(payload)} unpacked, {len(mapping)} keys)")

    return {
        "file": os.path.basename(archive_path),
        "archive_bytes": archive_bytes,
        "archive_sha256": hashlib.sha256(open(archive_path, "rb").read()).hexdigest(),
        "pickle_bytes": len(payload),
        "pickle_sha256": hashlib.sha256(payload).hexdigest(),
        "keys": len(mapping),
    }


# %%
# Main logic

def main():
    options = docopt(__doc__)

    release = options["--release"]
    tag = options["--tag"] if options["--tag"] != "auto" else f"maps-{release}"
    repository = resolve_repository(options["--repository"])
    output_dir = options["--output-dir"]
    if output_dir == "auto":
        output_dir = os.path.join(REPOSITORY_ROOT, "dist", "maps", release)
    os.makedirs(output_dir, exist_ok=True)

    log(f"Release:    {release}")
    log(f"Assets:     {output_dir}  (uploaded to {repository} under tag {tag})")

    assets = {
        "id2name": package("id2name", options["--id2name"], release, output_dir),
        "name2id": package("name2id", options["--name2id"], release, output_dir),
    }

    metadata = {
        "release": release,
        "repository": repository,
        "tag": tag,
        "assets": assets,
        "provenance": {
            "packaged_by": os.path.basename(__file__),
            "recovered_from": options["--recovered-from"],
            "gtf_path": None,
            "note": ("Repackaged from an already-built pickle. The GTF it was built "
                     "from is not recorded anywhere in that pickle and is unknown."),
        },
    }
    metadata_path = os.path.join(output_dir, "maps.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2, sort_keys=True)
        f.write("\n")
    log(f"  wrote {metadata_path}")

    log("Next: upload the archives and the metadata as release assets.")
    log(f"  gh release create {tag} --repo {repository} --title {tag} \\")
    log(f"      {' '.join(sorted(os.path.join(output_dir, name) for name in os.listdir(output_dir)))}")


if __name__ == "__main__":
    main()
