#!/usr/bin/env python3

"""
Measure the two routes from a GENCODE GTF to the lookup maps

bin/build_id_name_maps.py reads the input row by row. It could instead read it
through ngsutils.gtf.read_gtf, which is what id2name and name2id do when an
annotation file is given to them, and which is why the route matters: the packaged
map has to be the map those commands would have built. This measures what each
route costs and asserts they produce the same bytes.

Each route runs in its own process and reports its own peak resident set, so the
figure is that route's own rather than the high-water mark of both.

Usage:
  bench_build_routes.py --gtf <PATH> [--route <NAME>] [--output <PATH>]

Options:
  --gtf <PATH>     : Annotation file in GTF format, gzip allowed (required)
  --route <NAME>   : Measure row, read-gtf, or both of them [default: both]
  --output <PATH>  : Where a TSV of the measurements is written [default: none]

"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import pickle
import re
import resource
import subprocess
import sys
import time

from docopt import docopt

REPOSITORY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPOSITORY_ROOT not in sys.path:
    sys.path.insert(0, REPOSITORY_ROOT)

ID_PREFIX_LENGTH = 15
GENE_COLUMNS = ("gene_id", "gene_name")
TRANSCRIPT_COLUMNS = ("transcript_id", "transcript_name")
ATTRIBUTE_COLUMNS = list(GENE_COLUMNS) + list(TRANSCRIPT_COLUMNS)
READ_GTF_FEATURES = ["gene", "transcript"]

ATTRIBUTE_PATTERNS = {
    column: re.compile(rf'{column} "([^"]*)"') for column in ATTRIBUTE_COLUMNS
}

ROUTES = ("row", "read-gtf")


# %%
# The two routes

def pairs_row_by_row(gtf_path):
    """Every row, four attributes each -- what bin/build_id_name_maps.py does."""
    gene_pairs, transcript_pairs = set(), set()
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            attributes = line.split("\t", 8)[-1]
            values = {}
            for column, pattern in ATTRIBUTE_PATTERNS.items():
                match = pattern.search(attributes)
                values[column] = match.group(1) if match else ""
            gene_pairs.add(
                (values["gene_id"][:ID_PREFIX_LENGTH], values["gene_name"]))
            transcript_pairs.add(
                (values["transcript_id"][:ID_PREFIX_LENGTH], values["transcript_name"]))
    return sorted(gene_pairs), sorted(transcript_pairs)


def pairs_read_gtf(gtf_path):
    """Through ngsutils.gtf, the route id2name and name2id take for an annotation file."""
    import ngsutils.gtf as gtfparse
    import polars as pl

    gtf = gtfparse.read_gtf(
        gtf_path, features=READ_GTF_FEATURES, usecols=ATTRIBUTE_COLUMNS)

    def pairs_of(id_column, name_column):
        frame = gtf.select([
            pl.col(id_column).str.slice(0, ID_PREFIX_LENGTH).alias("id"),
            pl.col(name_column).alias("name"),
        ]).unique()
        return sorted(frame.iter_rows())

    return pairs_of(*GENE_COLUMNS), pairs_of(*TRANSCRIPT_COLUMNS)


def collapse(pairs, key_index):
    mapping = {}
    for pair in pairs:
        mapping.setdefault(pair[key_index], pair[1 - key_index])
    return mapping


def build(route, gtf_path):
    """The measurement of one route: seconds, peak RSS, and the bytes it produced."""
    started = time.perf_counter()
    gene_pairs, transcript_pairs = (
        pairs_row_by_row(gtf_path) if route == "row" else pairs_read_gtf(gtf_path))

    id2name = collapse(gene_pairs, 0) | collapse(transcript_pairs, 0)
    name2id = collapse(gene_pairs, 1) | collapse(transcript_pairs, 1)
    elapsed = time.perf_counter() - started

    return {
        "route": route,
        "seconds": round(elapsed, 1),
        "peak_rss_gigabytes": round(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2, 2),
        "gene_pairs": len(gene_pairs),
        "transcript_pairs": len(transcript_pairs),
        "id2name_keys": len(id2name),
        "name2id_keys": len(name2id),
        "id2name_sha256": hashlib.sha256(pickle.dumps(id2name)).hexdigest(),
        "name2id_sha256": hashlib.sha256(pickle.dumps(name2id)).hexdigest(),
    }


# %%
# Main logic

def run_child(route, gtf_path):
    """One route in its own process, which reports its own RUSAGE_SELF peak."""
    argv = [sys.executable, os.path.abspath(__file__),
            "--gtf", gtf_path, "--route", route]
    completed = subprocess.run(argv, capture_output=True, text=True,
                               cwd=REPOSITORY_ROOT)
    if completed.returncode != 0:
        raise SystemExit(f"the {route} route failed: {completed.stderr.strip()}")
    return json.loads(completed.stdout)


def main():
    options = docopt(__doc__)
    gtf_path = os.path.abspath(options["--gtf"])
    if not os.path.exists(gtf_path):
        raise SystemExit(f"No such GTF: {gtf_path}")

    route = options["--route"]

    if route in ROUTES:
        # A child run: report one route as JSON on stdout and let the parent print.
        print(json.dumps(build(route, gtf_path)))
        return

    if route != "both":
        raise SystemExit(f"Unknown route {route!r}, expected one of {ROUTES} or both")

    measurements = [run_child(name, gtf_path) for name in ROUTES]

    print(f"{'route':10s} {'seconds':>8s} {'peak RSS GB':>12s} {'id2name keys':>13s} "
          f"{'id2name sha256':>16s}")
    for measurement in measurements:
        print(f"{measurement['route']:10s} {measurement['seconds']:8.1f} "
              f"{measurement['peak_rss_gigabytes']:12.2f} "
              f"{measurement['id2name_keys']:13d} "
              f"{measurement['id2name_sha256'][:16]}")

    first, second = measurements
    for key in ("id2name_sha256", "name2id_sha256"):
        if first[key] != second[key]:
            raise SystemExit(
                f"the two routes disagree on {key}: {first[key]} against {second[key]}")
    print("\nboth routes pickle to the same bytes")

    if options["--output"] != "none":
        output_path = os.path.abspath(options["--output"])
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        columns = list(measurements[0])
        with open(output_path, "w") as f:
            f.write("\t".join(columns) + "\n")
            for measurement in measurements:
                f.write("\t".join(str(measurement[c]) for c in columns) + "\n")
        print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
