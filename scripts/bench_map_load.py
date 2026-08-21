#!/usr/bin/env python3

"""
Measure how long a lookup map takes to load, in three on-disk forms

The maps are shipped as gzipped pickles and cached unpacked, and this is the
measurement that decision rests on: a gzipped TSV is a quarter of the size and
inspectable, a pickle is faster. Rounds are interleaved, so machine drift lands on
all three forms equally rather than on whichever ran last.

Every form is asserted to reproduce the same dict before its timing is kept.

Unless it is switched off, the two commands are also timed end to end, together
with a bare interpreter start and the import of polars, so the load time can be
read against what the command spends elsewhere. Those runs need a resolvable map,
which means a warm cache or a vendored copy.

Usage:
  bench_map_load.py --map <PATH> [--rounds <INT>] [--output <PATH>] [--skip-commands]

Options:
  --map <PATH>      : A map to measure, a pickle or a gzipped pickle (required)
  --rounds <INT>    : How many interleaved rounds to run [default: 15]
  --output <PATH>   : Where a TSV of every timing is written [default: none]
  --skip-commands   : Do not time the two commands end to end

"""

from __future__ import annotations

import gzip
import os
import pickle
import statistics
import subprocess
import sys
import tempfile
import time

from docopt import docopt

REPOSITORY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPOSITORY_ROOT not in sys.path:
    sys.path.insert(0, REPOSITORY_ROOT)

GZIP_LEVEL = 6
COMMAND_ROUNDS = 7

# Three identifiers and three names that GENCODE has carried for years, used only to
# give the commands something to answer.
SAMPLE_IDS = ["ENSG00000141510.21", "ENSG00000012048", "ENSG00000146648"]
SAMPLE_NAMES = ["TP53", "BRCA1", "EGFR"]


# %%
# Helper functions

def read_map(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return pickle.load(f)
    with open(path, "rb") as f:
        return pickle.load(f)


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def load_tsv(path):
    with open(path) as f:
        return dict(line.split("\t") for line in f.read().splitlines())


def load_tsv_gz(path):
    with gzip.open(path, "rt") as f:
        return dict(line.split("\t") for line in f.read().splitlines())


def write_forms(mapping, directory):
    """The same mapping as a pickle, a TSV and a gzipped TSV."""
    body = "".join(f"{key}\t{value}\n" for key, value in sorted(mapping.items()))

    pickle_path = os.path.join(directory, "map.pkl")
    with open(pickle_path, "wb") as f:
        pickle.dump(mapping, f)

    tsv_path = os.path.join(directory, "map.tsv")
    with open(tsv_path, "w") as f:
        f.write(body)

    gz_path = os.path.join(directory, "map.tsv.gz")
    with gzip.open(gz_path, "wt", compresslevel=GZIP_LEVEL) as f:
        f.write(body)

    return [("pickle.load", load_pickle, pickle_path),
            ("tsv -> dict", load_tsv, tsv_path),
            ("tsv.gz -> dict", load_tsv_gz, gz_path)]


def time_forms(cases, reference, rounds):
    """Interleaved timings in milliseconds, keyed by form."""
    timings = {label: [] for label, _, _ in cases}
    for _ in range(rounds):
        for label, loader, path in cases:
            started = time.perf_counter()
            result = loader(path)
            timings[label].append((time.perf_counter() - started) * 1000)
            if result != reference:
                raise SystemExit(f"{label} did not reproduce the map")
            del result
    return timings


def time_command(argv, rounds=COMMAND_ROUNDS):
    times = []
    for _ in range(rounds):
        started = time.perf_counter()
        completed = subprocess.run(argv, capture_output=True, text=True,
                                   cwd=REPOSITORY_ROOT)
        times.append((time.perf_counter() - started) * 1000)
        if completed.returncode != 0:
            raise SystemExit(
                f"{' '.join(argv)} failed: {completed.stderr.strip()}\n"
                "The commands need a resolvable map -- a warm cache or a vendored copy.")
    return times


def command_cases(directory):
    ids_path = os.path.join(directory, "ids.txt")
    names_path = os.path.join(directory, "names.txt")
    with open(ids_path, "w") as f:
        f.write("".join(f"{value}\n" for value in SAMPLE_IDS))
    with open(names_path, "w") as f:
        f.write("".join(f"{value}\n" for value in SAMPLE_NAMES))

    return [
        ("id2name end to end",
         [sys.executable, "-m", "ngsutils.id2name", "--file", ids_path]),
        ("name2id end to end",
         [sys.executable, "-m", "ngsutils.name2id", "--file", names_path]),
        ("import polars, docopt, ngsutils.gtf",
         [sys.executable, "-c", "import polars, docopt, ngsutils.gtf"]),
        ("bare interpreter", [sys.executable, "-c", "pass"]),
    ]


def report(rows, output_path):
    print(f"{'what':38s} {'MB':>7s} {'min':>7s} {'p50':>7s} {'max':>7s}")
    for label, megabytes, times in rows:
        size = f"{megabytes:.1f}" if megabytes is not None else ""
        print(f"{label:38s} {size:>7s} {min(times):7.0f} "
              f"{statistics.median(times):7.0f} {max(times):7.0f}")

    if output_path == "none":
        return
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w") as f:
        f.write("what\tmegabytes\tround\tmilliseconds\n")
        for label, megabytes, times in rows:
            for index, value in enumerate(times, start=1):
                size = f"{megabytes:.3f}" if megabytes is not None else ""
                f.write(f"{label}\t{size}\t{index}\t{value:.3f}\n")
    print(f"\nwrote {output_path}")


# %%
# Main logic

def main():
    options = docopt(__doc__)
    rounds = int(options["--rounds"])

    reference = read_map(options["--map"])
    print(f"{options['--map']}: {len(reference)} entries, "
          f"python {sys.version.split()[0]}, {rounds} interleaved rounds\n")

    rows = []
    with tempfile.TemporaryDirectory() as directory:
        cases = write_forms(reference, directory)
        timings = time_forms(cases, reference, rounds)
        for label, _, path in cases:
            rows.append((label, os.path.getsize(path) / 1e6, timings[label]))

        if not options["--skip-commands"]:
            for label, argv in command_cases(directory):
                rows.append((label, None, time_command(argv)))

    report(rows, options["--output"])


if __name__ == "__main__":
    main()
