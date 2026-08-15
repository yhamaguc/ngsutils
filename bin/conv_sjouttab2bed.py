#!/usr/bin/env python3
"""Convert a STAR SJ.out.tab into BED records, one per splice-site boundary.

Every junction that survives the filters becomes *two* 1 bp BED records: the
donor boundary and the acceptor boundary. Both carry the same row number, so the
two ends of one junction can be paired again after any downstream interval
operation has shuffled them:

    <row>;s;<annotated|unannotated>    the intron_start boundary
    <row>;e;<annotated|unannotated>    the intron_end boundary

<row> is the junction's 1-based position among the *unfiltered* data rows, not
among the emitted ones, so it stays stable when --min-cov changes.

Coordinates: SJ.out.tab is 1-based inclusive and names the intron's first and
last base; BED is 0-based half-open. A boundary at 1-based position p is
therefore written as [p-1, p).

Columns of SJ.out.tab, which has no header:
  1 chromosome  2 intron_start  3 intron_end  4 strand  5 intron_motif
  6 annotation  7 n_unique_map  8 n_multi_map  9 max_splice_overhang

Output goes to stdout, tab separated and headerless, sorted by chromosome,
start, end and strand. Chromosome sorts lexicographically -- chr10 before chr2 --
which is what the BED consumers here already expect.

Usage:
  conv_sjouttab2bed.py [--min-cov=<n>] [--contigs=<list>] <sj-out-tab>
  conv_sjouttab2bed.py (-h | --help)

Options:
  --min-cov=<n>      Minimum uniquely-mapped read count for a junction to be
                     emitted [default: 3]
  --contigs=<list>   Comma-separated contigs to keep. Defaults to the canonical
                     human set chr1-chr22, chrX, chrY; pass an explicit list for
                     another assembly or to keep the scaffolds.
  -h --help          Show this message.
"""

import gzip
import sys

from docopt import docopt

# %%
# Constants

CANONICAL_CHROMOSOMES = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

N_COLUMNS = 9

# STAR encodes both of these as small integers in the file.
STRANDS = {0: ".", 1: "+", 2: "-"}
ANNOTATIONS = {0: "unannotated", 1: "annotated"}

BED_SCORE = 0


# %%
# Subs


def read_sj_records(path):
    """Yield (row_number, fields) for every non-comment data row.

    row_number is 1-based over the data rows only, so a leading comment block
    does not shift it.
    """
    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as handle:
        row_number = 0
        for line in handle:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) != N_COLUMNS:
                raise SystemExit(
                    f"conv_sjouttab2bed: {path} line {row_number + 1} has "
                    f"{len(fields)} fields, expected {N_COLUMNS}. "
                    "This does not look like a STAR SJ.out.tab."
                )

            row_number += 1
            yield row_number, fields


def build_bed(path, min_cov, contigs):
    """Return the BED rows as tuples, sorted the way the output is written."""
    kept = set(contigs)
    starts = []
    ends = []

    for row_number, fields in read_sj_records(path):
        chromosome = fields[0]
        if chromosome not in kept:
            continue

        if int(fields[6]) < min_cov:
            continue

        intron_start = int(fields[1])
        intron_end = int(fields[2])

        code = int(fields[3])
        if code not in STRANDS:
            raise SystemExit(
                f"conv_sjouttab2bed: unknown strand code {code} on data row "
                f"{row_number}; STAR writes 0, 1 or 2"
            )
        strand = STRANDS[code]

        code = int(fields[5])
        if code not in ANNOTATIONS:
            raise SystemExit(
                f"conv_sjouttab2bed: unknown annotation code {code} on data row "
                f"{row_number}; STAR writes 0 or 1"
            )
        annotation = ANNOTATIONS[code]

        starts.append((
            chromosome, intron_start - 1, intron_start,
            f"{row_number};s;{annotation}", BED_SCORE, strand,
        ))
        ends.append((
            chromosome, intron_end - 1, intron_end,
            f"{row_number};e;{annotation}", BED_SCORE, strand,
        ))

    # NOTE: every start record is placed before every end record, and only then
    #   sorted. The name field is deliberately absent from the sort key, so records
    #   that tie on position and strand keep this order -- Python's sort is stable.
    #
    #   This is not cosmetic. Two *different* junctions can put a start and an end
    #   on the same base and strand, and then the sort key alone does not decide
    #   which comes first. Building the list start-first reproduces the Julia
    #   implementation this replaces (`vcat(starts, ends)`) byte for byte;
    #   interleaving the pair per junction instead diverges on exactly those ties
    #   (4 lines out of 320524 on the reference SJ.out.tab).
    rows = starts + ends
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[5]))

    return rows


# %%
# Main


def main():
    args = docopt(__doc__)

    path = args["<sj-out-tab>"]
    min_cov = float(args["--min-cov"])

    if args["--contigs"]:
        contigs = [c.strip() for c in args["--contigs"].split(",") if c.strip()]
    else:
        contigs = CANONICAL_CHROMOSOMES

    rows = build_bed(path, min_cov, contigs)

    out = sys.stdout
    for row in rows:
        out.write("\t".join(str(field) for field in row))
        out.write("\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
