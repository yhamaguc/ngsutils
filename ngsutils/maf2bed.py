#!/usr/bin/env python3
"""Convert a MAF on stdin to BED6 on stdout.

One BED record per MAF data row. The four positional columns of a MAF -- Chromosome,
Start_Position, End_Position and Strand -- become the BED interval and its strand; every
remaining column is packed into the name field, so nothing is thrown away:

    chrom  start  end  "Hugo_Symbol";"Entrez_Gene_Id";...  0  strand

Values equal to the literal NA become empty, everything else is quoted, and the fields are
joined with ';' -- a semicolon rather than a comma because MAF free-text columns contain
commas.

Coordinates: MAF is 1-based inclusive, BED is 0-based half-open, so a variant at 1-based
[start, end] is written as [start - 1, end).

Padding: the slop option widens the interval on both sides, which is how you give a point
mutation some flanking context before intersecting it with anything -- a slop of 10 turns a
1 bp SNV into a 21 bp window. The lower edge is clamped at 0, since a negative coordinate is
not a legal BED start. The upper edge is NOT clamped, because that would need the contig
lengths and a MAF does not carry them, so a variant near the end of a contig can be padded
past it.

NOTE: that paragraph must not begin with a dash. docopt reads the whole docstring, and any
line whose first non-space character is '-' becomes an option definition -- this one used to
start with the option's own name and produced a junk '--' key. See the docopt entry under
Correctness traps in ~/.claude/CLAUDE.md.

Column positions are taken from the MAF specification (5 Chromosome, 6 Start_Position,
7 End_Position, 8 Strand), not from the header line, which is why a row with fewer than
8 columns is a hard error rather than something to skip quietly.

Usage:
  maf2bed [options]

Options:
  -s <INT>, --slop <INT>  Bases to pad the interval by on each side [default: 0]
  -h --help               Show this message.
"""

# NOTE: This script began as https://github.com/edawson/maf_to_bed and has been modified.
#   It now also carries the BED6 output (score and strand columns) that the Julia
#   conv_maf2bed.jl in the hs_tcga_pancancer_splicing repository produced; that script was
#   retired in favour of this one.

import sys

from docopt import docopt

# %%
# Constants

# 1-based MAF column numbers, per the MAF specification.
COLUMN_CHROMOSOME = 5
COLUMN_START = 6
COLUMN_END = 7
COLUMN_STRAND = 8

MIN_COLUMNS = COLUMN_STRAND

# Prefixes that mark a comment or a header rather than a data row.
SKIP_PREFIXES = ("#", "individual", "Hugo_Symbol")

# BED6 has no score to carry here, and 0 is the conventional placeholder.
BED_SCORE = 0

# BED accepts these three; the MAF specification says Strand is always '+', but a
# reheadered or hand-edited file is exactly the case worth catching.
VALID_STRANDS = ("+", "-", ".")

MISSING_VALUE = "NA"


# %%
# Subs


def die(message):
    print(f"maf2bed: error: {message}", file=sys.stderr)
    raise SystemExit(2)


def convert_line(line, line_number, slop):
    """Return one BED6 row as a list of fields, or None when the line is not data."""
    if line.startswith(SKIP_PREFIXES):
        return None

    line = line.strip()
    if not line:
        return None

    tokens = line.split("\t")
    if len(tokens) < MIN_COLUMNS:
        die(
            f"line {line_number} has {len(tokens)} columns, need at least {MIN_COLUMNS} "
            f"to reach Strand. This does not look like a MAF."
        )

    # NOTE: popping the same index repeatedly walks the four consecutive columns, and
    #   leaves `tokens` holding exactly the columns that are not positional.
    index = COLUMN_CHROMOSOME - 1
    chromosome = tokens.pop(index)
    start_raw = tokens.pop(index)
    end_raw = tokens.pop(index)
    strand = tokens.pop(index)

    try:
        start = int(start_raw)
        end = int(end_raw)
    except ValueError:
        die(
            f"line {line_number} has non-integer coordinates "
            f"(Start_Position={start_raw!r}, End_Position={end_raw!r})"
        )

    if strand not in VALID_STRANDS:
        die(
            f"line {line_number} has Strand={strand!r}, which BED does not accept; "
            f"expected one of {', '.join(VALID_STRANDS)}"
        )

    # MAF is 1-based inclusive, BED is 0-based half-open. Clamp the lower edge: --slop
    # larger than the start position would otherwise emit a negative coordinate.
    start_bed = max(0, start - slop - 1)
    end_bed = end + slop

    name = ";".join(
        f'"{value}"' if value != MISSING_VALUE else "" for value in tokens
    )

    return [chromosome, str(start_bed), str(end_bed), name, str(BED_SCORE), strand]


# %%
# Main


def main():
    options = docopt(__doc__)

    try:
        slop = int(options["--slop"])
    except ValueError:
        die(f"--slop must be an integer, got {options['--slop']!r}")

    if slop < 0:
        die(f"--slop must not be negative, got {slop}")

    out = sys.stdout
    for line_number, line in enumerate(sys.stdin, start=1):
        fields = convert_line(line, line_number, slop)
        if fields is None:
            continue
        out.write("\t".join(fields))
        out.write("\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
