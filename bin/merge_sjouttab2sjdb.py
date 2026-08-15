#!/usr/bin/env python3
"""Aggregate STAR SJ.out.tab files into one junction list for --sjdbFileChrStartEnd.

Emits one line per unique junction in the 4-column format the STAR manual
(2.7.11b) documents in section 2.2.4:

    Chr \t Start \t End \t Strand

Coordinates: STAR names a junction by its *intronic* bases, and both SJ.out.tab
(section 5.5) and --sjdbFileChrStartEnd (section 2.2.4) use the 1-based inclusive
first and last base of the intron. The two agree, so columns 1-3 pass through
unchanged. Never feed a BED interval here: BED is 0-based half-open, so its start
is one base short of what STAR expects, and STAR inserts the shifted junction
without complaining -- the intended junction is then missing from the index and
a bogus one takes its place.

Strand: SJ.out.tab column 4 encodes 0/1/2 while --sjdbFileChrStartEnd is
documented as +/-/., so this script converts.

Filtering: a junction is kept when its uniquely-mapped read count reaches
--min-cov in *at least one* input file. Counts are deliberately not pooled across
files -- n files contributing one read each is weaker evidence than one file
contributing n reads, and pooling would rank them alike. Keeping the per-file
maximum also means a junction real in a single sample survives, which matters
when sample-specific aberrant splicing is the object of study.

The threshold is applied to annotated and novel junctions alike. Dropping an
annotated junction here is harmless, because --sjdbGTFfile re-supplies it at
genome-generate time (section 9.3 step 2); applying the threshold uniformly keeps
this script from depending on the 1st-pass annotation being the same one used for
the build.

No motif filter is applied. STAR's --outSJfilter* defaults already require a
*novel* non-canonical junction to carry >= 3 unique reads and >= 30 nt overhang
before it reaches SJ.out.tab at all, and annotated junctions are exempt from
those filters and re-supplied by the GTF regardless. A motif filter here would
therefore act on an empty set in the novel case and on a redundant one otherwise.

Columns of SJ.out.tab, which has no header:
  1 chromosome  2 intron_start  3 intron_end  4 strand  5 intron_motif
  6 annotation  7 n_unique_map  8 n_multi_map  9 max_splice_overhang

Output goes to stdout, tab separated and headerless, sorted by chromosome, start
and end. Chromosome sorts lexicographically -- chr10 before chr2 -- matching
conv_sjouttab2bed.py. A summary of what was read, dropped and emitted goes to
stderr, so a run that silently kept far fewer junctions than expected is visible.

Usage:
  merge_sjouttab2sjdb.py [--min-cov=<n>] [--drop-contigs=<list>] <sj-out-tab>...
  merge_sjouttab2sjdb.py (-h | --help)

Options:
  --min-cov=<n>          Minimum uniquely-mapped read count a junction must reach
                         in at least one input file [default: 2]
  --drop-contigs=<list>  Comma-separated contigs to drop entirely. The default
                         drops the mitochondrion, whose junctions the STAR manual
                         names as likely false positives (section 9.3 step 2).
                         Pass an empty value to keep every contig.
                         [default: chrM]
  -h --help              Show this message.
"""

from __future__ import annotations

import gzip
import sys

from docopt import docopt

# %%
# Constants

N_COLUMNS = 9

# STAR writes the strand as a small integer; --sjdbFileChrStartEnd wants a symbol.
STRAND_SYMBOLS = {0: ".", 1: "+", 2: "-"}

UNDEFINED_STRAND_CODE = 0


# %%
# Subs


def read_sj_records(path):
    """Yield the split fields of every data row of one SJ.out.tab.

    This is the choke point for input shape: every row is checked to have
    N_COLUMNS fields here, and everything downstream trusts that.

    Parameters
    ----------
    path : str
        Path to an SJ.out.tab, optionally gzipped.

    Yields
    ------
    list[str]
        The nine tab-separated fields of one data row.
    """
    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) != N_COLUMNS:
                raise SystemExit(
                    f"merge_sjouttab2sjdb: {path} line {line_number} has "
                    f"{len(fields)} fields, expected {N_COLUMNS}. "
                    "This does not look like a STAR SJ.out.tab."
                )

            yield fields


def aggregate_junctions(paths, drop_contigs):
    """Collapse every input file into per-junction read counts and strand codes.

    Parameters
    ----------
    paths : list[str]
        SJ.out.tab paths, optionally gzipped.
    drop_contigs : set[str]
        Contigs to discard on sight.

    Returns
    -------
    max_unique : dict[tuple[str, int, int], int]
        Highest uniquely-mapped read count seen for the junction in any one file.
    strand_codes : dict[tuple[str, int, int], set[int]]
        Every strand code the junction was seen with.
    counts : dict[str, int]
        Row tallies for the stderr summary.
    """
    max_unique = {}
    strand_codes = {}
    counts = {"rows": 0, "rows_dropped_contig": 0}

    for path in paths:
        for fields in read_sj_records(path):
            counts["rows"] += 1

            chromosome = fields[0]
            if chromosome in drop_contigs:
                counts["rows_dropped_contig"] += 1
                continue

            key = (chromosome, int(fields[1]), int(fields[2]))

            code = int(fields[3])
            if code not in STRAND_SYMBOLS:
                raise SystemExit(
                    f"merge_sjouttab2sjdb: unknown strand code {code} for "
                    f"{chromosome}:{fields[1]}-{fields[2]} in {path}; "
                    "STAR writes 0, 1 or 2"
                )

            n_unique = int(fields[6])

            if key in max_unique:
                if n_unique > max_unique[key]:
                    max_unique[key] = n_unique
                strand_codes[key].add(code)
            else:
                max_unique[key] = n_unique
                strand_codes[key] = {code}

    return max_unique, strand_codes, counts


def resolve_strand(key, codes):
    """Collapse the strand codes a junction was seen with into one symbol.

    An undefined code (0) loses to a defined one, since STAR leaves the strand
    undefined only when it could not call it. Seeing both + and - for the same
    intron is a contradiction rather than something to average over, so it stops
    the run.

    Parameters
    ----------
    key : tuple[str, int, int]
        Chromosome, 1-based intron start, 1-based intron end.
    codes : set[int]
        Strand codes observed across the input files.

    Returns
    -------
    str
        One of "+", "-" or ".".
    """
    defined = codes - {UNDEFINED_STRAND_CODE}

    if len(defined) > 1:
        chromosome, start, end = key
        raise SystemExit(
            f"merge_sjouttab2sjdb: junction {chromosome}:{start}-{end} is on "
            "the + strand in one input file and the - strand in another. The "
            "inputs do not describe one assembly. Refusing to guess."
        )

    if not defined:
        return STRAND_SYMBOLS[UNDEFINED_STRAND_CODE]

    return STRAND_SYMBOLS[defined.pop()]


def build_sjdb_rows(max_unique, strand_codes, min_cov):
    """Return the surviving junctions as sorted 4-tuples, plus tallies.

    Parameters
    ----------
    max_unique : dict[tuple[str, int, int], int]
        Per-junction maximum uniquely-mapped read count.
    strand_codes : dict[tuple[str, int, int], set[int]]
        Per-junction observed strand codes.
    min_cov : int
        Minimum read count a junction must reach in at least one input file.

    Returns
    -------
    rows : list[tuple[str, int, int, str]]
        Chromosome, 1-based intron start, 1-based intron end, strand symbol.
    counts : dict[str, int]
        Junction tallies for the stderr summary.
    """
    counts = {
        "junctions": len(max_unique),
        "junctions_below_min_cov": 0,
        "junctions_undefined_strand": 0,
    }

    rows = []
    for key, n_unique in max_unique.items():
        if n_unique < min_cov:
            counts["junctions_below_min_cov"] += 1
            continue

        strand = resolve_strand(key, strand_codes[key])
        if strand == STRAND_SYMBOLS[UNDEFINED_STRAND_CODE]:
            counts["junctions_undefined_strand"] += 1

        rows.append((key[0], key[1], key[2], strand))

    rows.sort(key=lambda row: (row[0], row[1], row[2]))

    return rows, counts


def report(handle, n_files, min_cov, drop_contigs, row_counts, junction_counts, n_emitted):
    """Write the run summary, naming how much the threshold removed."""
    kept_fraction = 100 * n_emitted / junction_counts["junctions"] if junction_counts["junctions"] else 0

    handle.write(f"merge_sjouttab2sjdb: {n_files} input file(s)\n")
    handle.write(f"  rows read                     : {row_counts['rows']}\n")
    handle.write(
        f"  rows dropped by --drop-contigs : {row_counts['rows_dropped_contig']}"
        f"  ({','.join(sorted(drop_contigs)) if drop_contigs else 'none'})\n"
    )
    handle.write(f"  unique junctions              : {junction_counts['junctions']}\n")
    handle.write(
        f"  dropped below --min-cov={min_cov}      : "
        f"{junction_counts['junctions_below_min_cov']}\n"
    )
    handle.write(f"  emitted                       : {n_emitted}  ({kept_fraction:.1f}% of unique)\n")
    handle.write(
        f"  emitted with undefined strand : {junction_counts['junctions_undefined_strand']}\n"
    )


# %%
# Main


def main():
    args = docopt(__doc__)

    paths = args["<sj-out-tab>"]
    min_cov = int(args["--min-cov"])
    drop_contigs = {c.strip() for c in args["--drop-contigs"].split(",") if c.strip()}

    max_unique, strand_codes, row_counts = aggregate_junctions(paths, drop_contigs)

    if not max_unique:
        raise SystemExit(
            "merge_sjouttab2sjdb: no junctions survived reading the inputs. "
            f"{row_counts['rows']} row(s) were read and "
            f"{row_counts['rows_dropped_contig']} dropped by --drop-contigs."
        )

    rows, junction_counts = build_sjdb_rows(max_unique, strand_codes, min_cov)

    out = sys.stdout
    for row in rows:
        out.write("\t".join(str(field) for field in row))
        out.write("\n")

    report(
        sys.stderr, len(paths), min_cov, drop_contigs, row_counts, junction_counts, len(rows)
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
