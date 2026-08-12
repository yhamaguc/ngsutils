#!/usr/bin/env python3
"""write_coverage_track — per-base coverage as bedGraph, for a region or an Ensembl ID.

One job: resolve a set of target intervals, then write one bedGraph per BAM over exactly those
intervals. Nothing is summarised and nothing else is written.

Two ways to say what to measure:

  1. by coordinate      --region chr17:7673700-7673837
  2. by ID (+ flanks)   --id ENST00000269305 --lower 200 --upper 200

The ID is resolved against a GENCODE/Ensembl GTF, its exons become the query intervals, and
`samtools depth -a` supplies the depth.

    ngsutils write_coverage_track --output-dir tracks --gtf gencode.v50.gtf.gz --id ENST00000269305 s.bam
    ngsutils write_coverage_track --output-dir tracks --gtf ... --id TP53 --mode merged_exon --flank 500 *.bam
    ngsutils write_coverage_track --output-dir tracks --region chr17:7673700-7673837 a.bam b.bam

BAMs are positional and repeatable, as `samtools depth [options] in.bam [in.bam ...]` takes them,
so a file list composes directly:

    find data -name '*.bam' | sort | xargs ngsutils write_coverage_track --output-dir tracks \
      --gtf ... --id ENST00000269305 --label-from parent

Why the locus and not the exons. **`span` is the default since 2026-08-12**; it was `auto`, which
is exonic. A track that stops at the exons cannot show a retained intron, and intronic coverage is
a signal rather than noise: `samtools depth` does not count CIGAR `N`, so a read that spliced
across an intron contributes 0 there and only a read aligned *through* the intron leaves depth.
Intronic depth in a `span` track is therefore the intron-retaining population, and the RI ratio is
that depth over the flanking exons' — a query over this track, not a mode of this program.

The cost is real and was the old default's reason: TP53 ENST00000269305 spans 19,070 bp of chr17
but holds only 2,512 bp of exon, so a `span` track is about 7.6x the positions and mostly intron.
Use `merged_exon` or `exon` when that is not wanted. What depth alone cannot separate: a retained
intron from nascent pre-mRNA, from an overlapping feature on either strand (this program does not
split by strand), or from the primary alignment of a multi-mapping read in a repeat — `--min-mq`
is the only lever here for the last of those.

Sums and means are deliberately absent. They are a length-weighted sum over this track, so the
track is the primitive and a summary is a query over it:

    bedtools intersect -a tracks/s1.bedgraph -b exons.bed -wb \
    | awk -F'\t' '{s[$8] += ($3-$2)*$4} END {for (k in s) print k, s[k]}'

  Measured against `samtools bedcov -j` on a `10M20N10M` read: identical. Note that
  `bedtools map -c 4 -o sum` is **wrong** here and returns 2 where the answer is 20 — it adds the
  run-length-encoded values without weighting them by interval length. Clip first, then weight.

Zero is written, never omitted. A missing interval in a bedGraph means *no data*, so a zero-depth
base left out would read as unmeasured. `-a` forces every queried position to be reported and every
one is written, so an absent interval means one thing only: outside the requested regions.

Intron skips are excluded, and this was measured rather than assumed: for a `10M20N10M` read
`samtools depth` reports 0 inside the N. There is no flag to change that (`-J` covers deletions
only), so an RNA-seq read spanning an intron never inflates an exon it did not touch.

Labels name a sample once, and the name is both the bedGraph filename and the track name.
`stem` is the default: `s1/aligned.bam` becomes `aligned`. That collapses under the layout `find`
walks, where every sample directory holds a file of the same name, so `parent-stem` prepends the
directory (`s1_aligned`) the way MultiQC's `-d` does, and `parent` uses the directory alone.
Duplicate labels are a hard error: the second file would silently overwrite the first.

Python 3.8+ and `docopt`, both of which `ngsutils` already requires. Also needs `samtools` on
PATH, which `pyproject.toml` cannot express — it is an external program, not a Python package.

Modes:
  span         the locus, introns included                    (THE DEFAULT)
  merged_exon  union of every isoform's exons, introns cut out
  exon         one interval per exon
  auto         exonic, chosen by ID kind: merged_exon for a gene or gene name, exon otherwise

  `merged_exon` was called `merged` until 2026-08-12, which did not say what was merged.

Flanks:
  Lower and upper are plus-strand directions, NOT 5'/3'. On a minus-strand feature the upper
  flank is the 5' side. For an ID they extend the target's outer boundary; for a region the
  region itself is extended. Either way they are queried like any other base.
  (No line of this docstring may begin with a dash except the Options entries below — see the
  note above `options()`.)

NOTE: invoked as `ngsutils write_coverage_track`, but the Usage patterns below name ONE word.
docopt reads the first word of a pattern as the program name and every later word as a
required argument, so `ngsutils write_coverage_track ...` there would demand a literal
`write_coverage_track` argument that `cli.py` has already consumed. `--help` still passes,
because it short-circuits before the pattern is matched — so this fails only in real use.
This note is deliberately outside the Usage block: that block is parsed as patterns, and
prose in it is a DocoptLanguageError.

Usage:
  write_coverage_track --output-dir=<dir> [--id=<id>...] [--region=<spec>...] [--label=<name>...]
                          [options] <bam>...
  write_coverage_track (-h | --help | --version)

Options:
  -o <dir>, --output-dir=<dir>   Directory for the bedGraph files. Created if absent.
  -a <file>, --gtf=<file>        GENCODE/Ensembl GTF, .gz accepted. Required with --id.
  -t <id>, --id=<id>             ENSG/ENST/ENSE ID, or a gene name. Repeatable. Version optional.
  -r <spec>, --region=<spec>     1-based inclusive CONTIG:START-END, optionally NAME=CONTIG:START-END.
                          Repeatable.
  -M <mode>, --mode=<mode>       auto | exon | merged_exon | span [default: span]
  -L <bp>, --lower=<bp>          Flank toward lower coordinates [default: 0]
  -U <bp>, --upper=<bp>          Flank toward higher coordinates [default: 0]
  -F <bp>, --flank=<bp>          Shorthand for --lower BP --upper BP.
  -l <name>, --label=<name>      Name for one BAM, repeatable; as many as <bam>, in the same order.
  -f <source>, --label-from=<source>   Derive labels instead: basename | stem | parent | parent-stem
                          [default: stem]
  -T, --no-track-line            Omit the `track type=bedGraph` line, for a bedtools consumer.
  -Q <int>, --min-mq=<int>       Mapping-quality threshold (samtools -Q).
  -G <str>, --excl-flags=<str>   Flags to exclude, e.g. UNMAP,SECONDARY,QCFAIL,DUP (samtools -G).
  -@ <int>, --threads=<int>      Additional decompression threads (samtools depth -@).
  -q, --quiet                    Silence stderr progress.
  -h --help               Show this message.
  -V, --version                  Show the version.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Tuple

from docopt import docopt

# One interval. `start` is 0-based and `end` exclusive, as BED. `part` says what the interval is —
# an exon, a merged_exon block, a span, or a flank — and is reported rather than written, so a flank is
# never mistaken for exonic sequence in the run log.
Interval = Tuple[str, int, int, str, str, str]  # contig, start, end, name, strand, part

VERSION = "write_coverage_track 2.0"

MODES = ("auto", "exon", "merged_exon", "span")

LABEL_SOURCES = ("basename", "stem", "parent", "parent-stem")
LABEL_FROM_DEFAULT = "stem"
BAM_SUFFIXES = (".bam", ".cram", ".sam")
# A label becomes a filename, so these would either break the path or produce a file nobody asked
# for. Checked rather than sanitised: silently rewriting a name the caller chose is worse.
LABEL_FORBIDDEN = ("/", "\\", "\0")


def die(msg: str, code: int = 2) -> None:
    print(f"write_coverage_track: error: {msg}", file=sys.stderr)
    raise SystemExit(code)


def note(msg: str) -> None:
    print(f"  {msg}", file=sys.stderr)


# %%
# Targets

REGION_RE = re.compile(
    r"^(?:(?P<name>[^=]+)=)?(?P<contig>[^:]+):(?P<start>[\d,]+)-(?P<end>[\d,]+)$"
)
ENSEMBL_RE = re.compile(r"^(?P<prefix>ENS[A-Z]*)(?P<kind>[GTE])(?P<num>\d+)(?:\.\d+)?$")


def base_id(value: str) -> str:
    """Ensembl IDs carry a version suffix in GENCODE; queries usually do not."""
    return value.split(".", 1)[0]


def parse_region(spec: str) -> Interval:
    """`chr17:7673700-7673837`, optionally `name=chr17:...`. 1-based inclusive, as samtools."""
    m = REGION_RE.match(spec.strip())
    if not m:
        die(f"--region {spec!r} is not CONTIG:START-END (1-based inclusive)")
    start1 = int(m.group("start").replace(",", ""))
    end = int(m.group("end").replace(",", ""))
    if start1 < 1 or end < start1:
        die(f"--region {spec!r}: START must be >= 1 and END >= START")
    name = m.group("name") or f'{m.group("contig")}:{start1}-{end}'
    # BED is 0-based half-open; GTF and samtools regions are 1-based inclusive. Getting this
    # wrong shifts every interval by one base, silently.
    return (m.group("contig"), start1 - 1, end, name, ".", "region")


def id_kind(value: str) -> str:
    """`gene` / `transcript` / `exon` from the ID itself, or `gene_name` when it is not an ID."""
    m = ENSEMBL_RE.match(value)
    if not m:
        return "gene_name"
    return {"G": "gene", "T": "transcript", "E": "exon"}[m.group("kind")]


# %%
# GTF

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')
EXON_NUMBER_RE = re.compile(r"exon_number\s+\"?(\d+)")


def attributes(field: str) -> Dict[str, str]:
    """GTF column 9. Only quoted values; `exon_number` is unquoted and read separately."""
    return dict(ATTR_RE.findall(field))


def open_maybe_gzip(path: str):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "rt")


def scan_gtf(gtf: str, wanted: Dict[str, str], quiet: bool) -> Dict[str, List[Interval]]:
    """One pass over the GTF, collecting exons for every requested ID at once.

    The single pass is the point: the GENCODE v50 GTF holds 5,087,789 exon records and takes about
    a minute to read, and scanning it once for fifty IDs costs the same as scanning it once for
    one. It is also why every BAM is handled in one invocation rather than one run per sample.
    """
    if not os.path.exists(gtf):
        die(f"--gtf {gtf}: no such file")
    by_query: Dict[str, List[Interval]] = {q: [] for q in wanted}
    keys = {
        "gene": "gene_id",
        "transcript": "transcript_id",
        "exon": "exon_id",
        "gene_name": "gene_name",
    }
    # query -> (attribute key, value to compare). Compared on the base ID, so a query without a
    # version suffix still matches GENCODE's versioned one.
    probes = [
        (q, keys[k], base_id(q) if k != "gene_name" else q) for q, k in wanted.items()
    ]

    if not quiet:
        note(f"scanning {gtf}")
    n_exon = 0
    with open_maybe_gzip(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            n_exon += 1
            attrs = attributes(f[8])
            for query, key, target in probes:
                got = attrs.get(key)
                if got is None:
                    continue
                if (got if key == "gene_name" else base_id(got)) != target:
                    continue
                m = EXON_NUMBER_RE.search(f[8])
                num = m.group(1) if m else ""
                tid = attrs.get("transcript_id", "")
                eid = attrs.get("exon_id", "")
                if key == "exon_id":
                    # One physical exon, listed once per transcript that contains it. Naming it
                    # after the exon rather than after each transcript is what lets the duplicate
                    # rows collapse to the single interval they describe.
                    name = eid or query
                else:
                    name = f"{tid}.e{num}" if tid and num else (eid or query)
                by_query[query].append((f[0], int(f[3]) - 1, int(f[4]), name, f[6], "exon"))
    if not quiet:
        note(f"{n_exon} exon records read")
    return by_query


# %%
# Interval shaping

def dedupe(intervals: Sequence[Interval], label: str, quiet: bool) -> List[Interval]:
    """Drop exactly-repeated intervals within one target.

    GENCODE lists a shared exon once per transcript that contains it, so an `ENSE` query returns
    the same interval several times.
    """
    seen = set()
    out: List[Interval] = []
    for iv in intervals:
        key = (iv[0], iv[1], iv[2], iv[3])
        if key in seen:
            continue
        seen.add(key)
        out.append(iv)
    dropped = len(intervals) - len(out)
    if dropped and not quiet:
        note(f"{label}: {dropped} duplicate interval(s) collapsed")
    return out


def merge(intervals: Sequence[Interval], label: str) -> List[Interval]:
    """Union of overlapping or touching intervals — the collapsed exon model.

    Needed for a gene: its transcripts share constitutive exons, so unmerged intervals overlap.
    """
    out: List[Interval] = []
    for contig, start, end, _, strand, _ in sorted(intervals, key=lambda i: (i[0], i[1], i[2])):
        if out and out[-1][0] == contig and start <= out[-1][2]:
            prev = out[-1]
            out[-1] = (prev[0], prev[1], max(prev[2], end), prev[3], prev[4], "merged_exon")
        else:
            out.append((contig, start, end, f"{label}.m{len(out) + 1}", strand, "merged_exon"))
    return out


def span(intervals: Sequence[Interval], label: str) -> List[Interval]:
    """Min..max per contig. Includes introns; here only because it is sometimes what is wanted."""
    per: Dict[str, List[int]] = {}
    for contig, start, end, _, _, _ in intervals:
        if contig in per:
            lo, hi = per[contig]
            per[contig] = [min(lo, start), max(hi, end)]
        else:
            per[contig] = [start, end]
    strand = intervals[0][4] if intervals else "."
    return [(c, lo, hi, f"{label}.span", strand, "span") for c, (lo, hi) in sorted(per.items())]


def flanks(
    core: Sequence[Interval], label: str, lower: int, upper: int, lengths: Dict[str, int]
) -> List[Interval]:
    """Flanking intervals outside the target's outer boundary, one per side per contig.

    Separate intervals rather than extended exons, so the run log can say how many bases are
    flanking rather than folding them into the exonic total.

    `lower` extends toward lower coordinates and `upper` toward higher ones. Those are plus-strand
    directions, not 5'/3' — on a minus-strand feature `upper` is the 5' side.
    """
    out: List[Interval] = []
    if lower <= 0 and upper <= 0:
        return out
    per: Dict[str, List[int]] = {}
    strand = core[0][4] if core else "."
    for contig, start, end, _, _, _ in core:
        if contig in per:
            lo, hi = per[contig]
            per[contig] = [min(lo, start), max(hi, end)]
        else:
            per[contig] = [start, end]

    for contig, (lo, hi) in sorted(per.items()):
        if lower > 0:
            new_lo = max(0, lo - lower)
            if new_lo < lo:
                if new_lo != lo - lower:
                    note(f"{label}: lower flank clipped at the start of {contig}")
                out.append((contig, new_lo, lo, f"{label}.flank_lower", strand, "flank_lower"))
            else:
                note(f"{label}: no room for a lower flank on {contig}")
        if upper > 0:
            limit = lengths.get(contig)
            new_hi = hi + upper if limit is None else min(limit, hi + upper)
            if new_hi > hi:
                if new_hi != hi + upper:
                    note(f"{label}: upper flank clipped at the end of {contig} ({limit} bp)")
                out.append((contig, hi, new_hi, f"{label}.flank_upper", strand, "flank_upper"))
            else:
                note(f"{label}: no room for an upper flank on {contig}")
    return out


def extend(iv: Interval, lower: int, upper: int, lengths: Dict[str, int]) -> Interval:
    """Grow a single interval in place — what flanking means for an explicit `--region`."""
    contig, start, end, name, strand, part = iv
    new_start = max(0, start - lower)
    limit = lengths.get(contig)
    new_end = end + upper if limit is None else min(limit, end + upper)
    if lower and new_start != start - lower:
        note(f"{name}: lower flank clipped at the start of {contig}")
    if upper and new_end != end + upper:
        note(f"{name}: upper flank clipped at the end of {contig} ({limit} bp)")
    suffix = f"[-{lower},+{upper}]" if (lower or upper) else ""
    return (contig, new_start, new_end, name + suffix, strand, part)


def union(intervals: Sequence[Interval]) -> List[Tuple[str, int, int]]:
    """Overlapping intervals collapsed to disjoint `(contig, start, end)`, in coordinate order.

    A position's depth does not depend on which exon asked for it, but `--mode exon` on a gene
    yields one overlapping interval per isoform, and `samtools depth -b` would then report that
    position once per interval. Two bedGraph intervals covering the same base is malformed output,
    so the query is unioned first.

    Coordinate order because a BED is read that way and GENCODE lists a minus-strand transcript's
    exons 5'->3', which is descending. **It does not decide the output order**: an earlier version
    of this docstring said the order mattered because "samtools depth seeks the index per region",
    which was wrong twice over — `-b` alone does not seek at all, and the rows come back in BAM
    header order regardless. `seek_regions` is what puts the `-r` calls in that order.
    """
    out: List[Tuple[str, int, int]] = []
    for contig, start, end, _, _, _ in sorted(intervals, key=lambda i: (i[0], i[1], i[2])):
        if out and out[-1][0] == contig and start <= out[-1][2]:
            prev = out[-1]
            out[-1] = (prev[0], prev[1], max(prev[2], end))
        else:
            out.append((contig, start, end))
    return out


# %%
# Labels

def drop_bam_suffix(name: str) -> str:
    for suffix in BAM_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def resolve_labels(bams: Sequence[str], explicit: Sequence[str], source: str) -> List[str]:
    """One label per BAM: the bedGraph filename and the track name.

    Duplicates are an error, not a disambiguation. The second file would overwrite the first, so
    twenty BAMs would become fewer files with every number in them correct — and guessing a suffix
    would put a name in the output that appears in no argument.
    """
    if explicit:
        if len(explicit) != len(bams):
            die(
                f"{len(explicit)} --label for {len(bams)} BAM(s): give one per BAM in the same "
                f"order, or use --label-from to derive them"
            )
        labels = list(explicit)
    else:
        labels = []
        for bam in bams:
            base = os.path.basename(bam)
            # abspath first: the parent of `./x.bam` is `.`, which names nothing.
            parent = os.path.basename(os.path.dirname(os.path.abspath(bam)))
            if source == "basename":
                labels.append(base)
            elif source == "stem":
                labels.append(drop_bam_suffix(base))
            elif source == "parent":
                labels.append(parent)
            else:  # parent-stem
                labels.append(f"{parent}_{drop_bam_suffix(base)}")

    for label, bam in zip(labels, bams):
        if not label.strip():
            die(f"empty label for {bam}")
        for bad in LABEL_FORBIDDEN:
            if bad in label:
                die(f"label {label!r} for {bam} contains {bad!r}; it is used as a filename")

    seen: Dict[str, List[str]] = {}
    for label, bam in zip(labels, bams):
        seen.setdefault(label, []).append(bam)
    clashes = {k: v for k, v in seen.items() if len(v) > 1}
    if clashes:
        detail = "; ".join(f"{k!r} <- " + ", ".join(v) for k, v in sorted(clashes.items()))
        # The hint has to name the knob the caller actually used, or it sends them to a flag that
        # cannot affect the labels they supplied by hand.
        if explicit:
            hint = "these came from --label, so give distinct values"
        elif source == "stem":
            hint = "--label-from parent-stem distinguishes them"
        else:
            hint = f"--label-from {source} does not distinguish these; give --label explicitly"
        die(
            f"duplicate label(s): {detail}. Each label names one bedGraph file, so a repeat would "
            f"overwrite it. {hint}"
        )
    return labels


# %%
# samtools

def samtools_or_die() -> str:
    exe = shutil.which("samtools")
    if not exe:
        die("samtools not found on PATH")
    return exe  # type: ignore[return-value]


def contig_lengths(bams: Sequence[str]) -> Dict[str, int]:
    """@SQ lengths from the first BAM header, used to clip an upper flank.

    The GTF carries no contig lengths, so without this a flank past the end of a contig becomes a
    samtools error at the very end of the run.
    """
    exe = shutil.which("samtools")
    if not exe or not bams:
        return {}
    p = subprocess.run([exe, "view", "-H", bams[0]], capture_output=True, text=True)
    if p.returncode != 0:
        return {}
    out: Dict[str, int] = {}
    for line in p.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        name = length = None
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                length = int(field[3:])
        if name and length:
            out[name] = length
    return out


def check_bam(path: str) -> None:
    if not os.path.exists(path):
        die(f"{path}: no such file")
    stem = os.path.splitext(path)[0]
    for candidate in (path, stem):
        for suffix in (".bai", ".csi", ".crai"):
            if os.path.exists(candidate + suffix):
                return
    # Refused rather than worked around. Dropping `-r` when the index is missing would still be
    # correct and would be silently ~300x slower on a whole-BAM scan, which is the kind of quiet
    # degradation this program is not allowed to choose for the caller.
    die(
        f"{path} has no index — run `samtools index {path}` first. depth is given -r so that it "
        f"seeks to the target instead of reading the whole file, and -r requires the index"
    )


def seek_regions(
    blocks: Sequence[Tuple[str, int, int]], contig_order: Sequence[str]
) -> List[Tuple[str, str]]:
    """One `-r` region per contig, in BAM header order. `(contig, "contig:start-end")`, 1-based.

    **`samtools depth -b BED` does not use the index.** It reads every alignment in the file, and
    the `-b` list is applied to what streams past. Adding `-r` makes it seek. Verified here, on
    samtools 1.19.2, by what each form requires rather than by timing: `-b` alone runs to
    completion on a BAM with **no index at all**, while adding `-r` fails with `cannot load index`.
    A form that does not need the index cannot be using it.

    The speed this buys was measured elsewhere and is **not** measured in this repository, which
    holds no BAM: 448 BAMs totalling 6.77 TB took 16 h 28 m to reach 64 % under `-b` alone and
    34 s in full with `-r` (reported 2026-08-12). Treat the ratio as that setup's, not as a
    property of samtools.

    **One region per contig, never one span across them.** `-r` takes a single region, so a
    min-max span over a multi-contig BED keeps the first contig and drops the rest — measured on
    a two-contig fixture as 50 rows becoming 30, with **exit 0 and nothing on stderr**. That is a
    track that looks complete and is not, which is the failure this program exists to refuse.

    **Ordered by the BAM header, not by contig name.** `-b` alone emits in header order; a loop
    driven by the union's name sort would emit `chr10` before `chr2` where samtools emits `chr2`
    first, and a bedGraph whose intervals are reordered is a different file. Both orders were
    checked against a fixture whose header order and name order disagree.
    """
    extent: Dict[str, Tuple[int, int]] = {}
    for contig, start, end in blocks:
        lo, hi = extent.get(contig, (start, end))
        extent[contig] = (min(lo, start), max(hi, end))

    unknown = [c for c in extent if c not in contig_order]
    if unknown:
        die(
            f"contig(s) {', '.join(sorted(unknown))} are not in the BAM header, so no region can "
            f"be built for them. The target and the alignments disagree about contig naming "
            f"(`chr1` against `1`), or the BAM is not the one this target belongs to"
        )
    # BED is 0-based half-open; a samtools region is 1-based inclusive.
    return [
        (contig, f"{contig}:{extent[contig][0] + 1}-{extent[contig][1]}")
        for contig in contig_order
        if contig in extent
    ]


def run_depth(
    bed: str,
    bams: Sequence[str],
    args,
    blocks: Sequence[Tuple[str, int, int]],
    contig_order: Sequence[str],
) -> List[List[str]]:
    """`samtools depth -a` over the regions: one row per position, one depth column per BAM.

    `-a` is not optional. Without it a position no read overlaps is omitted, and an omitted
    position in a bedGraph means "not measured" — measured: inside a read's span the ref-skipped
    intron is reported as 0 either way, but past the read's end only `-a` reports anything.

    `-J` is deliberately not passed, so deletions are excluded. There is no flag for CIGAR N and
    none is needed: it is excluded already, which is what makes this track comparable with
    `samtools bedcov -j`.

    **The BAMs are not processed one after another.** One `depth` call takes all of them and reports
    one column per file at each position, so they are read together in a single pass over the
    queried positions. `--threads` adds decompression threads to that pass; it does not fan the
    files out across cores, and the pileup loop itself stays single-threaded. Whether it helps is a
    property of the data and was not measured here — this repository holds no BAM.

    **One call per contig**, each carrying `-r` so that samtools seeks rather than scanning the
    whole file — see [`seek_regions`], which is where that reasoning and its measurements live.
    `-b` is still passed to every call: `-r` selects what is *read*, `-b` selects what is
    *reported*, and the output is byte-identical to the single `-b`-only call this replaced.
    """
    rows: List[List[str]] = []
    for _, region in seek_regions(blocks, contig_order):
        cmd = [samtools_or_die(), "depth", "-a", "-b", bed, "-r", region]
        if args.min_mq is not None:
            cmd += ["-Q", str(args.min_mq)]
        if args.excl_flags is not None:
            # NOTE: depth's -G *adds* to a default filter-out list (UNMAP,SECONDARY,QCFAIL,DUP).
            cmd += ["-G", args.excl_flags]
        if args.threads is not None:
            cmd += ["-@", str(args.threads)]
        cmd += list(bams)
        if not args.quiet:
            note("running: " + " ".join(cmd))
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            sys.stderr.write(p.stderr)
            die(f"samtools depth exited {p.returncode}", code=1)
        if p.stderr.strip() and not args.quiet:
            sys.stderr.write(p.stderr)
        rows += [line.split("\t") for line in p.stdout.splitlines() if line.strip()]
    return rows


def write_bedgraph(
    directory: str, rows: Sequence[Sequence[str]], labels: Sequence[str], track_line: bool
) -> List[Tuple[str, int]]:
    """One bedGraph per BAM. Returns `(path, interval count)` per BAM, in the order given.

    bedGraph is 0-based half-open and single-track, so it is one file per sample rather than one
    file with a column each. Consecutive positions of equal depth are merged into one interval —
    that is what the format is for — but a run is broken by a coordinate gap as well as by a change
    in value, because merging across an intron would claim coverage for bases never queried.

    The file is named after the label, which `resolve_labels` has already checked is unique and
    usable as a filename, so nothing here can overwrite a file it just wrote.
    """
    os.makedirs(directory, exist_ok=True)
    written: List[Tuple[str, int]] = []
    for i, label in enumerate(labels):
        path = os.path.join(directory, f"{label}.bedgraph")
        n = 0
        with open(path, "w") as fh:
            if track_line:
                fh.write(f'track type=bedGraph name="{label}" visibility=full\n')
            run: Optional[List] = None  # [contig, start0, end0, depth]
            for r in rows:
                contig, pos1, depth = r[0], int(r[1]), int(r[2 + i])
                start0 = pos1 - 1
                if run is not None and run[0] == contig and run[2] == start0 and run[3] == depth:
                    run[2] = pos1
                    continue
                if run is not None:
                    fh.write(f"{run[0]}\t{run[1]}\t{run[2]}\t{run[3]}\n")
                    n += 1
                run = [contig, start0, pos1, depth]
            if run is not None:
                fh.write(f"{run[0]}\t{run[1]}\t{run[2]}\t{run[3]}\n")
                n += 1
        written.append((path, n))
    return written


# %%
# Main

def as_int(raw: Optional[str], flag: str) -> Optional[int]:
    """`raw` as an integer, or `None` when the flag was not given.

    docopt hands everything back as a string and validates nothing about the value, so every
    numeric flag passes through here — one place to be wrong rather than five.
    """
    if raw is None:
        return None
    try:
        return int(raw)
    except ValueError:
        die(f"{flag} expects an integer, got {raw!r}")
        return None  # unreachable


def options(argv: Optional[Sequence[str]] = None) -> SimpleNamespace:
    """The command line as the attribute names the rest of the script uses."""
    args = docopt(__doc__, argv=argv, version=VERSION)

    # docopt 0.6.2 parses **every** line of the docstring that begins with a dash as an option
    # definition — `parse_defaults` splits the whole doc, not an `Options:` section. So a prose
    # line starting with `--upper` redefines `--upper` and silently drops its `[default: 0]`,
    # leaving `None` where an integer is expected. This check turns that into a named error
    # instead of a TypeError several functions later; it fired once already.
    for flag in ("--mode", "--lower", "--upper", "--label-from"):
        if args[flag] is None:
            die(
                f"internal: {flag} lost its default — some line of this script's docstring "
                "begins with a dash and has redefined it"
            )

    mode = args["--mode"]
    if mode not in MODES:
        die(f"--mode {mode!r} is not one of: {', '.join(MODES)}")
    if args["--label-from"] not in LABEL_SOURCES:
        die(
            f"--label-from {args['--label-from']!r} is not one of: "
            f"{', '.join(LABEL_SOURCES)}"
        )
    if args["--label"] and args["--label-from"] != LABEL_FROM_DEFAULT:
        die("--label and --label-from set the same thing; give one of them")

    return SimpleNamespace(
        outdir=args["--output-dir"],
        gtf=args["--gtf"],
        id=args["--id"],
        region=args["--region"],
        bam=args["<bam>"],
        mode=mode,
        lower=as_int(args["--lower"], "--lower"),
        upper=as_int(args["--upper"], "--upper"),
        flank=as_int(args["--flank"], "--flank"),
        label=args["--label"],
        label_from=args["--label-from"],
        no_track_line=args["--no-track-line"],
        min_mq=as_int(args["--min-mq"], "--min-mq"),
        excl_flags=args["--excl-flags"],
        threads=as_int(args["--threads"], "--threads"),
        quiet=args["--quiet"],
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = options(argv)
    if args.flank is not None:
        args.lower = args.upper = args.flank
    if args.lower < 0 or args.upper < 0:
        die("--lower and --upper must be >= 0")
    if not args.id and not args.region:
        die("give at least one --id or --region")
    if args.id and not args.gtf:
        die("--id needs --gtf to resolve coordinates from")

    # Before any scan: a duplicate label is a naming error in the arguments, and there is no reason
    # to spend a minute of GTF reading discovering it.
    labels = resolve_labels(args.bam, args.label, args.label_from)
    for bam in args.bam:
        check_bam(bam)

    lengths = contig_lengths(args.bam)
    if args.upper and not lengths and not args.quiet:
        note("no contig lengths available; an upper flank past a contig end will not be clipped")

    # target name -> its intervals, in the order the targets were given.
    resolved: Dict[str, List[Interval]] = {}
    order: List[str] = []

    if args.id:
        wanted = {q: id_kind(q) for q in args.id}
        found = scan_gtf(args.gtf, wanted, args.quiet)
        for query, kind in wanted.items():
            hits = found[query]
            # An ID that matched nothing is an error. An empty query would run cleanly and report
            # zero coverage, which reads exactly like real data.
            if not hits:
                die(
                    f"--id {query!r} ({kind}) matched no exon in {args.gtf}. Check the ID, and "
                    f"that the GTF is the annotation release you think it is"
                )
            mode = args.mode
            if mode == "auto":
                mode = "merged_exon" if kind in ("gene", "gene_name") else "exon"
            if mode == "merged_exon":
                core = merge(hits, query)
            elif mode == "span":
                core = span(hits, query)
            else:
                core = dedupe(hits, query, args.quiet)
            resolved[query] = core + flanks(core, query, args.lower, args.upper, lengths)
            order.append(query)

    for spec in args.region:
        iv = parse_region(spec)
        iv = extend(iv, args.lower, args.upper, lengths)
        resolved[iv[3]] = [iv]
        order.append(iv[3])

    intervals = [iv for q in order for iv in resolved[q]]

    if not args.quiet:
        for q in order:
            core = [iv for iv in resolved[q] if not iv[5].startswith("flank")]
            fl = [iv for iv in resolved[q] if iv[5].startswith("flank")]
            bases = sum(e - s for _, s, e, _, _, _ in core)
            msg = f"{q}: {len(core)} interval(s), {bases} bp"
            if fl:
                msg += f" + {sum(e - s for _, s, e, _, _, _ in fl)} bp flanking"
            note(msg)

    blocks = union(intervals)
    handle = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
    handle.close()
    try:
        with open(handle.name, "w") as fh:
            for contig, start, end in blocks:
                fh.write(f"{contig}\t{start}\t{end}\n")
        rows = run_depth(handle.name, args.bam, args, blocks, list(lengths))
        # `-a` reports every queried base exactly once, so this equality is the check that the
        # regions were unioned and that nothing was dropped. A short count would otherwise surface
        # as a track with a hole in it.
        queried = sum(end - start for _, start, end in blocks)
        if len(rows) != queried:
            die(
                f"samtools depth returned {len(rows)} positions for {queried} queried bases — "
                f"the region list is not disjoint, or depth dropped positions"
            )
        for r in rows:
            if len(r) != 2 + len(args.bam):
                die(f"unexpected samtools depth output with {len(r)} columns: {r!r}")
        tracks = write_bedgraph(
            args.outdir, rows, labels, not args.no_track_line
        )
    finally:
        os.unlink(handle.name)

    if not args.quiet:
        for path, n in tracks:
            note(f"{n} interval(s) -> {path}")
        note(
            "zero-depth bases are written as 0-valued intervals, so a gap in a track means the "
            "position was not queried, never that it had no reads"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
