ngsutils
========
Utilities for NGS processing

## Installation

```bash
$ git clone https://github.com/yh549848/ngsutils.git
$ cd ngsutils && pip install .
```

For development, install in editable mode so that changes apply without reinstalling:

```bash
$ cd ngsutils && pip install -e .
```

## Usage

```bash
$ ngsutils --help              # List subcommands
$ ngsutils <subcommand> --help # Show help for a subcommand
```

## External programs

Most subcommands need only the Python dependencies in `pyproject.toml`. `write_coverage_track`
also needs **`samtools` on `PATH`**, which `pyproject.toml` cannot express because it is not a
Python package. It is checked at run time and a missing `samtools` is an error, not a fallback.

## Tests

```bash
$ python3 -m unittest discover -s tests   # stdlib only, no extra dependency
$ pytest tests/                           # equivalent, if pytest is installed
```

The tests build every input they need in a temporary directory — SAM converted to an indexed BAM,
and a small GTF — so nothing but `tests/data/sample.gff` is a fixture on disk. Tests needing
`samtools` skip when it is absent rather than failing.

## `write_coverage_track` — per-base coverage as bedGraph

Imported from another repository on 2026-08-12, where it was `tools/write_coverage_track.py`. One
job: resolve a set of target intervals, then write one bedGraph per BAM over exactly those
intervals. Nothing is summarised and nothing else is written.

```bash
# by ID (needs --gtf)
$ ngsutils write_coverage_track --output-dir tracks --gtf gencode.v50.gtf.gz \
    --id ENST00000269305 sample.bam

# by coordinate, 1-based inclusive as samtools writes regions
$ ngsutils write_coverage_track --output-dir tracks --region chr17:7673700-7673837 sample.bam
```

**BAMs are positional**, as `samtools depth [options] in.bam [in.bam ...]` takes them, so a file
list composes with no argument rewriting:

```bash
$ find data -name '*.bam' | sort \
| xargs ngsutils write_coverage_track --output-dir tracks --gtf ... \
        --id ENST00000269305 --label-from parent-stem
```

One `samtools depth` call per contig takes **all** the BAMs and reports one column per file at each
position, so they are read together. Every call carries `-r` so samtools seeks via the index
instead of reading the whole file — `depth -b BED` alone does **not** use the index, which is why a
BAM without one is a hard error naming `samtools index`.

### Modes

| `--mode` | Intervals | |
|---|---|---|
| `span` | the locus, **introns included** | **the default** |
| `merged_exon` | union of every isoform's exons, each base once | |
| `exon` | one per exon | |
| `auto` | exonic: `merged_exon` for a gene or gene name, `exon` otherwise | |

The default spans the locus because a track that stops at the exons cannot show a retained intron.
`samtools depth` does not count CIGAR `N`, so a read that spliced across an intron contributes 0
there and only a read aligned *through* it leaves depth. Measured on a three-read fixture: exon
depth 3, intron depth 1. `-J` does not change that — it admits `D` and not `N`.

The cost: TP53 `ENST00000269305` spans 19,070 bp with 2,512 bp of exon, so a `span` track is about
7.6x the positions and mostly intron. Use `merged_exon` or `exon` when that is not wanted.

**This is not a retained-intron caller.** Depth alone does not separate a retained intron from
nascent pre-mRNA, from an overlapping feature on either strand — this program does not split by
strand — or from the primary alignment of a multi-mapping read in a repeat, where `--min-mq` is the
only lever. The `N`/`M` distinction is the aligner's, not samtools': a read falling wholly inside
an intron carries no junction and is `M` whatever produced it.

### Sums and means are derived, not emitted

A per-exon mean is a length-weighted sum over this track, so the track is the primitive and a
summary is a query over it:

```bash
$ bedtools intersect -a tracks/s1.bedgraph -b exons.bed -wb \
| awk -F'\t' '{s[$8] += ($3-$2)*$4} END {for (k in s) print k, s[k]}'
```

Note that `bedtools map -c 4 -o sum` is **wrong** here and returns 2 where the answer is 20: it
sums the bedGraph values without weighting by interval length. The identity against
`samtools bedcov -j` is asserted in the tests.

Zero is written, never omitted. A missing interval in a bedGraph means *no data*, so a zero-depth
base left out would read as unmeasured. `ngsutils write_coverage_track --help` has the rest.
