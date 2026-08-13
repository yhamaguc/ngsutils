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

Every script in `bin/` needs **`docopts` on `PATH`** for the same reason, plus whatever tool it
wraps (STAR, salmon, fastp, the UCSC utilities, …). See below.

## `bin/` scripts and `docopts`

The shell scripts in `bin/` declare their command line in a docopt usage string and parse it with
[`docopts`](https://github.com/docopt/docopts). That is the **Go binary**, not the `docopt` Python
package listed in `pyproject.toml` — they share the docopt language but not the implementation, and
installing one does not provide the other.

There is no conda package for `docopts`, so it cannot be added to an environment file the way
`samtools` can. Install it by dropping the release binary on `PATH`:

```bash
$ VERSION=v0.6.4-with-no-mangle-double-dash
$ curl -fsSL -o ~/.local/bin/docopts \
    https://github.com/docopt/docopts/releases/download/${VERSION}/docopts_linux_amd64
$ chmod +x ~/.local/bin/docopts
```

The release publishes a `sha256sum.txt` beside the binaries; `docopts_linux_amd64` for that version
is `a2e565b6f2ca73103005bc62295ed3401b40f10772f205587778595aeaa2b8fb`. Verify against it rather than
trusting the download. With a Go toolchain, `go install github.com/docopt/docopts@${VERSION}` builds
the same program instead.

Where the compute nodes are not guaranteed to match the submit host, bake `docopts` into the
container image rather than installing it per node — that is the only option here that survives a
node the scripts have never run on.

Each script checks for `docopts` before it does anything else and exits 127 with a message naming
it, so a missing binary fails immediately instead of part way through a run:

```
Error: docopts not found on PATH. bin/*.sh parse their arguments with
       it; install it from https://github.com/docopt/docopts
```

A malformed command line gets the usage message and exit 64; `--help` prints the full interface.

Note that `pip install .` puts most of `bin/` on `PATH` through `script-files` in
`pyproject.toml`, but `trim_fastp.sh` and `conv_sjouttab2bed.py` are not listed there and have to be
run from a checkout.

### Gzipped references

Every `build_*.sh` takes a `.gz` FASTA or GTF and expands it before the indexer sees it, so a
GENCODE download can be used as it arrives:

```bash
$ build_star.sh --gtf=gencode.v50.annotation.gtf.gz --output-dir=index GRCh38.primary_assembly.fa.gz
```

The indexes are named after the decompressed files, so `GRCh38.fa.gz` and `GRCh38.fa` produce the
same output directory. This was measured, not assumed: STAR indexes built from the gzipped and
plain forms of the same reference compare byte-identical, as do salmon's.

The tools disagree about compressed input and only some of them say so. `salmon index` reads `.gz`
directly and `kallisto index` documents that it does; `STAR --runMode genomeGenerate` refuses one
with `Make sure the file is uncompressed (unzipped)`; `rsem-prepare-reference` fails with **`Number
of transcripts in the reference is less than 1!`**, which never mentions compression and is the
real reason this is handled in the scripts rather than left to the caller. `hisat2-build` does not
document `.gz` support either way. Rather than keep a per-tool table correct against five tools
that change independently, every script expands `.gz` the same way.

The expanded copy goes in a private directory beside `--output-dir`, not under `$TMPDIR` — an
uncompressed primary assembly is a few GB and `/tmp` on a compute node rarely holds one — and is
removed on exit, including on failure. Budget the space where the index is being written.

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
