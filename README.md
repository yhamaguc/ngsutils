ngsutils
========
Utilities for NGS processing

## Installation

```bash
$ git clone https://github.com/yhamaguc/ngsutils.git
$ cd ngsutils && pip install .
```

For development, install in editable mode so that changes apply without reinstalling:

```bash
$ cd ngsutils && pip install -e .
```

`pip install .` reuses `build/lib/` from earlier builds and copies whatever is there into the
wheel, so a module deleted from the repository keeps being installed. Measured on 2026-08-21:
`align_star.py` and `sort_star_bam.py`, removed in 72a8cca, were still in site-packages after a
fresh `pip install .`. Remove the tree when files have been deleted or renamed:

```bash
$ rm -rf build ngsutils.egg-info && pip install .
```

## Usage

```bash
$ ngsutils --help              # List subcommands
$ ngsutils <subcommand> --help # Show help for a subcommand
```

Every subcommand is installed under its own name as well, so the two forms below are the same
program and take the same arguments:

```bash
$ ngsutils id2name --gtf gencode.v50.gtf.gz < ids.tsv
$ id2name --gtf gencode.v50.gtf.gz < ids.tsv
```

The bare names are how these utilities were installed before 0.2.0; the move to `pyproject.toml`
dropped them for one release, which broke callers that used them. `[project.scripts]` and the
subcommand table in `ngsutils/cli.py` must list the same commands, and `tests/test_entry_points.py`
fails if they drift apart.

## Annotation maps

`id2name` and `name2id` take `--gtf`. Without it they look up an identifier or a name in a
prebuilt map, and that map is **not in the repository**: at 23 MB per file per GENCODE
release, committing it would put 9 MB of packed objects into every clone for every release
ever built. The maps are GitHub release assets instead, and what is committed is a 1.6 KB
manifest, `ngsutils/data/maps.json`:

| Field | What it is |
|---|---|
| `release` | `gencode.v50` — the annotation the maps were built from |
| `repository`, `tag` | where the assets live: `yhamaguc/ngsutils`, tag `maps-gencode.v50` |
| `assets.<map>.file` | `id2name_gencode.v50.pkl.gz`, 4.5 MB on the wire |
| `assets.<map>.archive_sha256`, `.pickle_sha256` | checked before anything is unpickled |
| `provenance` | source path and MD5 of the GTF, its header, row count, build route, git revision |

`ngsutils/maps.py` is the only code that decides what a map is. It resolves one in three
steps and says on stderr which it took:

1. a plain pickle vendored into `ngsutils/data` by hand, for an install with no network
2. `${XDG_CACHE_HOME:-~/.cache}/ngsutils/maps/<map>_<release>.pkl`, the unpacked cache
3. the release asset, downloaded once (4.5 MB), verified, and unpacked into that cache

**The SHA-256 in the manifest is checked before the bytes are unpickled**, and again every
time the cache is read. A pickle executes on load, so a substituted or truncated file has to
be an error rather than a wrong answer; it is, and the message names the expected and found
digest. There is no silent fallback: if the download fails, the command exits non-zero with
the URL it tried and tells the caller to pass `--gtf` instead.

Three properties of the maps are worth knowing before trusting a lookup:

- **Identifiers are truncated to 15 characters**, which drops the version suffix, so
  `ENSG00000141510.21` and `ENSG00000141510` both resolve. Nothing collides under that
  truncation in the primary assembly file, which carries no `PAR_Y` entries, but the full
  `gencode.vNN.annotation.gtf` does carry them and its chrY copies collapse onto their chrX
  counterparts. They share a name, so the map stays single-valued; the build fails rather
  than picks if two rows ever disagree on the value of one key.
- **`name2id` is not injective in reverse.** 484 gene names in v50 are shared by more than
  one gene identifier (`5S_rRNA`, `U1`, `Y_RNA` and the like); the map keeps the first
  identifier in sorted order and the build reports the count. Passing `--gtf` instead
  resolves the same name from row order, which is not defined, so the two routes can
  disagree on a reused name — measured with `DDX11L16`, which v50 puts at five loci. The
  prebuilt map is the reproducible of the two.
- **v50 is much larger than v45**, the release this replaced: 316,118 keys to 725,519.
  411,473 of those are new, and 410,002 of the 646,577 transcripts in the input are
  TAGENE — long-read derived. 2,072 identifiers the v45 map answered are absent, and 11,150
  names changed (`C2orf83` → `SLC19A4P`). Some of that is scope rather than annotation: the
  v45 map came from a different input file.

### Rebuilding for a newer release

```bash
$ bin/build_id_name_maps.py --gtf gencode.v51.primary_assembly.annotation.gtf.gz
$ gh release create maps-gencode.v51 --repo yhamaguc/ngsutils \
      dist/maps/id2name_gencode.v51.pkl.gz dist/maps/name2id_gencode.v51.pkl.gz
$ git add ngsutils/data/maps.json && git commit
```

The build reads the release label out of the GTF header, writes the two archives into
`dist/maps` (gitignored, they are assets) and the manifest into `ngsutils/data`, and prints
the `gh release create` line to run. It reads the input **one row at a time**: 46 s and
0.24 GB of resident memory for v50, against 25 s and 9.06 GB through `ngsutils.gtf.read_gtf`,
which is the route `--gtf` takes at run time. `--cross-check` runs that route too and
asserts the two agree pair for pair; they pickle to identical bytes.

Downloading the input is the slow part, not the build: 125 MB at 927 kB/s through a
proxied link is 134 s, so zero to ready is about three minutes.

### Earlier releases

Every release built stays downloadable under its own tag, so bumping the manifest does not
take an older map away from anyone who was using it:

| Tag | Maps | Provenance |
|---|---|---|
| [`maps-gencode.v50`](https://github.com/yhamaguc/ngsutils/releases/tag/maps-gencode.v50) | 725,519 and 723,889 keys | built from `gencode.v50.primary_assembly.annotation.gtf.gz`, cross-checked |
| [`maps-gencode.v45`](https://github.com/yhamaguc/ngsutils/releases/tag/maps-gencode.v45) | 316,118 and 314,509 keys | the bytes that used to be committed, recovered from git blob `bf793a4` |

Each release carries its own `maps.json` with the digests, so an asset can be checked without
this repository:

```bash
$ curl -fLO https://github.com/yhamaguc/ngsutils/releases/download/maps-gencode.v45/maps.json
$ curl -fLO https://github.com/yhamaguc/ngsutils/releases/download/maps-gencode.v45/id2name_gencode.v45.pkl.gz
$ sha256sum id2name_gencode.v45.pkl.gz   # against archive_sha256 in maps.json
```

**The v45 maps have no recorded provenance.** A pickle carries none, and the GTF those were
built from was never written down, so `maps.json` for that release says so rather than
guessing. They also predate the deterministic collision rule, so a reused gene name may
resolve to a different identifier than a rebuild would give. To put an old map on a release
without inventing a source for it, `scripts/package_map_archives.py` does the packaging half
alone.

### Measurements

Every figure above comes from a script in `scripts/`, and none of them are inline
one-liners. Run them to reproduce:

```bash
$ scripts/bench_map_load.py --map dist/maps/id2name_gencode.v50.pkl.gz \
      --output dist/measurements/map_load.tsv
$ scripts/bench_build_routes.py --gtf gencode.v50.primary_assembly.annotation.gtf.gz \
      --output dist/measurements/build_routes.tsv
$ git show HEAD:ngsutils/data/id2name_gencode.v45.pkl > /tmp/before.pkl
$ scripts/compare_id_name_maps.py --before /tmp/before.pkl \
      --after dist/maps/id2name_gencode.v50.pkl.gz \
      --output dist/measurements/id2name_v45_to_v50.tsv
```

Measured 2026-08-21, Python 3.12.12, 15 interleaved rounds for the load figures. The map is
kept as a pickle rather than a gzipped TSV because of the first of those: 244 ms against
408 ms at the median, on a command whose whole run is 520 ms, of which 484 ms is importing
polars — which is why both commands defer that import to the `--gtf` branch that needs it.

## External programs

Most subcommands need only the Python dependencies in `pyproject.toml`. `write_coverage_track`
also needs **`samtools` on `PATH`**, which `pyproject.toml` cannot express because it is not a
Python package. It is checked at run time and a missing `samtools` is an error, not a fallback.

## What moved out

Commit 72a8cca reduced this repository to what its name says: a Python package of annotation,
identifier and IGV utilities. Two things that were never that went elsewhere, and their
documentation went with them:

| Was here | Now | What it is |
|---|---|---|
| `bin/*.sh`, `bin/*.py`, `bin/*.R` | `ngstool-recipes` | Wrappers that shell out to an installed tool (STAR, salmon, fastp, DESeq2, …), each carrying a preset. They parse their arguments with [`docopts`](https://github.com/docopt/docopts) — the Go binary, not the `docopt` Python package this project depends on — and that repository documents installing it. |
| `ngsutils/align_star.py`, `ngsutils/sort_star_bam.py` | `rnaseq-snakemake` | Designed with the Snakemake rules that call them, and revised with those rules. They are called from that checkout by path, so the version that runs is the one beside the Snakefile. |

Nothing in this repository needs `docopts`. Every utility here is a `[project.scripts]` console
script. Two directories hold things that are not utilities:

- `bin/build_id_name_maps.py` builds the annotation maps. It is installed, through a
  `[tool.setuptools] script-files` entry — an allowlist, so a file left in `bin/` without an
  entry there is never installed.
- `scripts/` holds the measurement scripts the numbers in this README come from, and
  `package_map_archives.py`, which wraps an already-built map into a release asset. Nothing
  there is installed or packaged; they exist so a figure can be re-derived rather than trusted.

## Tests

```bash
$ python3 -m unittest discover -s tests   # stdlib only, no extra dependency
$ pytest tests/                           # equivalent, if pytest is installed
```

The tests build every input they need in a temporary directory — SAM converted to an indexed BAM,
and a small GTF — so nothing but `tests/data/sample.gff` is a fixture on disk. Tests needing
`samtools` skip when it is absent rather than failing.

`tests/test_id_name_maps.py` asserts the shape of `ngsutils/data/maps.json`, that the asset
names follow from the release it declares, that no map is *committed* (the cost the release
assets exist to avoid), and that a wrong-hashed file is refused rather than unpickled. **It
never downloads anything**: the tests that need a real map skip unless one is already vendored
or cached, so a clean checkout with no cache reports skips rather than failures. When a map is
resolvable it also checks four stable genes both ways and one gene GENCODE renamed after v45 —
the one check a relabelled old map cannot pass.

`tests/test_entry_points.py` reads `pyproject.toml` and the docstrings instead, and needs neither:
it asserts that the console scripts and the subcommand table agree, that every command's docopt
options and `[default:]` values are the intended ones, and the same for the standalone scripts in
`bin/` and `scripts/`, whose docstrings it reads off disk without importing them.

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
