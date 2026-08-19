#!/usr/bin/env python3
"""Align RNA-seq reads to a genome with STAR, from a BAM or from FASTQ.

One entry point, two inputs. An <input> ending in .bam is re-extracted to FASTQ with
biobambam2's bamtofastq, streamed through FIFOs so nothing lands on disk; anything else is
treated as read 1 FASTQ and handed to STAR with --readFilesCommand. Everything downstream
of that -- presets, sjdb mode, genome load, sorting, verification -- is identical for both,
and that identity is checked: on the same 40,000-record TCGA-UVM BAM the two paths produce
byte-identical alignment records.

This replaces a pair of shell scripts (a BAM-only align_star.sh and a FASTQ-only
align_star_fastq.sh) that had drifted apart in both features and style. The move to Python
was for two specific reasons, both of which had already cost a run:

  1. `docopts` drops an option's [default:] without a word of complaint when another
     declared option's name appears, with its leading --, inside the option's DESCRIPTION.
     Measured 2026-08-19: --sort-with and --sort-ram both parsed as EMPTY strings. Python's
     docopt does not have that failure mode, and the tests below pin the defaults anyway.
  2. The branching -- BAM or FASTQ, sjdb-only or aligned, STAR's sorter or samtools',
     per-sample 2-pass or not -- has no unit test in shell and cannot easily get one.
     tests/test_align_star.py covers the parts that decide what STAR is asked to do.

Sorting is a separate step
==========================
STAR's sorter is in-memory and bounded by --limitBAMsortRAM, its requirement scales with
read count, and it fails only AFTER the whole mapping phase has been paid for. Any fixed
ceiling is therefore one large sample away from an expensive crash. `--sort-with none`
writes Aligned.out.bam and stops; `ngsutils sort_star_bam` then sorts it with samtools,
which spills to disk. Splitting them means the memory-hungry half can be retried, given
different resources, or moved to another host without repeating the alignment -- and in a
workflow it becomes its own rule with its own resource request.

Two 2-pass schemes, which are not substitutes
=============================================
The --make-sjdb flag drives a COHORT-level 2-pass: run it over every sample, merge the SJ.out.tab
files, rebuild one index from them, then align against that. --two-pass is STAR's own
PER-SAMPLE 2-pass (--twopassMode Basic), as GDC runs it. They compose. The cohort scheme
requires that both of its passes use the SAME --preset: the 1st pass decides which
junctions enter the sjdb, so changing filters between passes silently changes what the
rebuilt index contains.

Usage:
  align_star [options] <input> <genome-dir> <out-prefix> [<threads>]
  align_star (-h | --help)

Arguments:
  <input>       Read 1 FASTQ, or a BAM. A .bam suffix selects the BAM path, where reads
                are re-extracted with bamtofastq; anything else is FASTQ.
  <genome-dir>  STAR genome index directory.
  <out-prefix>  Prefix prepended to every STAR output file. A trailing / makes it a
                directory.
  <threads>     Threads [default: 16]

Options:
  --fastq2=<fastq2>            Mate 2 FASTQ. FASTQ input only; omit for single-end. With
                               a BAM the layout is read from the file itself.
  --make-sjdb                  1st-pass mode: emit only <out-prefix>SJ.out.tab and skip
                               alignment output. No BAM is written, so nothing has to be
                               deleted afterwards.
  --two-pass                   STAR's own per-sample 2-pass. Roughly doubles runtime.
                               Forces a fallback to genome-load NoSharedMemory, and on the
                               BAM path COSTS DISK: STAR reads its input twice, which a
                               FIFO cannot serve, so the FASTQ is written out next to
                               <out-prefix> instead of streamed.
  --preset=<preset>            Parameter set: none, encode, gdc. [default: none]
  --genome-load=<genome-load>  STAR --genomeLoad: NoSharedMemory, LoadAndKeep,
                               LoadAndRemove, LoadAndExit, Remove. Falls back to
                               NoSharedMemory, loudly, when the run needs run-time
                               junction insertion. [default: LoadAndKeep]
  --sort-with=<sort-with>      star: STAR sorts internally, bounded by the sort-ram
                               ceiling. none: stop after writing Aligned.out.bam and
                               leave the coordinate sort to `ngsutils sort_star_bam`,
                               which spills to disk. Use none whenever the sample might
                               be large -- see "Sorting" above. [default: star]
  --sort-ram=<sort-ram>        Bytes for STAR's internal BAM sort. Ignored when sort-with
                               is none. [default: 160000000000]
  --layout=<layout>            paired, single, or auto. auto reads the layout out of the
                               BAM with one samtools flagstat, which decodes it whole -- once
                               per invocation, so twice per sample across a 1st and 2nd
                               pass, on top of the two extractions. A
                               caller that already knows the layout (a sample table has it)
                               should say so. FASTQ input ignores this: whether fastq2 was
                               given IS the layout there. [default: auto]
  --read-files-command=<cmd>   Decompressor for FASTQ input, e.g. 'zcat' or 'pigz -d -c'.
                               Auto-selected from the .gz suffix when not given; pass
                               'none' to force no decompression. FASTQ input only.
                               [default: auto]
  -h --help                    Show this message.

Defaults measured on CCLE data, 16 threads, 3 samples, matched pairs:
  pigz -d -c saves a median 4.0% of mapping time over zcat (range -4.5% to +9.7%).
  LoadAndKeep cuts the genome load from ~23 s to ~5 s from the second run onward.
"""

from __future__ import annotations

import os
import re
import signal
import shlex
import subprocess
import sys
import time
from pathlib import Path

from docopt import docopt

# %%
# Constants

DEFAULT_THREADS = 16

# pigz beat zcat on 4 of 6 matched pairs, median 4.0% of mapping time. It is absent from
# the STAR image itself and reaches the container through the bind mounts in the STAR
# shim; if that shim is replaced, this default breaks loudly (0 input reads), which
# assert_reads_were_mapped() catches.
DEFAULT_GZ_COMMAND = "pigz -d -c"

BAM_SUFFIX = ".bam"
GZIP_SUFFIX = ".gz"
NO_DECOMPRESSION = "none"
AUTO_DECOMPRESSION = "auto"

# "none" leaves the coordinate sort to ngsutils sort_star_bam, which owns it.
LAYOUTS = ("paired", "single", "auto")
AUTO_LAYOUT = "auto"

SORTERS = ("star", "none")
NO_SORT = "none"

# STAR refuses any shared-memory genome when junctions have to be inserted at run time,
# because the segment is sized to the on-disk index with no headroom (measured:
# 30,946,456,755 bytes of shm against 30,946,456,451 bytes of Genome + SA + SAindex) and is
# shared with other processes that would see it mutate.
SHARED_MEMORY_LOADS = ("LoadAndKeep", "LoadAndRemove", "LoadAndExit")
NO_SHARED_MEMORY = "NoSharedMemory"

SORTED_BAM_SUFFIX = "Aligned.sortedByCoord.out.bam"
UNSORTED_BAM_SUFFIX = "Aligned.out.bam"
TRANSCRIPTOME_BAM_SUFFIX = "Aligned.toTranscriptome.out.bam"
SJ_SUFFIX = "SJ.out.tab"
LOG_FINAL_SUFFIX = "Log.final.out"
STAR_TMP_SUFFIX = "_STARtmp"

# --quantMode TranscriptomeSAM opens this file, which only exists in a GTF-built index.
# Without it STAR dies at exit 109 after the FASTQ has already been streamed.
GENE_INFO = "geneInfo.tab"

INPUT_READS_PATTERN = re.compile(r"Number of input reads\s*\|\s*(\d+)")

# How often run_watched() checks on the two processes. Small enough that a dead
# extraction is caught in about a second, large enough to cost nothing over hours.
WATCHDOG_INTERVAL_SECONDS = 1.0

# samtools flagstat's first field is QC-passed; both counts come from one decode.
#
# NOTE: the denominator is PRIMARY, not "in total". Measured 2026-08-20 on a real TCGA-UVM
#   BAM: in total 40,000 = primary 15,295 + secondary 24,705, while "paired in sequencing" is
#   15,295 -- it counts primary records only. Comparing paired against "in total" therefore
#   reports every multimapping paired BAM as a mix of paired and unpaired reads, which would
#   have refused perfectly ordinary GDC input. Against primary the two agree exactly.
# NOTE: BOTH fields are captured and summed. bamtofastq does not exclude QC-failed records
#   (0x200), so the counts have to describe what it will actually extract -- reading the
#   QC-passed field alone would call a BAM whose paired records are all QC-failed
#   single-end, and then extract and map its mates as independent reads.
PRIMARY_PATTERN = re.compile(r"^(\d+) \+ (\d+) primary$", re.MULTILINE)
PAIRED_PATTERN = re.compile(r"^(\d+) \+ (\d+) paired in sequencing", re.MULTILINE)

# %%
# Presets
#
# A preset carries alignment parameters only. Anything this module owns -- --runThreadN,
# --genomeDir, --readFilesIn, --readFilesCommand, --outFileNamePrefix, --outSAMtype,
# --genomeLoad, --limitBAMsortRAM, --twopassMode -- is never in here.
#
# `common` shapes WHICH junctions STAR reports, so it must be identical in both passes of a
# cohort-level 2-pass: an sjdb built under one filter regime and consumed under another is
# not the sjdb the second pass thinks it is. `second` is output formatting, inert or fatal
# under --outSAMtype None, so it is withheld from the sjdb-only pass.

ENCODE_COMMON = (
    # Verbatim from the STAR manual 3.3.2 "ENCODE options" (p.9). All nine are
    # filter/alignment parameters, so none of them belong in `second`.
    "--outFilterType", "BySJout",
    "--outFilterMultimapNmax", "20",
    "--alignSJoverhangMin", "8",
    "--alignSJDBoverhangMin", "1",
    "--outFilterMismatchNmax", "999",
    "--outFilterMismatchNoverReadLmax", "0.04",
    "--alignIntronMin", "20",
    "--alignIntronMax", "1000000",
    "--alignMatesGapMax", "1000000",
)

# Transcribed from the @PG CL line of a GDC-realigned TCGA BAM
# (STAR 2.7.5c, star-2.7.5c_GRCh38.d1.vd1_gencode.v36).
#
# Deliberately NOT carried over, with the reason for each:
#   --runThreadN, --genomeDir, --readFilesIn, --outFileNamePrefix
#                            arguments of this command
#   --readFilesCommand zcat  --read-files-command owns it, and the BAM path's FIFOs carry
#                            uncompressed FASTQ
#   --outSAMtype BAM Unsorted
#                            owned by --make-sjdb and --sort-with
#   --genomeLoad NoSharedMemory
#                            --genome-load owns it; its default differs deliberately
#   --outSAMattrRGline       per-sample (ID/SM/LB/PU), not a preset
#   --twopassMode Basic      moved out to --two-pass, so the per-sample 2-pass is chosen
#                            independently of the GDC parameter set
#   --quantMode GeneCounts   not needed here; TranscriptomeSAM is kept
#   --chim*                  --chimOutType WithinBAM is a fatal parameter error under
#                            --outSAMtype None (exit 102), and chimeric detection is out
#                            of scope
GDC_COMMON = (
    "--limitSjdbInsertNsj", "1200000",
    "--outFilterType", "BySJout",
    "--outFilterMultimapNmax", "20",
    "--outFilterScoreMinOverLread", "0.33",
    "--outFilterMatchNminOverLread", "0.33",
    "--outFilterMismatchNmax", "999",
    "--outFilterMismatchNoverLmax", "0.1",
    "--outFilterIntronMotifs", "None",
    "--alignIntronMin", "20",
    "--alignIntronMax", "1000000",
    "--alignMatesGapMax", "1000000",
    "--alignSJoverhangMin", "8",
    "--alignSJDBoverhangMin", "1",
    "--alignSoftClipAtReferenceEnds", "Yes",
)

GDC_SECOND = (
    "--outSAMstrandField", "intronMotif",
    "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "ch",
    "--outSAMunmapped", "Within",
    "--quantMode", "TranscriptomeSAM",
)

PRESETS = {
    "none": {"common": (), "second": (), "needs_gene_info": False},
    "encode": {"common": ENCODE_COMMON, "second": (), "needs_gene_info": False},
    "gdc": {"common": GDC_COMMON, "second": GDC_SECOND, "needs_gene_info": True},
}


class AlignError(RuntimeError):
    """A refusal or a failure that main() turns into a message and a non-zero exit."""


# %%
# Decisions
#
# Everything in this section is a pure function of its arguments, which is the point: these
# are what decide what STAR is asked to do, and they are what tests/test_align_star.py
# pins down.


def input_kind(input_path: str | Path) -> str:
    """"bam" or "fastq", from the suffix alone.

    The suffix is the whole decision. A BAM named without .bam would be handed to STAR as
    FASTQ and fail on binary input -- name files honestly.
    """
    return "bam" if str(input_path).endswith(BAM_SUFFIX) else "fastq"


def resolve_genome_load(genome_load: str, two_pass: bool) -> tuple[str, str | None]:
    """(the value to use, a note to print or None).

    Falls back off shared memory rather than letting STAR abort, because --twopassMode
    needs run-time junction insertion. The note is not optional decoration: the run then
    costs ~18 s more of genome load, and silence would make that unattributable.
    """
    if two_pass and genome_load in SHARED_MEMORY_LOADS:
        return NO_SHARED_MEMORY, (
            f"--two-pass needs run-time junction insertion, which STAR forbids on a "
            f"shared-memory genome. Falling back from {genome_load} to {NO_SHARED_MEMORY}."
        )
    return genome_load, None


def resolve_read_files_command(input_path: str | Path, requested: str) -> list[str]:
    """The --readFilesCommand argv fragment, empty when no decompression is wanted.

    'auto' inspects the suffix only. A gzipped file named without .gz would be read as
    plain text and STAR would fail on binary input -- name files honestly.
    """
    if requested == AUTO_DECOMPRESSION:
        command = (
            DEFAULT_GZ_COMMAND if str(input_path).endswith(GZIP_SUFFIX)
            else NO_DECOMPRESSION
        )
    else:
        command = requested

    if command == NO_DECOMPRESSION:
        return []
    # shlex so 'pigz -d -c' reaches STAR as three arguments.
    return ["--readFilesCommand", *shlex.split(command)]


def expected_outputs(
    out_prefix: str, make_sjdb: bool, needs_gene_info: bool, sort_with: str = "star"
) -> list[str]:
    """What a successful run must have written.

    --outSAMtype None still writes SJ.out.tab and Log.final.out, so the sjdb-only pass
    never materialises a BAM in the first place. Under --sort-with none the genome BAM is
    the UNSORTED one, because sorting is a later step.
    """
    if make_sjdb:
        return [f"{out_prefix}{SJ_SUFFIX}"]
    outputs = [
        f"{out_prefix}{UNSORTED_BAM_SUFFIX}" if sort_with == NO_SORT
        else f"{out_prefix}{SORTED_BAM_SUFFIX}"
    ]
    if needs_gene_info:
        outputs.append(f"{out_prefix}{TRANSCRIPTOME_BAM_SUFFIX}")
    return outputs


def build_star_argv(
    *,
    genome_dir: str,
    out_prefix: str,
    read_files: list[str],
    threads: int,
    genome_load: str,
    preset: str,
    make_sjdb: bool,
    two_pass: bool,
    sort_with: str,
    sort_ram: str,
    read_files_command: list[str],
) -> list[str]:
    """The full STAR command line.

    The sorted BAM has the same filename whichever sorter produced it, so design sheets and
    rMATS input lists do not have to know.
    """
    argv = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", genome_dir,
        "--genomeLoad", genome_load,
        "--readFilesIn", *read_files,
        "--outFileNamePrefix", out_prefix,
        *read_files_command,
    ]

    if make_sjdb:
        argv += ["--outSAMtype", "None"]
    elif sort_with == NO_SORT:
        argv += ["--outSAMtype", "BAM", "Unsorted"]
    else:
        argv += ["--outSAMtype", "BAM", "SortedByCoordinate"]
        # --limitBAMsortRAM defaults to 0, which STAR reads as "the size of the genome
        # index" -- about 30 GB for a human index. A CCLE run of ~90M read pairs needs
        # ~40 GB, so the default aborts AFTER the whole mapping phase has been paid for.
        # History of this default, which is the point: the index-size default failed on
        # 40 GB requests; 64 GB failed on a sample that asked for 132,227,075,545 bytes
        # after 83 min of mapping. Any fixed ceiling is one large sample away from
        # breaking, because STAR's sort is in-memory and scales with read count. The
        # durable answer is --sort-with none plus ngsutils sort_star_bam.
        argv += ["--limitBAMsortRAM", str(sort_ram)]

    if two_pass:
        argv += ["--twopassMode", "Basic"]

    argv += list(PRESETS[preset]["common"])
    if not make_sjdb:
        argv += list(PRESETS[preset]["second"])
    return argv


def parse_input_reads(log_final: str) -> int | None:
    """`Number of input reads` out of Log.final.out, or None when it is not there."""
    match = INPUT_READS_PATTERN.search(log_final)
    return int(match.group(1)) if match else None


# %%
# Validation


def validate(
    *,
    input_path: Path,
    genome_dir: Path,
    fastq2: str | None,
    read_files_command: str,
    sort_with: str,
    preset: str,
    layout: str = AUTO_LAYOUT,
) -> str:
    """Refuse everything that would otherwise fail late or silently. Returns the kind."""
    if preset not in PRESETS:
        raise AlignError(
            f"unknown --preset {preset!r} (expected: {', '.join(sorted(PRESETS))})"
        )
    if sort_with not in SORTERS:
        raise AlignError(
            f"unknown --sort-with {sort_with!r} (expected: {', '.join(SORTERS)})"
        )
    if layout not in LAYOUTS:
        raise AlignError(f"unknown --layout {layout!r} (expected: {', '.join(LAYOUTS)})")
    if not input_path.exists():
        raise AlignError(f"input not found: {input_path}")
    if not genome_dir.is_dir():
        raise AlignError(f"genome dir not found: {genome_dir}")

    kind = input_kind(input_path)

    # The two FASTQ-only options. Accepting them silently on the BAM path would be worse
    # than refusing: the FIFOs carry plain uncompressed FASTQ and the layout comes from the
    # BAM, so either option would have no effect and the caller would believe otherwise.
    if kind == "bam":
        if fastq2:
            raise AlignError(
                f"--fastq2 is FASTQ-only; the BAM path reads the layout from {input_path}"
            )
        if read_files_command != AUTO_DECOMPRESSION:
            raise AlignError(
                "--read-files-command is FASTQ-only; the FIFOs from bamtofastq carry "
                "uncompressed FASTQ and STAR is given no decompressor for them"
            )
    elif fastq2 and not Path(fastq2).exists():
        raise AlignError(f"mate 2 FASTQ not found: {fastq2}")

    return kind


def assert_gene_info(genome_dir: Path, preset: str) -> None:
    """--quantMode TranscriptomeSAM needs a GTF-built index. Say so before mapping."""
    if not PRESETS[preset]["needs_gene_info"]:
        return
    if not (genome_dir / GENE_INFO).exists():
        raise AlignError(
            f"preset {preset} needs --quantMode TranscriptomeSAM, which requires an index "
            f"built with --sjdbGTFfile. Missing: {genome_dir / GENE_INFO}"
        )


def output_directory(out_prefix: str) -> Path:
    """The directory to create for <out-prefix>.

    <out-prefix> is a prefix, not a path, so the directory is its leading component --
    except when it ends in "/", where STAR takes the whole string as a directory and
    dirname would drop that last component, leaving it uncreated.
    """
    return Path(out_prefix) if out_prefix.endswith("/") else Path(out_prefix).parent


# %%
# Read extraction (BAM input only)


def bamtofastq_argv(bam: str, targets: list[str], paired: bool) -> list[str]:
    """The bamtofastq command line for one BAM.

    collate=1 with F=/F2= for paired data. For single-end, collate=0 and stdout: S= is
    silently ignored under collate=0 -- bamtofastq writes every category to stdout
    regardless. Verified against biobambam2 2.0.183; do not "fix" this to S=. collate=0 is
    deliberate for single-end, which needs no collation pass.
    """
    if paired:
        return [
            "bamtofastq",
            f"filename={bam}",
            "collate=1",
            f"F={targets[0]}",
            f"F2={targets[1]}",
        ]
    return ["bamtofastq", f"filename={bam}", "collate=0"]


def resolve_layout(requested: str, bam: Path | None) -> bool:
    """True when the input is paired. Reads the BAM only when the caller says nothing.

    NOTE: `auto` is the expensive path -- see the --layout note. It is kept as the default
      because a wrong layout is worse than a slow one: it maps every mate as an independent
      read, silently.
    """
    if requested == "paired":
        return True
    if requested == "single":
        return False
    return detect_bam_layout(bam)


def detect_bam_layout(bam: Path) -> bool:
    """True when the BAM holds paired reads, refusing a BAM that holds both.

    NOTE: THIS IS NOT THE CHOKE POINT FOR THE WORKFLOWS. Both of them pass --layout=paired
      from the sample table, so this function is never reached there and nothing on that path
      validates that a BAM really is purely paired. The owner of that invariant is the
      generator of processed/metadata/samples.tsv, which does not exist yet (see the
      gdcdata-prep progress log): it must assert layout purity per BAM, because the
      collate=1 F=/F2= extraction below captures no unpaired category and would drop those
      reads silently -- and assert_reads_were_mapped only refuses ZERO reads, not fewer.
      Until that generator exists, the guarantee rests on GDC's layout field being right.

    samtools is called on its own and its exit code checked, because a failed count must
    not fall through to "single-end" -- that would map every mate as an independent read.
    """
    # NOTE: flagstat, not two `view -c` calls. Both numbers come from ONE decode of the
    #   whole BAM; counting twice doubled the cost of the very thing the --layout option
    #   exists to let a caller skip.
    result = subprocess.run(
        ["samtools", "flagstat", str(bam)], capture_output=True, text=True, check=False
    )
    if result.returncode != 0 or not result.stdout.strip():
        raise AlignError(f"samtools flagstat failed on {bam}: {result.stderr.strip()}")

    counts = {}
    for label, pattern in (("primary", PRIMARY_PATTERN), ("paired", PAIRED_PATTERN)):
        match = pattern.search(result.stdout)
        if match is None:
            raise AlignError(
                f"could not read the {label} count out of samtools flagstat on {bam}"
            )
        counts[label] = int(match.group(1)) + int(match.group(2))

    # NOTE: "any record is paired" is not enough. A BAM holding both kinds would be treated
    #   as paired, and the collate=1 F/F2 extraction captures no unpaired category (no O=,
    #   O2= or S=), so those reads would vanish without a word. Refuse instead.
    if 0 < counts["paired"] < counts["primary"]:
        raise AlignError(
            f"{bam} mixes paired ({counts['paired']}) and unpaired "
            f"({counts['primary'] - counts['paired']}) primary records. The extraction here "
            f"handles one or the other; split the BAM, or pass --layout to state which to "
            f"treat it as."
        )
    return counts["paired"] > 0


def start_extraction(bam: Path, targets: list[str], paired: bool) -> subprocess.Popen:
    """Start bamtofastq writing into `targets`, which may be FIFOs.

    NOTE: for single-end the output is stdout, and stdout has to be redirected into the
      FIFO by a process that is NOT this one -- opening a FIFO for writing blocks until a
      reader appears, and the reader is STAR, which has not started yet. A one-line
      `sh -c` child does the redirection so the block happens over there. The paths are
      shlex-quoted; nothing else reaches the shell.
    """
    argv = bamtofastq_argv(str(bam), targets, paired)
    if paired:
        return subprocess.Popen(argv)
    command = f"exec {shlex.join(argv)} > {shlex.quote(targets[0])}"
    return subprocess.Popen(["sh", "-c", command])


# %%
# Verification


def assert_reads_were_mapped(out_prefix: str, kind: str, read_files_command: list[str]) -> int:
    """The read count out of Log.final.out, refusing zero.

    A non-empty BAM is NOT sufficient. When --readFilesCommand cannot be executed STAR does
    not fail: the reads command dies inside its own subshell, STAR sees an empty stream,
    maps zero reads, writes a small but structurally valid BAM and exits 0. Observed with
    "pigz: command not found" inside the apptainer image, which produced a 30 KB sorted BAM
    and Number of input reads = 0. The read count is the only thing that separates that
    from a real run, and it covers a bamtofastq that produced nothing in the same way.
    """
    log_path = Path(f"{out_prefix}{LOG_FINAL_SUFFIX}")
    try:
        reads = parse_input_reads(log_path.read_text())
    except OSError as error:
        raise AlignError(f"could not read {log_path}: {error}") from error

    if reads is None:
        raise AlignError(f"could not find the input read count in {log_path}")
    if reads == 0:
        detail = (
            f"--readFilesCommand {' '.join(read_files_command[1:])!r} most likely failed; "
            f"check {out_prefix}{STAR_TMP_SUFFIX}/readsCommand_read1"
            if kind == "fastq"
            else "bamtofastq most likely produced nothing"
        )
        raise AlignError(f"STAR processed 0 input reads. {detail}")
    return reads


def assert_outputs_exist(outputs: list[str]) -> None:
    for path in outputs:
        target = Path(path)
        if not target.exists() or target.stat().st_size == 0:
            raise AlignError(f"expected output missing or empty: {path}")


def discard_partial_outputs(out_prefix: str, outputs: list[str]) -> None:
    """A partial SJ.out.tab silently poisons the next genome index, and a structurally
    valid but truncated BAM is worse than no BAM. Remove both rather than leave them for a
    consumer to find."""
    tmp_dir = Path(f"{out_prefix}{STAR_TMP_SUFFIX}")
    if tmp_dir.is_dir():
        subprocess.run(["rm", "-rf", str(tmp_dir)], check=False)
    for path in outputs:
        Path(path).unlink(missing_ok=True)


# %%
# Main


def install_signal_handlers() -> None:
    """Turn the signals that actually kill this into exceptions, so the cleanup runs.

    NOTE: Python's default SIGTERM handling terminates the process WITHOUT unwinding -- no
      finally, no except. Snakemake cancels jobs with SIGTERM, so without this the reaping
      of bamtofastq and the removal of partial output never happen on the most common
      abnormal exit. SIGKILL cannot be caught and is why mkfifo tolerates a leftover FIFO.

    Two windows remain, knowingly: a second signal arriving during the cleanup itself, and
    SIGKILL.
    """
    def raise_system_exit(number, frame):
        raise SystemExit(128 + number)

    for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        # NOTE: never override an inherited SIG_IGN. `nohup snakemake ... &` sets SIGHUP to
        #   ignore and every child inherits that; installing a handler re-enables the signal,
        #   so closing the terminal would kill every in-flight alignment and discard its
        #   output -- the handler would turn a deliberately ignored signal into a job killer.
        #   Same reasoning protects SIGINT for a backgrounded run.
        if signal.getsignal(number) is not signal.SIG_IGN:
            signal.signal(number, raise_system_exit)


def run_watched(
    argv: list[str], extraction: subprocess.Popen | None, input_path: Path
) -> int:
    """Run STAR, watching the extraction it reads from, and return STAR's exit status.

    NOTE: a poll taken just once, before STAR starts, cannot see the failure that matters.
      bamtofastq has to open the BAM, initialise libmaus and parse the header before it can
      die on a corrupt one -- tens of milliseconds -- by which time STAR is already blocked
      in open() on a FIFO nothing will ever write to, and this process is blocked in wait().
      A hang is invisible to a workflow that only watches exit codes, and it holds a slot
      until a person notices. So both processes are polled for as long as STAR runs.
    NOTE: an extraction that exits 0 while STAR is still mapping is the normal end of the
      stream, not a failure. Only a non-zero exit kills STAR here.
    """
    star = subprocess.Popen(argv)
    try:
        while True:
            status = star.poll()
            if status is not None:
                return status

            if extraction is not None:
                extraction_status = extraction.poll()
                if extraction_status is not None and extraction_status != 0:
                    star.kill()
                    star.wait()
                    raise AlignError(
                        f"bamtofastq exited {extraction_status} while STAR was running; "
                        f"STAR was killed rather than left to map a truncated stream or to "
                        f"block on a FIFO that will never be written: {input_path}"
                    )

            time.sleep(WATCHDOG_INTERVAL_SECONDS)
    except BaseException:
        if star.poll() is None:
            star.kill()
            star.wait()
        raise


def align(
    *,
    input_path: Path,
    genome_dir: Path,
    out_prefix: str,
    threads: int,
    fastq2: str | None,
    make_sjdb: bool,
    two_pass: bool,
    preset: str,
    genome_load: str,
    sort_with: str,
    sort_ram: str,
    read_files_command_requested: str,
    layout: str = AUTO_LAYOUT,
) -> int:
    kind = validate(
        input_path=input_path,
        genome_dir=genome_dir,
        fastq2=fastq2,
        read_files_command=read_files_command_requested,
        sort_with=sort_with,
        preset=preset,
        layout=layout,
    )
    if not make_sjdb:
        assert_gene_info(genome_dir, preset)

    genome_load, note = resolve_genome_load(genome_load, two_pass)
    if note:
        print(f"Note    : {note}")

    output_directory(out_prefix).mkdir(parents=True, exist_ok=True)

    # A FIFO can only be read once, so --twopassMode -- which makes STAR read its input
    # twice -- cannot stream. Verified: an identical run hangs on a FIFO and exits 0 on a
    # regular file. It therefore spills real FASTQ next to <out-prefix>, tens of GB per
    # TCGA RNA-seq sample times the batch width. Streaming is the default to avoid that.
    streaming = kind == "bam" and not two_pass

    if kind == "bam":
        paired = resolve_layout(layout, input_path)
        read_files = (
            # Not ".R1.fq": with an out-prefix ending in "/" a leading dot makes the FIFO
            # invisible to a bare ls, which is exactly when someone is looking for a
            # leftover from a killed run.
            [f"{out_prefix}R1.fq", f"{out_prefix}R2.fq"] if paired
            else [f"{out_prefix}reads.fq"]
        )
        temp_read_files = list(read_files)
        read_files_command = []
    else:
        paired = bool(fastq2)
        read_files = [str(input_path)] + ([fastq2] if fastq2 else [])
        # Never the caller's own files: deleting these would destroy the input.
        temp_read_files = []
        read_files_command = resolve_read_files_command(
            input_path, read_files_command_requested
        )

    outputs = expected_outputs(
        out_prefix, make_sjdb, PRESETS[preset]["needs_gene_info"], sort_with
    )

    print(f"Input   : {input_path} ({kind})")
    print(f"Mode    : {'PE' if paired else 'SE'}")
    print(f"Pass    : {'1st (sjdb only)' if make_sjdb else '2nd (aligned BAM)'}")
    print(f"Preset  : {preset}")
    print(f"Load    : {genome_load}")
    if kind == "fastq":
        print(f"Decomp  : {' '.join(read_files_command[1:]) or NO_DECOMPRESSION}")
    else:
        print(f"Stream  : {streaming} (bamtofastq through FIFOs)")
    print(f"Sorter  : {sort_with}" + (
        "  (sort separately: ngsutils sort_star_bam)" if sort_with == NO_SORT else ""
    ))
    print(f"2-pass  : {two_pass} (STAR per-sample --twopassMode)")

    extraction: subprocess.Popen | None = None
    try:
        if kind == "bam":
            if streaming:
                for target in read_files:
                    # A SIGKILLed run leaves the FIFO behind, and mkfifo would then die on
                    # FileExistsError -- a retry must not need a second retry.
                    Path(target).unlink(missing_ok=True)
                    os.mkfifo(target)
            else:
                print("Note    : --two-pass cannot stream; extracting FASTQ to disk first")
            extraction = start_extraction(input_path, read_files, paired)
            if not streaming and extraction.wait() != 0:
                raise AlignError(f"bamtofastq failed: {input_path}")

        argv = build_star_argv(
            genome_dir=str(genome_dir),
            out_prefix=out_prefix,
            read_files=read_files,
            threads=threads,
            genome_load=genome_load,
            preset=preset,
            make_sjdb=make_sjdb,
            two_pass=two_pass,
            sort_with=sort_with,
            sort_ram=sort_ram,
            read_files_command=read_files_command,
        )
        print(f"CMD: {shlex.join(argv)}")
        status = run_watched(argv, extraction if streaming else None, input_path)

        if kind == "bam" and streaming:
            if status != 0:
                # NOTE: poll BEFORE killing. When both die, STAR's death is usually observed
                #   first (the watchdog polls it first) while bamtofastq's own non-zero exit
                #   is the likelier root cause -- a truncated BAM. Killing first replaces that
                #   status with the signal and loses the only diagnosis this tool can offer.
                extraction_status = extraction.poll()
                if extraction_status is None:
                    # STAR may have died before opening the FIFOs, leaving bamtofastq blocked
                    # on open() -- kill it or the wait never returns.
                    extraction.kill()
                    extraction.wait()
                elif extraction_status != 0:
                    raise AlignError(
                        f"bamtofastq exited {extraction_status} and STAR then failed (exit "
                        f"{status}); the extraction is the likelier cause: {input_path}"
                    )
            elif extraction.wait() != 0:
                # A bamtofastq that died mid-stream looks like EOF to STAR, which then
                # exits 0 on truncated input. This check is the only thing separating the
                # two cases.
                raise AlignError(f"bamtofastq failed, alignment is truncated: {input_path}")

        if status != 0:
            raise AlignError(f"STAR failed (exit {status}): {input_path}")

        reads = assert_reads_were_mapped(out_prefix, kind, read_files_command)
        assert_outputs_exist(outputs)

    except BaseException:
        # NOTE: not just AlignError. A missing STAR shim raises FileNotFoundError and a
        #   a signal we installed a handler for raises SystemExit (see
        #   install_signal_handlers), and either would otherwise leave a partial SJ.out.tab
        #   or BAM behind for a consumer to find. Snakemake removes a failed job's DECLARED
        #   outputs, but standalone use has no such backstop.
        discard_partial_outputs(out_prefix, outputs)
        raise
    finally:
        # NOTE: unlinking a FIFO does NOT unblock a writer already sleeping in open() -- the
        #   inode lives on while bamtofastq holds it -- so the extraction has to be reaped
        #   here rather than only on the normal path. Without this, every abnormal exit
        #   leaves one bamtofastq blocked forever.
        if extraction is not None and extraction.poll() is None:
            extraction.kill()
            extraction.wait()
        for target in temp_read_files:
            Path(target).unlink(missing_ok=True)

    print(f"Reads   : {reads} input reads")
    for path in outputs:
        print(f"Output  : {path}")
    return 0


def main() -> int:
    opts = docopt(__doc__)
    install_signal_handlers()
    try:
        return align(
            input_path=Path(opts["<input>"]),
            genome_dir=Path(opts["<genome-dir>"]),
            out_prefix=opts["<out-prefix>"],
            threads=int(opts["<threads>"] or DEFAULT_THREADS),
            fastq2=opts["--fastq2"],
            make_sjdb=opts["--make-sjdb"],
            two_pass=opts["--two-pass"],
            preset=opts["--preset"],
            genome_load=opts["--genome-load"],
            sort_with=opts["--sort-with"],
            sort_ram=opts["--sort-ram"],
            read_files_command_requested=opts["--read-files-command"],
            layout=opts["--layout"],
        )
    except AlignError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
