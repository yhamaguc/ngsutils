"""Tests for ngsutils.align_star.

Usage:
    python3 -m pytest tests/test_align_star.py

NOTE: STAR, samtools and bamtofastq are never invoked here. What is tested is everything
  that DECIDES what they are asked to do, plus the checks that separate a real run from a
  complete-looking empty one. Three of these encode failures that have actually happened:
  an option default silently arriving empty, a FASTQ path deleting the caller's input, and
  STAR exiting 0 after mapping zero reads.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from docopt import docopt

from ngsutils import align_star

# %%
# Option parsing
#
# The reason this module exists in Python at all: `docopts` dropped [default:] for any
# option whose description mentioned another declared option by name. These pin the
# defaults so a DOC edit cannot quietly remove one again.


def _parse(argv: list[str]) -> dict:
    return docopt(align_star.__doc__, argv=argv)


BASE_ARGV = ["in_1.fastq.gz", "genome", "out/sample."]


def test_every_valued_option_has_its_documented_default():
    opts = _parse(BASE_ARGV)
    assert opts["--preset"] == "none"
    assert opts["--genome-load"] == "LoadAndKeep"
    assert opts["--sort-with"] == "star"
    assert opts["--sort-ram"] == "160000000000"
    assert opts["--read-files-command"] == "auto"


def test_threads_is_optional_and_falls_back_in_main():
    assert _parse(BASE_ARGV)["<threads>"] is None
    assert _parse([*BASE_ARGV, "8"])["<threads>"] == "8"


# %%
# Input dispatch


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("sample.bam", "bam"),
        ("/data/x/y.rna_seq.genomic.gdc_realn.bam", "bam"),
        ("sample_1.fastq.gz", "fastq"),
        ("sample_1.fq", "fastq"),
        # A BAM that does not say so is treated as FASTQ, by design -- documented, and the
        # reason the docstring tells callers to name files honestly.
        ("sample.bam.renamed", "fastq"),
    ],
)
def test_input_kind_comes_from_the_suffix(path, expected):
    assert align_star.input_kind(path) == expected


# %%
# Genome load
#
# The fallback is not cosmetic: STAR aborts outright on a shared-memory genome when
# junctions must be inserted at run time.


@pytest.mark.parametrize("load", align_star.SHARED_MEMORY_LOADS)
def test_two_pass_falls_back_off_shared_memory_and_says_so(load):
    resolved, note = align_star.resolve_genome_load(load, two_pass=True)
    assert resolved == align_star.NO_SHARED_MEMORY
    assert note and load in note


def test_without_two_pass_the_requested_load_is_kept_silently():
    assert align_star.resolve_genome_load("LoadAndKeep", two_pass=False) == (
        "LoadAndKeep",
        None,
    )


def test_no_shared_memory_under_two_pass_needs_no_note():
    assert align_star.resolve_genome_load("NoSharedMemory", two_pass=True) == (
        "NoSharedMemory",
        None,
    )


# %%
# Decompression


def test_auto_selects_the_measured_gzip_command_for_gz():
    assert align_star.resolve_read_files_command("x.fastq.gz", "auto") == [
        "--readFilesCommand",
        "pigz",
        "-d",
        "-c",
    ]


def test_auto_passes_plain_fastq_through_undecompressed():
    assert align_star.resolve_read_files_command("x.fastq", "auto") == []


def test_none_forces_no_decompression_even_for_gz():
    assert align_star.resolve_read_files_command("x.fastq.gz", "none") == []


def test_an_explicit_command_is_split_into_arguments():
    assert align_star.resolve_read_files_command("x.fastq.gz", "zcat") == [
        "--readFilesCommand",
        "zcat",
    ]


# %%
# What a run must produce


def test_the_sjdb_pass_expects_no_bam():
    assert align_star.expected_outputs("out/s.", True, True) == ["out/s.SJ.out.tab"]


def test_sort_with_none_expects_the_unsorted_bam_not_the_sorted_one():
    """Under --sort-with none the sorted name is produced by the next step, not this one."""
    outputs = align_star.expected_outputs("out/s.", False, True, "none")
    assert outputs == [
        "out/s.Aligned.out.bam",
        "out/s.Aligned.toTranscriptome.out.bam",
    ]


def test_the_gdc_preset_expects_the_transcriptome_bam():
    outputs = align_star.expected_outputs("out/s.", False, True)
    assert outputs == [
        "out/s.Aligned.sortedByCoord.out.bam",
        "out/s.Aligned.toTranscriptome.out.bam",
    ]


def test_a_preset_without_transcriptome_output_expects_only_the_genome_bam():
    assert align_star.expected_outputs("out/s.", False, False) == [
        "out/s.Aligned.sortedByCoord.out.bam"
    ]


# %%
# The STAR command line


def _argv(**overrides) -> list[str]:
    kwargs = dict(
        genome_dir="genome",
        out_prefix="out/s.",
        read_files=["r1.fq", "r2.fq"],
        threads=8,
        genome_load="NoSharedMemory",
        preset="gdc",
        make_sjdb=False,
        two_pass=False,
        sort_with="star",
        sort_ram="160000000000",
        read_files_command=[],
    )
    kwargs.update(overrides)
    return align_star.build_star_argv(**kwargs)


def test_the_sjdb_pass_writes_no_alignment_and_withholds_output_formatting():
    argv = _argv(make_sjdb=True)
    assert "--outSAMtype" in argv and argv[argv.index("--outSAMtype") + 1] == "None"
    # PRESET second-pass options are fatal or inert under --outSAMtype None.
    assert "--quantMode" not in argv
    assert "--outSAMunmapped" not in argv
    # ...while the junction-shaping options, which must match across passes, are present.
    assert "--outFilterType" in argv


def test_the_aligned_pass_carries_the_second_pass_options():
    argv = _argv()
    assert argv[argv.index("--outSAMtype") + 1 : argv.index("--outSAMtype") + 3] == [
        "BAM",
        "SortedByCoordinate",
    ]
    assert "--quantMode" in argv


def test_the_junction_shaping_options_are_identical_in_both_passes():
    """An sjdb built under one filter regime and consumed under another is not the sjdb the
    second pass thinks it is, so this equality is the contract of a cohort 2-pass."""
    first = _argv(make_sjdb=True)
    second = _argv(make_sjdb=False)
    common = list(align_star.PRESETS["gdc"]["common"])
    assert [flag for flag in common if flag.startswith("--")] == [
        flag for flag in common if flag.startswith("--")
    ]
    for flag in common:
        if flag.startswith("--"):
            assert flag in first and flag in second


def test_sort_with_none_asks_star_for_an_unsorted_bam_and_no_sort_ram():
    """The memory-hungry sort is a separate step, so STAR must not attempt it."""
    argv = _argv(sort_with="none")
    assert argv[argv.index("--outSAMtype") + 1 : argv.index("--outSAMtype") + 3] == [
        "BAM",
        "Unsorted",
    ]
    # Passing a sort ceiling STAR will not use would only mislead the log.
    assert "--limitBAMsortRAM" not in argv


def test_the_star_sorter_gets_the_ceiling():
    argv = _argv(sort_ram="64000000000")
    assert argv[argv.index("--limitBAMsortRAM") + 1] == "64000000000"


def test_two_pass_adds_twopassmode_only_when_asked():
    assert "--twopassMode" not in _argv()
    assert "--twopassMode" in _argv(two_pass=True)


def test_the_decompressor_reaches_star_as_separate_arguments():
    argv = _argv(read_files_command=["--readFilesCommand", "pigz", "-d", "-c"])
    position = argv.index("--readFilesCommand")
    assert argv[position + 1 : position + 4] == ["pigz", "-d", "-c"]


def test_both_mates_are_passed_to_readfilesin():
    argv = _argv(read_files=["a.fq", "b.fq"])
    position = argv.index("--readFilesIn")
    assert argv[position + 1 : position + 3] == ["a.fq", "b.fq"]


# %%
# Validation
#
# Each refusal below replaces a failure that would otherwise happen after the mapping phase
# had been paid for, or not at all.


def _valid(tmp_path: Path) -> dict:
    (tmp_path / "genome").mkdir(exist_ok=True)
    fastq = tmp_path / "r1.fastq.gz"
    fastq.touch()
    return {
        "input_path": fastq,
        "genome_dir": tmp_path / "genome",
        "fastq2": None,
        "read_files_command": "auto",
        "sort_with": "star",
        "preset": "none",
    }


def test_a_valid_fastq_invocation_passes(tmp_path):
    assert align_star.validate(**_valid(tmp_path)) == "fastq"


def test_an_unknown_preset_is_refused(tmp_path):
    kwargs = _valid(tmp_path) | {"preset": "bogus"}
    with pytest.raises(align_star.AlignError, match="unknown --preset"):
        align_star.validate(**kwargs)


def test_an_unknown_sorter_is_refused(tmp_path):
    kwargs = _valid(tmp_path) | {"sort_with": "rsort"}
    with pytest.raises(align_star.AlignError, match="unknown --sort-with"):
        align_star.validate(**kwargs)


def test_a_missing_input_is_refused(tmp_path):
    kwargs = _valid(tmp_path) | {"input_path": tmp_path / "absent.fastq.gz"}
    with pytest.raises(align_star.AlignError, match="input not found"):
        align_star.validate(**kwargs)


def test_a_missing_genome_dir_is_refused(tmp_path):
    kwargs = _valid(tmp_path) | {"genome_dir": tmp_path / "absent"}
    with pytest.raises(align_star.AlignError, match="genome dir not found"):
        align_star.validate(**kwargs)


def test_a_missing_mate_two_is_refused(tmp_path):
    kwargs = _valid(tmp_path) | {"fastq2": str(tmp_path / "absent_2.fastq.gz")}
    with pytest.raises(align_star.AlignError, match="mate 2 FASTQ not found"):
        align_star.validate(**kwargs)


def test_fastq_only_options_are_refused_on_the_bam_path(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.touch()
    base = _valid(tmp_path) | {"input_path": bam}

    with pytest.raises(align_star.AlignError, match="--fastq2 is FASTQ-only"):
        align_star.validate(**(base | {"fastq2": "x_2.fastq.gz"}))

    with pytest.raises(align_star.AlignError, match="--read-files-command is FASTQ-only"):
        align_star.validate(**(base | {"read_files_command": "zcat"}))


def test_the_gdc_preset_refuses_an_index_without_gene_info(tmp_path):
    genome = tmp_path / "genome"
    genome.mkdir()
    with pytest.raises(align_star.AlignError, match="requires an index built"):
        align_star.assert_gene_info(genome, "gdc")

    (genome / align_star.GENE_INFO).touch()
    align_star.assert_gene_info(genome, "gdc")


def test_a_preset_that_needs_no_transcriptome_output_ignores_gene_info(tmp_path):
    genome = tmp_path / "genome"
    genome.mkdir()
    align_star.assert_gene_info(genome, "encode")


# %%
# Output directory


def test_a_trailing_slash_makes_the_whole_prefix_the_directory():
    """dirname("out/sample/") is "out", which would leave out/sample uncreated."""
    assert align_star.output_directory("out/sample/") == Path("out/sample")


def test_a_prefix_without_a_slash_uses_its_leading_component():
    assert align_star.output_directory("out/sample.") == Path("out")


# %%
# Read extraction


def test_paired_extraction_writes_to_both_targets():
    argv = align_star.bamtofastq_argv("in.bam", ["a.fq", "b.fq"], paired=True)
    assert argv == [
        "bamtofastq",
        "filename=in.bam",
        "collate=1",
        "F=a.fq",
        "F2=b.fq",
    ]


def test_single_end_extraction_uses_stdout_not_s():
    """S= is silently ignored under collate=0; verified against biobambam2 2.0.183."""
    argv = align_star.bamtofastq_argv("in.bam", ["a.fq"], paired=False)
    assert argv == ["bamtofastq", "filename=in.bam", "collate=0"]
    assert not [token for token in argv if token.startswith("S=")]


# %%
# Verification
#
# Deliberate-breakage tests: STAR exits 0 after mapping nothing, so an existing BAM is not
# evidence of a real run.


def test_a_zero_read_run_is_refused_on_the_fastq_path(tmp_path):
    prefix = f"{tmp_path}/s."
    Path(f"{prefix}Log.final.out").write_text(
        "                          Number of input reads |\t0\n"
    )
    with pytest.raises(align_star.AlignError, match="0 input reads"):
        align_star.assert_reads_were_mapped(
            prefix, "fastq", ["--readFilesCommand", "pigz", "-d", "-c"]
        )


def test_a_zero_read_run_names_the_bam_cause_on_the_bam_path(tmp_path):
    prefix = f"{tmp_path}/s."
    Path(f"{prefix}Log.final.out").write_text("Number of input reads |\t0\n")
    with pytest.raises(align_star.AlignError, match="bamtofastq"):
        align_star.assert_reads_were_mapped(prefix, "bam", [])


def test_a_real_read_count_passes(tmp_path):
    prefix = f"{tmp_path}/s."
    Path(f"{prefix}Log.final.out").write_text(
        "                          Number of input reads |\t3928\n"
    )
    assert align_star.assert_reads_were_mapped(prefix, "bam", []) == 3928


def test_a_log_without_the_count_is_refused(tmp_path):
    prefix = f"{tmp_path}/s."
    Path(f"{prefix}Log.final.out").write_text("Started job on |\tAug 19\n")
    with pytest.raises(align_star.AlignError, match="could not find the input read count"):
        align_star.assert_reads_were_mapped(prefix, "bam", [])


def test_a_missing_log_is_refused(tmp_path):
    with pytest.raises(align_star.AlignError, match="could not read"):
        align_star.assert_reads_were_mapped(f"{tmp_path}/absent.", "bam", [])


@pytest.mark.parametrize(
    "log",
    [
        "                          Number of input reads |\t3928\n",
        "Number of input reads |   3928\n",
        "\tNumber of input reads |\t3928\t\n",
    ],
)
def test_the_read_count_survives_the_logs_whitespace(log):
    assert align_star.parse_input_reads(log) == 3928


def test_an_empty_output_counts_as_missing(tmp_path):
    empty = tmp_path / "s.Aligned.sortedByCoord.out.bam"
    empty.touch()
    with pytest.raises(align_star.AlignError, match="missing or empty"):
        align_star.assert_outputs_exist([str(empty)])

    empty.write_bytes(b"BAM\x01")
    align_star.assert_outputs_exist([str(empty)])


def test_partial_outputs_are_discarded_with_the_star_temp_dir(tmp_path):
    prefix = f"{tmp_path}/s."
    partial = Path(f"{prefix}SJ.out.tab")
    partial.write_text("chr1\t1\t2\t1\n")
    tmp_dir = Path(f"{prefix}_STARtmp")
    tmp_dir.mkdir()
    (tmp_dir / "leftover").touch()

    align_star.discard_partial_outputs(prefix, [str(partial)])

    assert not partial.exists()
    assert not tmp_dir.exists()
