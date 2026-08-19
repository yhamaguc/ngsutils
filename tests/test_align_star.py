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


def test_the_docstring_declares_exactly_the_intended_options():
    """docopt reads ANY line whose first non-space character is "-" as an option definition,
    so a section underline of dashes or a prose line beginning with a flag name silently
    becomes an option. Measured 2026-08-20: the docstring had produced the keys "--",
    "---------------------------------------------" and "--preset:" this way.
    """
    # NOTE: startswith("-"), not "--". A continuation line beginning with a SHORT dash token
    #   -- "-T", "-m", both plausible in prose about samtools -- creates a single-dash option
    #   that a "--" filter would never see.
    declared = {key for key in _parse(BASE_ARGV) if key.startswith("-")}
    assert declared == {
        "--fastq2",
        "--genome-load",
        "--help",
        "--layout",
        "--make-sjdb",
        "--preset",
        "--read-files-command",
        "--sort-ram",
        "--sort-with",
        "--two-pass",
    }  # -h folds into --help; a stray short option would show up as its own key


def test_no_option_description_names_another_declared_option():
    """The failure this prevents is silent: naming a declared option inside another option's
    DESCRIPTION drops the latter's [default:]. It cost --sort-with and --sort-ram under
    docopts, and --layout under Python docopt before this test existed.
    """
    declared = {
        key for key in _parse(BASE_ARGV) if key.startswith("-")
    } - {"--help", "-h"}
    doc = align_star.__doc__
    options = doc[doc.index("Options:"):]
    offenders = []
    for line in options.splitlines():
        stripped = line.strip()
        own = stripped.split()[0].split("=")[0] if stripped.startswith("-") else None
        for flag in declared:
            if flag != own and flag in stripped:
                offenders.append((own or "(continuation)", flag, stripped[:60]))
    assert not offenders, f"option descriptions naming other options: {offenders}"


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


def _flag_values(argv: list[str]) -> dict[str, list[str]]:
    """{flag: its values} out of a STAR command line."""
    parsed, current = {}, None
    for token in argv:
        if token.startswith("--"):
            current = token
            parsed[current] = []
        elif current is not None:
            parsed[current].append(token)
    return parsed


def test_the_junction_shaping_options_are_identical_in_both_passes():
    """An sjdb built under one filter regime and consumed under another is not the sjdb the
    second pass thinks it is, so this equality is the contract of a cohort 2-pass.

    The earlier version of this test compared a list comprehension to a copy of itself and
    then checked only that each flag NAME appeared in both -- it would have passed while the
    two passes disagreed on every value.
    """
    first = _flag_values(_argv(make_sjdb=True))
    second = _flag_values(_argv(make_sjdb=False))
    common = _flag_values(["STAR", *align_star.PRESETS["gdc"]["common"]])

    for flag, values in common.items():
        assert first.get(flag) == values, f"{flag} differs in the sjdb pass"
        assert second.get(flag) == values, f"{flag} differs in the aligned pass"


def test_the_sjdb_pass_carries_no_second_pass_flag():
    """--outSAMtype None makes the output-formatting options inert or fatal (exit 102)."""
    first = _flag_values(_argv(make_sjdb=True))
    second_only = _flag_values(["STAR", *align_star.PRESETS["gdc"]["second"]])
    assert not (set(first) & set(second_only))


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


# %%
# Layout
#
# `auto` costs a full decode of the BAM, twice per sample across the two passes. A caller
# that already knows the layout -- a sample table carries it -- must be able to say so, and
# a wrong value must be refused rather than silently mapping mates as independent reads.


def test_the_layout_default_is_auto():
    assert _parse(BASE_ARGV)["--layout"] == "auto"


def test_an_explicit_layout_skips_reading_the_bam():
    """None as the BAM proves nothing was read: detect_bam_layout would raise on it."""
    assert align_star.resolve_layout("paired", None) is True
    assert align_star.resolve_layout("single", None) is False


def test_an_unknown_layout_is_refused(tmp_path):
    kwargs = _valid(tmp_path) | {"layout": "pared"}
    with pytest.raises(align_star.AlignError, match="unknown --layout"):
        align_star.validate(**kwargs)


def test_the_fifos_are_not_hidden_files():
    """With an out-prefix ending in "/", a leading dot hides the FIFO from a bare ls -- which
    is when someone is hunting a leftover from a killed run."""
    source = (
        Path(align_star.__file__).read_text()
        if hasattr(align_star, "__file__")
        else ""
    )
    assert '{out_prefix}.R1.fq' not in source
    assert '{out_prefix}R1.fq' in source


# %%
# Mixed layouts
#
# NOTE: the collate=1 F/F2 extraction captures no unpaired category (no O=, O2=, S=), so a
# BAM holding both kinds would lose its unpaired reads silently if "any record is paired"
# decided the layout. Refusing is the only safe answer, and --layout is the way past it.


def _flagstat(primary: int, paired: int, secondary: int = 0, qc_failed: int = 0) -> str:
    """samtools flagstat output, with the real field set and ordering.

    NOTE: `in total` includes secondary alignments while `paired in sequencing` counts only
    primary records -- measured on a real TCGA-UVM BAM as 40,000 total = 15,295 primary +
    24,705 secondary against 15,295 paired. The secondary argument here exists so the tests
    would fail if the code ever compared paired against `in total` again.
    """
    return (
        f"{primary + secondary} + {qc_failed} in total (QC-passed reads + QC-failed reads)\n"
        f"{primary} + {qc_failed} primary\n"
        f"{secondary} + 0 secondary\n"
        f"0 + 0 supplementary\n"
        f"{paired} + {qc_failed} paired in sequencing\n"
    )


class _Ran:
    """A stand-in for subprocess.run's result."""

    def __init__(self, stdout="", returncode=0, stderr=""):
        self.stdout, self.returncode, self.stderr = stdout, returncode, stderr


def test_a_bam_mixing_paired_and_unpaired_records_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(
        align_star.subprocess, "run", lambda *a, **k: _Ran(_flagstat(100, 40))
    )
    with pytest.raises(align_star.AlignError, match="mixes paired"):
        align_star.detect_bam_layout(tmp_path / "mixed.bam")


def test_a_fully_paired_bam_is_paired(tmp_path, monkeypatch):
    monkeypatch.setattr(
        align_star.subprocess, "run", lambda *a, **k: _Ran(_flagstat(100, 100))
    )
    assert align_star.detect_bam_layout(tmp_path / "pe.bam") is True


def test_multimapping_does_not_make_a_paired_bam_look_mixed(tmp_path, monkeypatch):
    """The regression this pins: comparing `paired in sequencing` against `in total` reports
    every multimapping paired BAM as a mix, and would have refused ordinary GDC input. Real
    numbers from a TCGA-UVM BAM: 15,295 primary, 24,705 secondary, 15,295 paired."""
    monkeypatch.setattr(
        align_star.subprocess,
        "run",
        lambda *a, **k: _Ran(_flagstat(15295, 15295, secondary=24705)),
    )
    assert align_star.detect_bam_layout(tmp_path / "multimapped.bam") is True


def test_a_fully_unpaired_bam_is_single(tmp_path, monkeypatch):
    monkeypatch.setattr(
        align_star.subprocess, "run", lambda *a, **k: _Ran(_flagstat(100, 0))
    )
    assert align_star.detect_bam_layout(tmp_path / "se.bam") is False


def test_an_empty_bam_reads_as_single_end(tmp_path, monkeypatch):
    """0 paired of 0 total is not a mix. It reaches STAR and assert_reads_were_mapped
    refuses it there -- pinned deliberately so the fall-through stays intentional."""
    monkeypatch.setattr(
        align_star.subprocess, "run", lambda *a, **k: _Ran(_flagstat(0, 0))
    )
    assert align_star.detect_bam_layout(tmp_path / "empty.bam") is False


def test_a_failing_samtools_is_refused_rather_than_read_as_single(tmp_path, monkeypatch):
    """A failed count must not fall through to single-end: that would map every mate as an
    independent read."""
    monkeypatch.setattr(
        align_star.subprocess,
        "run",
        lambda *a, **k: _Ran("", returncode=1, stderr="truncated file"),
    )
    with pytest.raises(align_star.AlignError, match="flagstat failed"):
        align_star.detect_bam_layout(tmp_path / "broken.bam")


def test_flagstat_without_the_expected_lines_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(
        align_star.subprocess, "run", lambda *a, **k: _Ran("unexpected output\n")
    )
    with pytest.raises(align_star.AlignError, match="could not read"):
        align_star.detect_bam_layout(tmp_path / "odd.bam")


# %%
# The watchdog
#
# NOTE: this is the fix for a defect that a single pre-flight poll could not catch. A
# bamtofastq dying on a corrupt header needs tens of milliseconds to get there, by which
# point STAR is already blocked in open() on a FIFO nothing will write to -- and a hang is
# invisible to a workflow that watches exit codes.


class _FakeProcess:
    """A Popen stand-in whose exit status appears after a given number of polls."""

    def __init__(self, returncode=0, polls_until_exit=0):
        self._returncode = returncode
        self._remaining = polls_until_exit
        self.killed = False

    def poll(self):
        if self._remaining > 0:
            self._remaining -= 1
            return None
        return self._returncode

    def kill(self):
        self.killed = True

    def wait(self):
        return self._returncode


def test_a_dead_extraction_kills_star_instead_of_letting_it_hang(monkeypatch, tmp_path):
    star = _FakeProcess(returncode=0, polls_until_exit=99)   # STAR would run forever
    extraction = _FakeProcess(returncode=1, polls_until_exit=0)  # died already
    monkeypatch.setattr(align_star.subprocess, "Popen", lambda *a, **k: star)
    monkeypatch.setattr(align_star.time, "sleep", lambda seconds: None)

    with pytest.raises(align_star.AlignError, match="bamtofastq exited 1"):
        align_star.run_watched(["STAR"], extraction, tmp_path / "in.bam")
    assert star.killed, "STAR must be killed, not left blocked on the FIFO"


def test_an_extraction_that_finished_cleanly_does_not_kill_star(monkeypatch, tmp_path):
    """Exit 0 while STAR still maps is the normal end of the stream."""
    star = _FakeProcess(returncode=0, polls_until_exit=3)
    extraction = _FakeProcess(returncode=0, polls_until_exit=0)
    monkeypatch.setattr(align_star.subprocess, "Popen", lambda *a, **k: star)
    monkeypatch.setattr(align_star.time, "sleep", lambda seconds: None)

    assert align_star.run_watched(["STAR"], extraction, tmp_path / "in.bam") == 0
    assert not star.killed


def test_stars_own_exit_status_is_returned(monkeypatch, tmp_path):
    star = _FakeProcess(returncode=102, polls_until_exit=0)
    monkeypatch.setattr(align_star.subprocess, "Popen", lambda *a, **k: star)
    monkeypatch.setattr(align_star.time, "sleep", lambda seconds: None)
    assert align_star.run_watched(["STAR"], None, tmp_path / "in.bam") == 102


def test_an_interruption_kills_a_running_star(monkeypatch, tmp_path):
    star = _FakeProcess(returncode=0, polls_until_exit=99)
    monkeypatch.setattr(align_star.subprocess, "Popen", lambda *a, **k: star)

    def interrupt(seconds):
        raise KeyboardInterrupt

    monkeypatch.setattr(align_star.time, "sleep", interrupt)
    with pytest.raises(KeyboardInterrupt):
        align_star.run_watched(["STAR"], None, tmp_path / "in.bam")
    assert star.killed


def test_qc_failed_records_are_counted_because_bamtofastq_extracts_them(tmp_path, monkeypatch):
    """A BAM whose paired records are all QC-failed reads as single-end if only the QC-passed
    field is read -- and then its mates are extracted and mapped as independent reads."""
    monkeypatch.setattr(
        align_star.subprocess,
        "run",
        lambda *a, **k: _Ran(_flagstat(0, 0, qc_failed=100)),
    )
    assert align_star.detect_bam_layout(tmp_path / "qcfail.bam") is True
