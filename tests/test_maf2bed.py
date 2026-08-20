"""Tests for ngsutils.maf2bed.

The conversion this command performs is the one named in the coordinate trap: MAF is 1-based
inclusive, BED is 0-based half-open. An off-by-one here is silent -- every downstream
intersection still runs and every answer shifts by a base -- so the coordinate tests assert
exact numbers rather than shapes.
"""

from __future__ import annotations

import io

import pytest
from docopt import docopt

from ngsutils import maf2bed


# %%
# Fixtures
#
# A MAF row with the four positional columns where the specification puts them: 5 Chromosome,
# 6 Start_Position, 7 End_Position, 8 Strand. The others carry recognisable values so the
# name field can be checked for order.


def _row(chromosome="chr7", start="140753336", end="140753336", strand="+", extra=None):
    before = ["BRAF", "673", "hgnc", "GRCh38"]
    after = extra if extra is not None else ["SNP", "A", "T"]
    return "\t".join([*before, chromosome, start, end, strand, *after])


def _convert(line, slop=0, line_number=1):
    return maf2bed.convert_line(line, line_number, slop)


# %%
# Coordinates


def test_a_point_mutation_becomes_a_one_base_interval():
    fields = _convert(_row(start="140753336", end="140753336"))
    assert fields[0] == "chr7"
    # 1-based inclusive [336, 336] is 0-based half-open [335, 336): one base.
    assert fields[1] == "140753335"
    assert fields[2] == "140753336"
    assert int(fields[2]) - int(fields[1]) == 1


def test_a_span_keeps_its_length():
    fields = _convert(_row(start="100", end="109"))
    assert (fields[1], fields[2]) == ("99", "109")
    assert int(fields[2]) - int(fields[1]) == 10


def test_slop_pads_both_edges():
    fields = _convert(_row(start="1000", end="1000"), slop=10)
    assert (fields[1], fields[2]) == ("989", "1010")
    assert int(fields[2]) - int(fields[1]) == 21


def test_slop_larger_than_the_start_clamps_at_zero():
    # Without the clamp this is -4, which is not a legal BED start.
    fields = _convert(_row(start="5", end="5"), slop=10)
    assert fields[1] == "0"
    assert fields[2] == "15"


def test_the_upper_edge_is_deliberately_not_clamped():
    # A MAF carries no contig lengths, so there is nothing to clamp against. Documented
    # rather than fixed; this test exists so the behaviour cannot change unnoticed.
    fields = _convert(_row(start="1000", end="1000"), slop=1000000)
    assert fields[2] == "1001000"


# %%
# BED6 shape


def test_the_output_is_bed6():
    fields = _convert(_row())
    assert len(fields) == 6


def test_the_score_column_is_a_placeholder_zero():
    assert _convert(_row())[4] == "0"


def test_the_strand_is_carried_through():
    assert _convert(_row(strand="-"))[5] == "-"
    assert _convert(_row(strand="."))[5] == "."


# %%
# The name field


def test_every_non_positional_column_survives_in_order():
    fields = _convert(_row(extra=["SNP", "A", "T"]))
    assert fields[3] == '"BRAF";"673";"hgnc";"GRCh38";"SNP";"A";"T"'


def test_na_becomes_empty_rather_than_the_string_na():
    fields = _convert(_row(extra=["NA", "A"]))
    assert fields[3].endswith(';;"A"')


def test_values_are_joined_with_semicolons_not_commas():
    # MAF free-text columns contain commas, which is why the separator is not one.
    fields = _convert(_row(extra=["one, two", "three"]))
    assert fields[3].count(";") == 5
    assert '"one, two"' in fields[3]


# %%
# What is skipped


@pytest.mark.parametrize("line", ["#version 2.4\n", "Hugo_Symbol\tEntrez\n", "individual\tx\n"])
def test_comments_and_headers_are_skipped(line):
    assert _convert(line) is None


def test_a_blank_line_is_skipped():
    assert _convert("\n") is None
    assert _convert("   \n") is None


# %%
# What is refused
#
# Each of these would otherwise emit a plausible BED line, which is the failure mode worth
# spending a test on.


def test_a_row_too_short_to_reach_strand_is_a_hard_error():
    with pytest.raises(SystemExit):
        _convert("a\tb\tc\td\tchr1\t100\t100\n")


def test_non_integer_coordinates_are_a_hard_error():
    with pytest.raises(SystemExit):
        _convert(_row(start="not-a-number"))


def test_a_strand_bed_does_not_accept_is_a_hard_error():
    # The specification says Strand is always '+', so this is a reheadered or edited file.
    with pytest.raises(SystemExit):
        _convert(_row(strand="1"))


def test_a_negative_slop_is_refused(monkeypatch):
    monkeypatch.setattr("sys.argv", ["maf2bed", "--slop", "-1"])
    monkeypatch.setattr("sys.stdin", io.StringIO(""))
    with pytest.raises(SystemExit):
        maf2bed.main()


def test_a_non_integer_slop_is_refused(monkeypatch):
    monkeypatch.setattr("sys.argv", ["maf2bed", "--slop", "ten"])
    monkeypatch.setattr("sys.stdin", io.StringIO(""))
    with pytest.raises(SystemExit):
        maf2bed.main()


# %%
# End to end


def test_main_writes_one_bed_line_per_data_row(monkeypatch, capsys):
    maf = "#version 2.4\n" + "Hugo_Symbol\tx\n" + _row() + "\n" + _row(chromosome="chrX") + "\n"
    monkeypatch.setattr("sys.argv", ["maf2bed"])
    monkeypatch.setattr("sys.stdin", io.StringIO(maf))

    assert maf2bed.main() == 0

    lines = capsys.readouterr().out.strip().split("\n")
    assert len(lines) == 2
    assert [line.split("\t")[0] for line in lines] == ["chr7", "chrX"]
    assert all(len(line.split("\t")) == 6 for line in lines)


# %%
# The docopt pins -- see the docopt entry in ~/.claude/CLAUDE.md for why these exist


def test_the_docstring_declares_exactly_the_intended_options():
    parsed = docopt(maf2bed.__doc__, argv=[])
    assert sorted(key for key in parsed if key.startswith("-")) == ["--help", "--slop"]


def test_the_slop_default_survives_parsing():
    assert docopt(maf2bed.__doc__, argv=[])["--slop"] == "0"


def test_no_prose_line_begins_with_a_dash():
    """A line starting with '-' anywhere in the docstring becomes an option definition.

    This file had one -- a paragraph opening with the option's own name -- and it produced a
    junk '--' key. Measured 2026-08-20.
    """
    offenders = [
        line
        for line in maf2bed.__doc__.splitlines()
        if line.strip().startswith("-") and "  " not in line.strip()[:24]
    ]
    assert offenders == [], f"these lines will parse as options: {offenders}"


def test_no_option_description_names_another_declared_option():
    declared = {"--slop", "--help"}
    for line in maf2bed.__doc__.splitlines():
        stripped = line.strip()
        if not stripped.startswith("-"):
            continue
        name = stripped.split()[0].split("=")[0].split(",")[0]
        description = stripped[len(stripped.split("  ")[0]):]
        for other in declared - {name}:
            assert other not in description, f"{name}'s description names {other}"
