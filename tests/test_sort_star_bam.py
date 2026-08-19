"""Tests for ngsutils.sort_star_bam.

Usage:
    python3 -m pytest tests/test_sort_star_bam.py

NOTE: samtools is never invoked. What is tested is the refusals and the path arithmetic --
  the record-count comparison itself needs a real BAM and is covered by the end-to-end run
  recorded in the module docstring.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from docopt import docopt

from ngsutils import sort_star_bam


def test_the_documented_defaults_are_what_docopt_gives():
    opts = docopt(sort_star_bam.__doc__, argv=["in.bam", "out.bam"])
    assert opts["--threads"] == "8"
    assert opts["--memory-per-thread"] == "4G"
    assert opts["--keep-unsorted"] is False
    assert opts["--no-index"] is False


def test_the_temp_prefix_sits_next_to_the_output():
    """-T must be on a filesystem with room for the BAM, not in a small /tmp."""
    assert sort_star_bam.temp_prefix("/data/x/Aligned.sortedByCoord.out.bam") == (
        "/data/x/Aligned.sortedByCoord.out.sorttmp"
    )


def test_a_target_without_the_bam_suffix_still_gets_a_prefix():
    assert sort_star_bam.temp_prefix("/data/x/out") == "/data/x/out.sorttmp"


def test_a_missing_input_is_refused(tmp_path):
    with pytest.raises(sort_star_bam.SortError, match="missing or empty"):
        sort_star_bam.sort_bam(str(tmp_path / "absent.bam"), str(tmp_path / "out.bam"))


def test_an_empty_input_is_refused(tmp_path):
    """An empty unsorted BAM means the alignment did not produce one; sorting it would
    produce a valid-looking empty result."""
    empty = tmp_path / "Aligned.out.bam"
    empty.touch()
    with pytest.raises(sort_star_bam.SortError, match="missing or empty"):
        sort_star_bam.sort_bam(str(empty), str(tmp_path / "out.bam"))


def test_leftover_temporaries_are_cleared(tmp_path):
    """Fragments from a died sort must not be read by a retry."""
    prefix = str(tmp_path / "Aligned.sortedByCoord.out.sorttmp")
    for n in range(3):
        Path(f"{prefix}.{n:04d}.bam").touch()
    unrelated = tmp_path / "keep.bam"
    unrelated.touch()

    sort_star_bam.clear_temporaries(prefix)

    assert not list(tmp_path.glob("*.sorttmp.*.bam"))
    assert unrelated.exists()


# %%
# Output naming
#
# The sample lives in the directory, so the derived name must keep the directory AND the
# prefix and change only the part that states how the file is ordered.


def test_stars_own_name_is_reproduced_exactly():
    """So nothing downstream can tell which sorter ran."""
    assert sort_star_bam.derive_sorted_name(
        "results/align_star/TCGA-UVM/f1/a1/Aligned.out.bam"
    ) == "results/align_star/TCGA-UVM/f1/a1/Aligned.sortedByCoord.out.bam"


def test_a_non_empty_star_prefix_is_kept():
    assert sort_star_bam.derive_sorted_name(
        "out/sample.Aligned.out.bam"
    ) == "out/sample.Aligned.sortedByCoord.out.bam"


def test_a_plain_bam_gets_the_infix():
    assert sort_star_bam.derive_sorted_name("x/sample.bam") == "x/sample.sortedByCoord.bam"


def test_a_name_without_a_bam_suffix_still_gets_one():
    assert sort_star_bam.derive_sorted_name("x/sample") == "x/sample.sortedByCoord.bam"


def test_the_derived_name_is_never_the_input(tmp_path):
    """A derivation that returned the input would have samtools sort read and write the
    same file. Refused rather than attempted."""
    same = tmp_path / "x.sortedByCoord.bam"
    same.write_bytes(b"BAM\x01")
    with pytest.raises(sort_star_bam.SortError, match="would overwrite the input"):
        sort_star_bam.sort_bam(str(same), str(same))


def test_a_bare_number_cannot_be_mistaken_for_the_output_name():
    """Regression: with threads as a trailing positional, `sort_star_bam in.bam 8` bound 8
    to <sorted-bam> and wrote a file called "8". Threads is an option now, so the same
    command line names the output explicitly or not at all."""
    opts = docopt(sort_star_bam.__doc__, argv=["in.bam", "--threads", "8"])
    assert opts["<sorted-bam>"] is None
    assert opts["--threads"] == "8"


def test_an_output_name_that_is_not_a_bam_is_refused(tmp_path):
    """The destructive one: `sort_star_bam in.bam 8` binds 8 to <sorted-bam>, and this
    command deletes its input after a verified sort -- so without this guard the unsorted
    BAM is gone and the result is a file named "8". Seen twice on real data."""
    unsorted = tmp_path / "Aligned.out.bam"
    unsorted.write_bytes(b"BAM\x01")
    with pytest.raises(sort_star_bam.SortError, match="does not end in .bam"):
        sort_star_bam.sort_bam(str(unsorted), "8")
    assert unsorted.exists(), "the input must survive a refused run"
