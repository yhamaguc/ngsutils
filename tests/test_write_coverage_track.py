#!/usr/bin/env python3
"""test_write_coverage_track — the parts that fail plausibly rather than loudly.

Two failure modes are worth a test each, and neither raises on its own:

  **A label that collapses.** Two BAMs reduced to one label give one bedGraph where two were
  expected, with every number in it correct — twenty samples in, fewer files out. This is what
  `find` over per-sample directories produces.

  **A track that lies by omission.** A missing interval means *no data* in bedGraph, so a
  zero-depth base left out reads as unmeasured, and a run merged across a coordinate gap claims
  coverage for bases never queried. Both look like a valid track.

The track's own correctness is pinned to an external number: `Σ (end−start)×value` over a region
must equal `samtools bedcov -j`'s sum for it. That is the identity which makes the removed
per-exon summary derivable rather than lost, so it is asserted here against bedcov directly.

Run it directly; `unittest`, no fixture on disk — every input is built in a temporary directory by
the test that needs it. The BAM is a two-exon read written as SAM and converted, so a CIGAR N is
present to test against; the tests needing it are skipped where `samtools` is absent.

    test_write_coverage_track.py            # quiet unless something fails
    test_write_coverage_track.py -v

Python 3.8+. Run from anywhere: `python3 -m unittest discover tests`, or `pytest tests/`.
Needs `ngsutils` importable (an editable install, or the repository root on PYTHONPATH — the
tests add it themselves) and `samtools` on PATH for all but the label tests.
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

# %%
# Constants

# Invoked through the dispatcher, not as a neighbouring script: `ngsutils write_coverage_track`
# is the only supported entry point, so that is what the CLI tests must exercise. `-m ngsutils.cli`
# rather than the installed `ngsutils` console script, so an editable install is not required and
# the tests read the working tree.
DISPATCH = [sys.executable, "-m", "ngsutils.cli", "write_coverage_track"]
REPO_ROOT = Path(__file__).resolve().parent.parent
SAMTOOLS = shutil.which("samtools")

CONTIG = "chr1"
CONTIG_LENGTH = 100

# One read, 10M20N10M at 1-based 11: exon 11-20, intron 21-40 skipped, exon 41-50. The N is the
# point — it is what separates "depth" from "bases the read spans".
READ_START = 11
EXON1 = (11, 20)
INTRON = (21, 40)
EXON2 = (41, 50)
EXONIC_BASES = (EXON1[1] - EXON1[0] + 1) + (EXON2[1] - EXON2[0] + 1)  # 20
SAM = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    f"@SQ\tSN:{CONTIG}\tLN:{CONTIG_LENGTH}\n"
    f"r1\t0\t{CONTIG}\t{READ_START}\t60\t10M20N10M\t*\t0\t0\t"
    "ACGTACGTACACGTACGTAC\tIIIIIIIIIIIIIIIIIIII\n"
)

# The whole read footprint as a 1-based inclusive region, intron included.
WHOLE_REGION = f"{CONTIG}:{EXON1[0]}-{EXON2[1]}"

# A region running past the read, so it contains positions no read overlaps at all.
#
# This distinction is the only thing that tests `-a`, and it was measured rather than assumed:
# inside the read's span, `samtools depth` emits the ref-skipped intron as depth 0 with or without
# `-a`, because the read overlaps those positions. Past the read's end nothing overlaps, and only
# `-a` reports them. A fixture whose region equals the read footprint therefore cannot tell the two
# apart — an earlier version of these tests could not, and passed with `-a` removed.
UNCOVERED = (EXON2[1] + 1, EXON2[1] + 10)
OVERRUN_REGION = f"{CONTIG}:{EXON1[0]}-{UNCOVERED[1]}"


# %%
# Helpers

def load_module():
    """The module itself, for the functions that need no subprocess."""
    sys.path.insert(0, str(REPO_ROOT))
    try:
        import ngsutils.write_coverage_track as module
    finally:
        sys.path.pop(0)
    return module


def run(*args) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO_ROOT), env["PYTHONPATH"]] if env.get("PYTHONPATH") else [str(REPO_ROOT)]
    )
    return subprocess.run(
        [*DISPATCH, *args], capture_output=True, text=True, env=env, cwd=str(REPO_ROOT)
    )


def refuse(test, module, *args):
    """`resolve_labels` expected to exit 2, with its message kept off the test's own stderr.

    `die()` writes to stderr by design and these tests provoke it on purpose; without the redirect
    the suite prints a screenful of correct error messages and stops being quiet-unless-failing.
    The wording is asserted through the CLI instead, in `LabelMessageTests`.
    """
    captured = io.StringIO()
    with contextlib.redirect_stderr(captured):
        with test.assertRaises(SystemExit) as caught:
            module.resolve_labels(*args)
    test.assertEqual(caught.exception.code, 2, captured.getvalue())
    return captured.getvalue()


def make_bam(directory: Path, name: str = "sample.bam") -> Path:
    """The two-exon read as an indexed BAM at `directory/name`."""
    directory.mkdir(parents=True, exist_ok=True)
    sam = directory / "in.sam"
    sam.write_text(SAM)
    bam = directory / name
    view = subprocess.run(
        [SAMTOOLS, "view", "-b", str(sam)], capture_output=True, check=True
    )
    subprocess.run(
        [SAMTOOLS, "sort", "-o", str(bam), "-"],
        input=view.stdout,
        capture_output=True,
        check=True,
    )
    subprocess.run([SAMTOOLS, "index", str(bam)], capture_output=True, check=True)
    return bam


def bedgraph(path: Path):
    """`(track_line_or_None, [(contig, start, end, value)])`."""
    track = None
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("track"):
            track = line
            continue
        if not line.strip():
            continue
        contig, start, end, value = line.split("\t")
        out.append((contig, int(start), int(end), int(value)))
    return track, out


def weighted_sum(rows, start: int, end: int) -> int:
    """`Σ (overlap length)×value` over `[start, end)` — the derivation the summary was replaced by.

    Length-weighted on purpose. `bedtools map -c 4 -o sum` adds the run-length-encoded values
    without weighting and returns 2 where this returns 20; the docstring of the script says so, and
    `test_the_track_total_equals_bedcov` is what holds it to the right answer.
    """
    total = 0
    for _, run_start, run_end, value in rows:
        lo, hi = max(run_start, start), min(run_end, end)
        if hi > lo:
            total += (hi - lo) * value
    return total


# %%
# Labels — no BAM, no samtools

class LabelTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_each_source_names_the_same_bam_differently(self):
        bams = ["s1/aligned.bam"]
        for source, expected in (
            ("stem", "aligned"),
            ("basename", "aligned.bam"),
            ("parent", "s1"),
            ("parent-stem", "s1_aligned"),
        ):
            with self.subTest(source=source):
                self.assertEqual(
                    self.module.resolve_labels(bams, [], source), [expected]
                )

    def test_a_relative_path_still_has_a_parent(self):
        """`./x.bam`'s parent is `.`, which names nothing — abspath first, or the label is unusable."""
        labels = self.module.resolve_labels(["./x.bam"], [], "parent")
        self.assertNotIn(labels[0], (".", ""), "the parent of a relative path must resolve")

    def test_the_find_layout_collapses_under_the_default_and_is_refused(self):
        """The whole reason labels exist: one file name per sample directory."""
        refuse(self, self.module, ["s1/aligned.bam", "s2/aligned.bam"], [], "stem")

    def test_parent_stem_separates_what_stem_collapses(self):
        self.assertEqual(
            self.module.resolve_labels(
                ["s1/aligned.bam", "s2/aligned.bam"], [], "parent-stem"
            ),
            ["s1_aligned", "s2_aligned"],
        )

    def test_explicit_labels_are_used_in_order(self):
        self.assertEqual(
            self.module.resolve_labels(
                ["s1/aligned.bam", "s2/aligned.bam"], ["tumour", "normal"], "stem"
            ),
            ["tumour", "normal"],
        )

    def test_a_label_that_would_be_a_path_is_refused(self):
        """The label is a filename, so `/` would write outside the directory or fail."""
        for bad in ("a/b", "a\\b"):
            with self.subTest(label=bad):
                refuse(self, self.module, ["x.bam"], [bad], "stem")

    def test_an_empty_label_is_refused(self):
        refuse(self, self.module, ["x.bam"], ["  "], "stem")


class LabelMessageTests(unittest.TestCase):
    """The message has to name the flag the caller used, or it sends them somewhere useless."""

    def test_a_stem_clash_points_at_parent_stem(self):
        r = run("--output-dir", "out", "--region", WHOLE_REGION, "s1/aligned.bam", "s2/aligned.bam")
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("duplicate label", r.stderr)
        self.assertIn("parent-stem", r.stderr)
        self.assertIn("s1/aligned.bam", r.stderr, "the message must name both paths")
        self.assertIn("s2/aligned.bam", r.stderr)
        self.assertEqual(r.stdout, "", "stdout carries no data at all, and least of all an error")

    def test_an_explicit_clash_does_not_point_at_label_from(self):
        r = run(
            "--output-dir", "out", "--region", WHOLE_REGION,
            "--label", "same", "--label", "same", "s1/a.bam", "s2/b.bam",
        )
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("--label", r.stderr)
        self.assertNotIn(
            "parent-stem", r.stderr, "--label-from cannot fix labels given by hand"
        )

    def test_a_parent_clash_asks_for_explicit_labels(self):
        r = run(
            "--output-dir", "out", "--region", WHOLE_REGION,
            "--label-from", "parent", "d/one.bam", "d/two.bam",
        )
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("--label", r.stderr)

    def test_a_label_count_mismatch_is_an_error(self):
        r = run(
            "--output-dir", "out", "--region", WHOLE_REGION,
            "--label", "only-one", "a.bam", "b.bam",
        )
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("--label", r.stderr)

    def test_label_and_label_from_together_are_an_error(self):
        r = run(
            "--output-dir", "out", "--region", WHOLE_REGION,
            "--label", "x", "--label-from", "parent", "a.bam",
        )
        self.assertEqual(r.returncode, 2, r.stderr)

    def test_an_unknown_label_source_lists_the_valid_ones(self):
        r = run(
            "--output-dir", "out", "--region", WHOLE_REGION, "--label-from", "bogus", "a.bam"
        )
        self.assertEqual(r.returncode, 2, r.stderr)
        for source in ("basename", "stem", "parent", "parent-stem"):
            self.assertIn(source, r.stderr)

    def test_the_docopt_defaults_survived_the_docstring(self):
        """docopt 0.6.2 reads any dash-initial docstring line as an option definition, silently
        discarding a `[default:]`. The script asserts its own defaults; this asserts the assertion
        has not been bypassed by a prose line added since."""
        r = run("--output-dir", "out", "--region", WHOLE_REGION, "missing.bam")
        self.assertNotIn("lost its default", r.stderr, r.stderr)


# %%
# The interface, which is positional on purpose

class InterfaceTests(unittest.TestCase):
    def test_bams_are_positional_and_repeatable(self):
        """`samtools depth [options] in.bam [in.bam ...]` takes them this way, and it is what makes
        `find | xargs` work with no argument rewriting. The BAMs are missing here, so the run fails
        at the file check — which is already past parsing."""
        r = run("--output-dir", "out", "--region", WHOLE_REGION, "a.bam", "b.bam", "c.bam")
        self.assertNotIn("Usage", r.stdout + r.stderr, "three BAMs must parse")
        self.assertIn("a.bam", r.stderr, "and fail on the first missing file, not on the syntax")

    def test_no_bam_is_a_usage_error(self):
        r = run("--output-dir", "out", "--region", WHOLE_REGION)
        self.assertEqual(r.returncode, 1)
        self.assertIn("Usage", r.stdout + r.stderr)

    def test_outdir_is_required(self):
        """The only output is files, so there is no run without somewhere to put them."""
        r = run("--region", WHOLE_REGION, "a.bam")
        self.assertEqual(r.returncode, 1)
        self.assertIn("Usage", r.stdout + r.stderr)

    def test_no_target_is_an_error(self):
        r = run("--output-dir", "out", "a.bam")
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("--id or --region", r.stderr)

    def test_an_id_without_a_gtf_is_an_error(self):
        r = run("--output-dir", "out", "--id", "ENST00000269305", "a.bam")
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("--gtf", r.stderr)


# %%
# The track, on a real BAM

@unittest.skipUnless(SAMTOOLS, "samtools not on PATH")
class TrackTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.bam = make_bam(self.dir)
        self.out = self.dir / "tracks"

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, *args, region: str = WHOLE_REGION):
        r = run("--output-dir", str(self.out), "--region", region, *args, str(self.bam))
        self.assertEqual(r.returncode, 0, r.stderr)
        return r

    def track(self, label: str = "sample"):
        return bedgraph(self.out / f"{label}.bedgraph")

    def test_the_track_reproduces_the_read_structure(self):
        self.write()
        line, rows = self.track()
        self.assertIsNotNone(line, "a track line by default")
        self.assertEqual(
            rows,
            [
                (CONTIG, EXON1[0] - 1, EXON1[1], 1),
                (CONTIG, INTRON[0] - 1, INTRON[1], 0),
                (CONTIG, EXON2[0] - 1, EXON2[1], 1),
            ],
            "three intervals: exon, intron at zero, exon",
        )

    def test_a_one_based_region_becomes_zero_based_half_open(self):
        """The off-by-one this invites: `--region` is 1-based inclusive, bedGraph is BED-like."""
        self.write()
        _, rows = self.track()
        self.assertEqual(rows[0][1], EXON1[0] - 1)
        self.assertEqual(rows[-1][2], EXON2[1])

    def test_an_intron_skip_is_not_counted_as_depth(self):
        """A junction-spanning read must not inflate an exon it never touched."""
        self.write()
        _, rows = self.track()
        inside = [r for r in rows if r[1] == INTRON[0] - 1]
        self.assertEqual(inside[0][3], 0, "the ref-skipped intron is depth 0")
        self.assertEqual(
            weighted_sum(rows, EXON1[0] - 1, EXON2[1]),
            EXONIC_BASES,
            "20 exonic bases at depth 1, and nothing from the intron",
        )

    def test_zero_depth_is_written_rather_than_omitted(self):
        """An absent interval means *no data*, so a zero left out is unmeasured. The intron alone
        does not test `-a` — see `UNCOVERED` — but it is the zero a reader will look for."""
        self.write()
        _, rows = self.track()
        zeros = [r for r in rows if r[3] == 0]
        self.assertEqual(len(zeros), 1)
        self.assertEqual((zeros[0][1], zeros[0][2]), (INTRON[0] - 1, INTRON[1]))

    def test_a_position_no_read_overlaps_is_still_reported(self):
        """The test that `-a` is passed at all. Past the read's end nothing overlaps, so without
        `-a` those positions vanish and the track stops early with nothing saying they were
        measured and empty."""
        self.write(region=OVERRUN_REGION)
        _, rows = self.track()
        self.assertEqual(
            rows[-1],
            (CONTIG, EXON2[1], UNCOVERED[1], 0),
            f"the uncovered tail must be one 0-valued interval, got {rows}",
        )

    def test_the_track_has_no_hole(self):
        self.write()
        _, rows = self.track()
        for earlier, later in zip(rows, rows[1:]):
            self.assertEqual(earlier[2], later[1], "abutting, so no base is unreported")

    def test_a_run_breaks_at_a_coordinate_gap(self):
        """Two disjoint regions of equal depth must stay two intervals: merging them would claim
        coverage for the bases between, which were never queried."""
        r = run(
            "--output-dir", str(self.out),
            "--region", f"a={CONTIG}:{EXON1[0]}-{EXON1[0] + 4}",
            "--region", f"b={CONTIG}:{EXON2[0] + 4}-{EXON2[1]}",
            str(self.bam),
        )
        self.assertEqual(r.returncode, 0, r.stderr)
        _, rows = self.track()
        self.assertEqual(len(rows), 2, f"one interval per region, got {rows}")
        self.assertEqual({value for *_, value in rows}, {1}, "same depth, still not merged")

    def test_overlapping_regions_are_unioned(self):
        """Overlapping input would make samtools depth report a position twice, and two bedGraph
        intervals over one base is malformed."""
        r = run(
            "--output-dir", str(self.out),
            "--region", f"{CONTIG}:{EXON1[0]}-{INTRON[1]}",
            "--region", f"{CONTIG}:{INTRON[0]}-{EXON2[1]}",
            str(self.bam),
        )
        self.assertEqual(r.returncode, 0, r.stderr)
        _, rows = self.track()
        for earlier, later in zip(rows, rows[1:]):
            self.assertLessEqual(earlier[2], later[1], f"intervals overlap: {rows}")

    def test_the_track_total_equals_bedcov(self):
        """The identity that makes the removed per-exon summary derivable rather than lost.

        `samtools bedcov -j` is asked directly, so this is a comparison against another
        implementation of the same quantity and not against a number written down here.
        """
        self.write()
        _, rows = self.track()
        bed = self.dir / "region.bed"
        bed.write_text(f"{CONTIG}\t{EXON1[0] - 1}\t{EXON2[1]}\twhole\t0\t+\n")
        p = subprocess.run(
            [SAMTOOLS, "bedcov", "-j", str(bed), str(self.bam)],
            capture_output=True, text=True, check=True,
        )
        bedcov_sum = int(p.stdout.split("\t")[-1])
        self.assertEqual(weighted_sum(rows, EXON1[0] - 1, EXON2[1]), bedcov_sum)

    def test_the_file_and_the_track_name_are_the_label(self):
        self.write("--label", "tumour")
        self.assertTrue((self.out / "tumour.bedgraph").exists())
        self.assertFalse((self.out / "sample.bedgraph").exists())
        line, _ = self.track("tumour")
        self.assertIn('name="tumour"', line)

    def test_no_track_line_omits_it(self):
        self.write("--no-track-line")
        line, rows = self.track()
        self.assertIsNone(line)
        self.assertTrue(rows, "the data is still there")

    def test_two_same_named_bams_get_one_file_each_under_parent_stem(self):
        """`parent-stem` on the layout `find` produces. This does **not** exercise the duplicate
        guard — parent-stem never collides here — see `test_the_default_refuses_to_collapse`."""
        first = make_bam(self.dir / "s1", "aligned.bam")
        second = make_bam(self.dir / "s2", "aligned.bam")
        r = run(
            "--output-dir", str(self.out), "--region", WHOLE_REGION,
            "--label-from", "parent-stem", str(first), str(second),
        )
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(
            sorted(p.name for p in self.out.glob("*.bedgraph")),
            ["s1_aligned.bedgraph", "s2_aligned.bedgraph"],
        )

    def test_the_default_refuses_to_collapse(self):
        """The collapse the labels exist to prevent, with real files: one bedGraph would be written
        where two were asked for, and every number in it would be correct."""
        first = make_bam(self.dir / "s1", "aligned.bam")
        second = make_bam(self.dir / "s2", "aligned.bam")
        r = run("--output-dir", str(self.out), "--region", WHOLE_REGION, str(first), str(second))
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("duplicate label", r.stderr)
        self.assertFalse(
            self.out.exists(), "no directory is created before the labels are known good"
        )

    def test_each_bam_gets_its_own_depth_column(self):
        """One samtools pass, one column per BAM, and column i must reach file i."""
        other = self.dir / "empty"
        other.mkdir()
        # A BAM with no reads over the region: its track must be all zeros while the first is not.
        sam = other / "in.sam"
        sam.write_text(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CONTIG}\tLN:{CONTIG_LENGTH}\n")
        view = subprocess.run([SAMTOOLS, "view", "-b", str(sam)], capture_output=True, check=True)
        empty_bam = other / "empty.bam"
        subprocess.run(
            [SAMTOOLS, "sort", "-o", str(empty_bam), "-"],
            input=view.stdout, capture_output=True, check=True,
        )
        subprocess.run([SAMTOOLS, "index", str(empty_bam)], capture_output=True, check=True)

        r = run(
            "--output-dir", str(self.out), "--region", WHOLE_REGION,
            str(self.bam), str(empty_bam),
        )
        self.assertEqual(r.returncode, 0, r.stderr)
        _, first = self.track("sample")
        _, second = self.track("empty")
        self.assertEqual(
            weighted_sum(first, EXON1[0] - 1, EXON2[1]), EXONIC_BASES, "the read's own track"
        )
        self.assertEqual(
            weighted_sum(second, EXON1[0] - 1, EXON2[1]), 0, "the empty BAM's track is all zeros"
        )
        self.assertEqual(len(second), 1, "one 0-valued interval across the whole region")

    def test_threads_reaches_samtools_without_changing_the_answer(self):
        """A passthrough flag that is silently dropped looks exactly like one that is honoured, so
        the command line is asserted — and the depth values must not move with the thread count."""
        plain = self.write()
        _, without = self.track()
        (self.out / "sample.bedgraph").unlink()
        threaded = self.write("--threads", "2")
        _, with_threads = self.track()
        self.assertIn("-@ 2", threaded.stderr, "the flag must reach samtools depth")
        self.assertNotIn("-@", plain.stderr, "and must be absent when not asked for")
        self.assertEqual(without, with_threads, "thread count is not part of the answer")

    def test_a_non_integer_thread_count_is_refused(self):
        r = run("--output-dir", str(self.out), "--region", WHOLE_REGION, "--threads", "abc",
                str(self.bam))
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("--threads", r.stderr)

    def test_an_unindexed_bam_is_an_error_naming_the_fix(self):
        unindexed = self.dir / "noindex.bam"
        shutil.copy(self.bam, unindexed)
        r = run("--output-dir", str(self.out), "--region", WHOLE_REGION, str(unindexed))
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("samtools index", r.stderr)


# %%
# Seeking — the `-r` path, and the two ways it silently goes wrong

# A second contig whose NAME sorts before the first but whose HEADER position is after it. That
# disagreement is the whole point: `samtools depth -b` emits rows in header order, so a `-r` loop
# driven by the union's coordinate sort would emit these two the other way round and quietly
# rewrite the bedGraph.
SEEK_SAM = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr2\tLN:2000\n"
    "@SQ\tSN:chr10\tLN:2000\n"
    "a\t0\tchr2\t301\t60\t50M\t*\t0\t0\t" + "A" * 50 + "\t" + "I" * 50 + "\n"
    "b\t0\tchr10\t101\t60\t50M\t*\t0\t0\t" + "A" * 50 + "\t" + "I" * 50 + "\n"
)


@unittest.skipUnless(SAMTOOLS, "samtools not on PATH")
class SeekTests(unittest.TestCase):
    """`-r` must buy speed without changing a byte of the answer.

    `samtools depth -b BED` does not use the index — verified by the fact that it runs on a BAM
    with none, while `-r` refuses to — so it reads the whole file. Adding `-r` makes it seek. The
    two failure modes below are what makes that substitution dangerous, and neither announces
    itself: both would exit 0 and write a track that looks complete.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.out = self.dir / "tracks"
        sam = self.dir / "seek.sam"
        sam.write_text(SEEK_SAM)
        self.bam = self.dir / "seek.bam"
        view = subprocess.run(
            [SAMTOOLS, "view", "-b", str(sam)], capture_output=True, check=True
        )
        subprocess.run(
            [SAMTOOLS, "sort", "-o", str(self.bam), "-"],
            input=view.stdout,
            capture_output=True,
            check=True,
        )
        subprocess.run([SAMTOOLS, "index", str(self.bam)], capture_output=True, check=True)

    def tearDown(self):
        self.tmp.cleanup()

    def depth_reference(self) -> str:
        """`samtools depth -a -b BED` with no `-r` at all: the answer the program must reproduce."""
        bed = self.dir / "ref.bed"
        bed.write_text("chr10\t100\t130\nchr2\t300\t320\n")
        p = subprocess.run(
            [SAMTOOLS, "depth", "-a", "-b", str(bed), str(self.bam)],
            capture_output=True,
            text=True,
            check=True,
        )
        return p.stdout

    def test_both_contigs_survive(self):
        """A single `-r` span across contigs keeps the first and drops the rest, exit 0.

        Measured on this fixture: 50 rows become 30. The regions must therefore be one per contig,
        and the check is that every queried base comes back — the program's own row-count assertion
        would also fire, but this names the reason.
        """
        r = run("--output-dir", str(self.out), "--region", "chr10:101-130",
                "--region", "chr2:301-320", str(self.bam))
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = self.depth_reference().splitlines()
        self.assertEqual(len(rows), 50, "fixture changed; the reference is no longer 50 positions")
        contigs = {line.split("\t")[0] for line in rows}
        self.assertEqual(contigs, {"chr2", "chr10"}, "the reference must span both contigs")
        written = sorted(p.name for p in self.out.glob("*.bedgraph"))
        self.assertEqual(written, ["seek.bedgraph"])
        body = [
            line for line in (self.out / "seek.bedgraph").read_text().splitlines()
            if not line.startswith("track")
        ]
        covered = sum(int(f[2]) - int(f[1]) for f in (line.split("\t") for line in body))
        self.assertEqual(covered, 50, f"positions were dropped:\n{body}")

    def test_rows_come_back_in_header_order_not_name_order(self):
        """`chr2` before `chr10`, because that is what `-b` alone does and the file must not move.

        Name order would put `chr10` first — '1' < '2' — so this fixture fails loudly if the `-r`
        loop is ever driven by the union's sort instead of the BAM header.
        """
        r = run("--output-dir", str(self.out), "--region", "chr10:101-130",
                "--region", "chr2:301-320", str(self.bam))
        self.assertEqual(r.returncode, 0, r.stderr)
        body = [
            line for line in (self.out / "seek.bedgraph").read_text().splitlines()
            if not line.startswith("track")
        ]
        seen = []
        for line in body:
            contig = line.split("\t")[0]
            if not seen or seen[-1] != contig:
                seen.append(contig)
        self.assertEqual(seen, ["chr2", "chr10"], "bedGraph interval order changed")

    def test_every_call_carries_r(self):
        """The point of the change. Without `-r` samtools reads the whole BAM."""
        r = run("--output-dir", str(self.out), "--region", "chr10:101-130",
                "--region", "chr2:301-320", str(self.bam))
        self.assertEqual(r.returncode, 0, r.stderr)
        calls = [line for line in r.stderr.splitlines() if "running:" in line]
        self.assertEqual(len(calls), 2, f"one depth call per contig expected:\n{r.stderr}")
        for call in calls:
            self.assertIn(" -r ", call, "a call without -r scans the whole file")
            self.assertIn(" -b ", call, "-r selects what is read; -b still selects what is reported")

    def test_a_contig_absent_from_the_bam_is_refused(self):
        """Not silently empty. `chr1` here is a naming mismatch, which is the common real cause."""
        r = run("--output-dir", str(self.out), "--region", "chr1:1-10", str(self.bam))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not in the BAM header", r.stderr)


# %%
# Modes — which positions a target resolves to, and what the default is

# Two isoforms of one gene sharing an intron. Isoform A: 11-20, 41-50. Isoform B: 11-20, 61-70.
# So the exon union is 11-20, 41-50, 61-70 (30 bp) and the locus span is 11-70 (60 bp) — the two
# differ by exactly the introns, which is what these tests are about.
GENE_GTF = "\n".join(
    f'{CONTIG}\tsrc\texon\t{start}\t{end}\t.\t+\t.\t'
    f'gene_id "ENSG00000000001"; transcript_id "{tx}"; gene_name "TESTG";'
    for tx, start, end in [
        ("ENST00000000001", 11, 20),
        ("ENST00000000001", 41, 50),
        ("ENST00000000002", 11, 20),
        ("ENST00000000002", 61, 70),
    ]
) + "\n"

EXON_UNION_BP = 30   # 11-20, 41-50, 61-70
LOCUS_SPAN_BP = 60   # 11-70, introns included


@unittest.skipUnless(SAMTOOLS, "samtools not on PATH")
class ModeTests(unittest.TestCase):
    """`--mode` decides whether introns are in the track, and nothing tested it until 2026-08-12.

    That gap is why changing the default from `auto` to `span` broke no test. The default is the
    one behaviour every caller gets without asking for it, so it is pinned first.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.out = self.dir / "tracks"
        self.bam = make_bam(self.dir)
        self.gtf = self.dir / "test.gtf"
        self.gtf.write_text(GENE_GTF)

    def tearDown(self):
        self.tmp.cleanup()

    def covered(self, *extra) -> int:
        """Positions in the written track, which is what a mode ultimately decides."""
        r = run("--output-dir", str(self.out), "--gtf", str(self.gtf), "--id",
                "ENSG00000000001", *extra, str(self.bam))
        self.assertEqual(r.returncode, 0, r.stderr)
        body = [
            line for line in (self.out / "sample.bedgraph").read_text().splitlines()
            if not line.startswith("track")
        ]
        return sum(int(f[2]) - int(f[1]) for f in (line.split("\t") for line in body))

    def test_the_default_is_the_locus_with_introns(self):
        """A default track must be able to show a retained intron, so it spans the locus.

        `samtools depth` does not count CIGAR N, so a spliced read contributes 0 inside an intron
        and only an intron-retaining read leaves depth there. An exonic default cannot show that
        at all — the positions are simply not queried.
        """
        self.assertEqual(self.covered(), LOCUS_SPAN_BP)

    def test_merged_exon_cuts_the_introns_out(self):
        self.assertEqual(self.covered("--mode", "merged_exon"), EXON_UNION_BP)

    def test_auto_is_still_exonic_for_a_gene(self):
        """`auto` did not change meaning; it stopped being the default."""
        self.assertEqual(self.covered("--mode", "auto"), EXON_UNION_BP)

    def test_span_and_the_default_agree(self):
        self.assertEqual(self.covered("--mode", "span"), LOCUS_SPAN_BP)

    def test_the_old_mode_name_is_refused(self):
        """`merged` was renamed to `merged_exon`. A stale script must fail, not fall back.

        Silently reinterpreting it would be the worse outcome of the two: the caller would get a
        track and no indication that the name they used no longer exists.
        """
        r = run("--output-dir", str(self.out), "--gtf", str(self.gtf), "--id",
                "ENSG00000000001", "--mode", "merged", str(self.bam))
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("merged_exon", r.stderr, "the error must name the replacement")


if __name__ == "__main__":
    unittest.main()
