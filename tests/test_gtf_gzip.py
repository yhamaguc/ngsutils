#!/usr/bin/env python3
"""test_gtf_gzip — a gzipped GTF must behave exactly as the plain file it decompresses to.

Two silent failures, one per class below.

  **A truncated read.** polars decompresses gzip itself, so `read_gtf("x.gtf.gz")` needs no
  code of ours -- but a single-member gzip decoder reads only the first member of a
  concatenated or bgzip-compressed file and returns a SHORT TABLE, not an error. Every
  downstream count, interval and mapping would then be computed from part of the
  annotation. The test builds a two-member file and asserts the row count.

  **An output named after the wrong stem.** `os.path.splitext("x.gtf.gz")` strips only
  ".gz", so the subcommands that name their output after their input wrote `x.gtf.bed` and
  `x.gtf.sqlite` for a gzipped input and `x.bed` / `x.sqlite` for the plain one. Nothing
  failed; the caller simply got a file under a name it did not look for.
  `ngsutils.gtf.gtf_stem` is the one place that decides, and these tests run the commands
  both ways and assert the same filename and the same bytes.

    python3 -m unittest discover tests
    pytest tests/test_gtf_gzip.py

Needs `ngsutils` importable. Every input is built in a temporary directory; nothing in the
repository is read or written.
"""

from __future__ import annotations

import gzip
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ngsutils.gtf import gtf_stem, read_gtf  # noqa: E402

# NOTE: one transcript, two exons, with the attributes gtf2bed and gene2tx require. Written
#   out rather than read from tests/data so that a change to a shared fixture cannot move
#   the row count this file asserts.
GTF_LINES = [
    "##description: test_gtf_gzip fixture",
    'chr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG00000000001.1"; gene_type "protein_coding"; gene_name "AAA";',
    'chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding"; transcript_name "AAA-201";',
    'chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1"; exon_number 1; exon_id "ENSE00000000001.1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding"; transcript_name "AAA-201";',
    'chr1\tHAVANA\tCDS\t100\t200\t.\t+\t0\tgene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1"; exon_number 1; exon_id "ENSE00000000001.1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding"; transcript_name "AAA-201";',
    'chr1\tHAVANA\tstart_codon\t100\t102\t.\t+\t0\tgene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1"; exon_number 1; exon_id "ENSE00000000001.1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding"; transcript_name "AAA-201";',
    'chr1\tHAVANA\texon\t400\t500\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1"; exon_number 2; exon_id "ENSE00000000002.1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding"; transcript_name "AAA-201";',
    'chr1\tHAVANA\tstop_codon\t498\t500\t.\t+\t0\tgene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1"; exon_number 2; exon_id "ENSE00000000002.1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding"; transcript_name "AAA-201";',
]
GTF_TEXT = "\n".join(GTF_LINES) + "\n"

# The comment line is not a record, so read_gtf must return one row per remaining line.
EXPECTED_ROWS = len(GTF_LINES) - 1

# NOTE: a second gene with gene/transcript/exon rows and NO CDS and no codon rows, which
#   is what an lncRNA, an Mt_tRNA or a retained_intron transcript looks like in GENCODE.
#   Kept out of GTF_LINES so that the single-gene expectations above stay single-gene;
#   write_mixed() below is the fixture that carries both.
NONCODING_ATTRIBUTES = (
    'gene_id "ENSG00000000002.1"; transcript_id "ENST00000000002.1"; '
    'gene_type "lncRNA"; gene_name "BBB"; transcript_type "lncRNA"; '
    'transcript_name "BBB-201";'
)
NONCODING_LINES = [
    'chr1\tHAVANA\tgene\t1000\t1500\t.\t-\t.\tgene_id "ENSG00000000002.1"; gene_type "lncRNA"; gene_name "BBB";',
    f'chr1\tHAVANA\ttranscript\t1000\t1500\t.\t-\t.\t{NONCODING_ATTRIBUTES}',
    f'chr1\tHAVANA\texon\t1000\t1100\t.\t-\t.\t{NONCODING_ATTRIBUTES} exon_number 1; exon_id "ENSE00000000003.1";',
    f'chr1\tHAVANA\texon\t1400\t1500\t.\t-\t.\t{NONCODING_ATTRIBUTES} exon_number 2; exon_id "ENSE00000000004.1";',
]

CODING_TRANSCRIPT = "ENST00000000001.1"
NONCODING_TRANSCRIPT = "ENST00000000002.1"

# NOTE: a GTF with no exon_id attribute ANYWHERE, which is what StringTie and several
#   other assemblers emit. gtf2bed has a branch that supplies synthetic ids for that case;
#   one attribute left on one row would create the column and skip the branch, so every
#   line here has to be free of it.
NO_EXON_ID_ATTRIBUTES = (
    'gene_id "GENE1"; transcript_id "TX1"; gene_type "protein_coding"; '
    'gene_name "CCC"; transcript_type "protein_coding"; transcript_name "CCC-201";'
)
NO_EXON_ID_LINES = [
    'chr2\tStringTie\tgene\t100\t500\t.\t+\t.\tgene_id "GENE1"; gene_type "protein_coding"; gene_name "CCC";',
    f'chr2\tStringTie\ttranscript\t100\t500\t.\t+\t.\t{NO_EXON_ID_ATTRIBUTES}',
    f'chr2\tStringTie\texon\t100\t200\t.\t+\t.\t{NO_EXON_ID_ATTRIBUTES} exon_number 1;',
    f'chr2\tStringTie\texon\t400\t500\t.\t+\t.\t{NO_EXON_ID_ATTRIBUTES} exon_number 2;',
]

FIXTURE_STEM = "probe"


def write_plain(directory, name=f"{FIXTURE_STEM}.gtf"):
    path = os.path.join(directory, name)
    with io.open(path, "w") as f:
        f.write(GTF_TEXT)
    return path


def write_gzipped(directory, name=f"{FIXTURE_STEM}.gtf.gz"):
    path = os.path.join(directory, name)
    with gzip.open(path, "wt") as f:
        f.write(GTF_TEXT)
    return path


def write_mixed(directory, name=f"{FIXTURE_STEM}.gtf"):
    """The coding gene plus the non-coding one, which is what a real GENCODE file holds."""
    path = os.path.join(directory, name)
    with io.open(path, "w") as f:
        f.write(GTF_TEXT + "\n".join(NONCODING_LINES) + "\n")
    return path


def write_without_exon_id(directory, name=f"{FIXTURE_STEM}.gtf"):
    """A GTF whose exon rows carry no exon_id attribute at all."""
    path = os.path.join(directory, name)
    with io.open(path, "w") as f:
        f.write("\n".join(NO_EXON_ID_LINES) + "\n")
    return path


def write_multi_member_gzip(directory, name=f"{FIXTURE_STEM}.gtf.gz"):
    """One gzip file of two members, the shape bgzip produces.

    Built with the gzip module rather than by calling bgzip so that the test does not
    depend on htslib being installed. The decoder path is the same: a reader that stops at
    the first member's end-of-stream sees only part of the file.
    """
    path = os.path.join(directory, name)
    head, tail = GTF_TEXT.split("\n", 1)
    with io.open(path, "wb") as out:
        out.write(gzip.compress((head + "\n").encode()))
        out.write(gzip.compress(tail.encode()))
    return path


def run_subcommand(module, argv, work_dir):
    """Run one subcommand out of the repository, not out of whatever is installed."""
    environment = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
    return subprocess.run(
        [sys.executable, "-m", f"ngsutils.{module}"] + argv,
        capture_output=True,
        text=True,
        cwd=work_dir,
        env=environment,
    )


# %%
# The stem, which is what names an output after its input

class GtfStemTests(unittest.TestCase):
    def test_the_stem_is_the_same_whether_the_input_is_compressed(self):
        """The regression this file exists for: ".gtf" survived splitext on "x.gtf.gz"."""
        self.assertEqual(gtf_stem("x.gtf"), gtf_stem("x.gtf.gz"))

    def test_the_declared_suffixes_are_stripped(self):
        cases = {
            "probe.gtf": "probe",
            "probe.gtf.gz": "probe",
            "probe.gtf.bgz": "probe",
            "probe.GTF.GZ": "probe",
            "/a/b/probe.gtf.gz": "probe",
            "gencode.v50.annotation.gtf.gz": "gencode.v50.annotation",
        }
        for path, expected in sorted(cases.items()):
            with self.subTest(path=path):
                self.assertEqual(gtf_stem(path), expected)

    def test_a_dot_that_is_not_a_suffix_is_kept(self):
        """"gencode.v50" must not become "gencode" -- the release is part of the name."""
        self.assertEqual(gtf_stem("gencode.v50.gtf.gz"), "gencode.v50")


# %%
# read_gtf, which polars decompresses for

class ReadGtfCompressionTests(unittest.TestCase):
    def test_a_gzipped_gtf_reads_as_the_plain_file_does(self):
        with tempfile.TemporaryDirectory() as directory:
            plain = read_gtf(write_plain(directory))
            gzipped = read_gtf(write_gzipped(directory))
        self.assertEqual(plain.shape, gzipped.shape)
        self.assertTrue(plain.equals(gzipped), "gzipped and plain GTF gave different tables")

    def test_the_row_count_is_every_record(self):
        """Guards the truncation this file's docstring names: a short table is not an error."""
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(read_gtf(write_plain(directory)).height, EXPECTED_ROWS)

    def test_a_multi_member_gzip_is_read_to_its_end(self):
        """bgzip output is concatenated gzip members. Reading one member loses the rest."""
        with tempfile.TemporaryDirectory() as directory:
            frame = read_gtf(write_multi_member_gzip(directory))
        self.assertEqual(
            frame.height,
            EXPECTED_ROWS,
            "a multi-member gzip was read short; every downstream count would be wrong",
        )

    def test_a_missing_file_is_an_error_and_not_an_empty_table(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                read_gtf(os.path.join(directory, "absent.gtf.gz"))


# %%
# The subcommands, run both ways

class SubcommandCompressionTests(unittest.TestCase):
    """Each command below names its output after its input, which is where the stem matters.

    extract_splice_sites and gtf2tsv write to stdout and so name nothing; they are covered
    by comparing their output between the two input forms.
    """

    def both_ways(self, module, argv_for, produced_name):
        """Run `module` on the plain and the gzipped fixture; return the two output paths."""
        outputs = {}
        for label, writer in (("plain", write_plain), ("gzip", write_gzipped)):
            directory = tempfile.mkdtemp()
            self.addCleanup(lambda d=directory: __import__("shutil").rmtree(d, ignore_errors=True))
            gtf_path = writer(directory)
            result = run_subcommand(module, argv_for(gtf_path), directory)
            self.assertEqual(
                result.returncode,
                0,
                f"{module} failed on {label} input:\n{result.stderr}",
            )
            expected = os.path.join(directory, produced_name)
            made = sorted(p for p in os.listdir(directory) if p != os.path.basename(gtf_path))
            self.assertTrue(
                os.path.exists(expected),
                f"{module} on {label} input wrote {made} and not {produced_name!r}",
            )
            outputs[label] = expected
        return outputs

    def assert_same_bytes(self, outputs):
        with io.open(outputs["plain"], "rb") as f:
            plain = f.read()
        with io.open(outputs["gzip"], "rb") as f:
            gzipped = f.read()
        self.assertEqual(plain, gzipped, "the two input forms produced different output")

    def test_gtf2bed_names_and_writes_the_same_bed(self):
        outputs = self.both_ways("gtf2bed", lambda p: [p], f"{FIXTURE_STEM}.bed")
        self.assert_same_bytes(outputs)

    def test_gene2tx_names_and_writes_the_same_table(self):
        outputs = self.both_ways(
            "gene2tx",
            lambda p: ["--output-dir", os.path.dirname(p), p],
            f"{FIXTURE_STEM}.gene2tx.txt",
        )
        self.assert_same_bytes(outputs)

    def test_abs2rel_names_the_same_bed(self):
        # NOTE: abs2rel did not run at all until 2026-09-09 -- pandas rejects
        #   `on=` together with `left_index=True` -- so this is the first coverage the
        #   command has. What it writes is checked by Abs2relOutputTests below.
        outputs = self.both_ways("abs2rel", lambda p: [p], f"{FIXTURE_STEM}_2.bed")
        self.assert_same_bytes(outputs)

    def test_gtf2sqlite_names_the_same_database(self):
        self.both_ways("gtf2sqlite", lambda p: [p], f"{FIXTURE_STEM}.sqlite")

    def test_the_stdout_commands_agree_between_the_two_input_forms(self):
        for module in ("extract_splice_sites", "gtf2tsv"):
            with self.subTest(module=module):
                written = {}
                for label, writer in (("plain", write_plain), ("gzip", write_gzipped)):
                    with tempfile.TemporaryDirectory() as directory:
                        gtf_path = writer(directory)
                        result = run_subcommand(module, [gtf_path], directory)
                        self.assertEqual(
                            result.returncode,
                            0,
                            f"{module} failed on {label} input:\n{result.stderr}",
                        )
                        written[label] = result.stdout
                self.assertEqual(written["plain"], written["gzip"])
                self.assertTrue(written["plain"].strip(), f"{module} wrote nothing")


# %%
# What abs2rel writes, now that it runs

class Abs2relOutputTests(unittest.TestCase):
    """abs2rel is BED12 in gene-relative coordinates, and the fixture pins what that means.

    The gene spans 100..500 in the 1-based inclusive GTF frame, which is 401 bases, and
    carries two 101-base exons at 100..200 and 400..500. In the 0-based half-open BED frame
    that record is [0, 401), its blocks are 101,101 at offsets 0,300.
    """

    GENE_LENGTH = 401
    BLOCK_SIZES = "101,101,"
    BLOCK_STARTS = "0,300,"

    def bed_record(self):
        directory = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(directory, ignore_errors=True))
        gtf_path = write_plain(directory)
        result = run_subcommand("abs2rel", [gtf_path], directory)
        self.assertEqual(result.returncode, 0, result.stderr)
        with io.open(os.path.join(directory, f"{FIXTURE_STEM}_2.bed")) as f:
            rows = [line.rstrip("\n") for line in f if not line.startswith("#")]
        self.assertEqual(len(rows), 1, f"expected one gene record, got {rows}")
        return rows[0].split("\t")

    def test_the_blocks_are_the_two_exons(self):
        """The part that was already right, and the reference the record end is checked against."""
        record = self.bed_record()
        self.assertEqual(record[9], "2")
        self.assertEqual(record[10], self.BLOCK_SIZES)
        self.assertEqual(record[11], self.BLOCK_STARTS)

    def test_the_record_ends_at_the_gene_length(self):
        """abs2rel wrote `ends_ - starts_ + 1` until 2026-09-09, which is end - start + 2.

        `starts_` is already `start - 1`, so adding one again counted one base twice and
        the record ended at 402 for a 401-base gene. Nothing raised; the file was written
        and a reader took the span at face value.
        """
        record = self.bed_record()
        self.assertEqual(int(record[2]), self.GENE_LENGTH)

    def test_the_last_block_ends_at_the_record_end(self):
        """The BED12 invariant a viewer relies on, and the one the off-by-one broke.

        blockStart + blockSize of the final block must equal end - start. It read 401
        against a record end of 402, so the record was not valid BED12.
        """
        record = self.bed_record()
        start, end = int(record[1]), int(record[2])
        last_start = int(self.BLOCK_STARTS.rstrip(",").split(",")[-1])
        last_size = int(self.BLOCK_SIZES.rstrip(",").split(",")[-1])
        self.assertEqual(last_start + last_size, end - start)

    def test_the_thick_range_lies_inside_the_record(self):
        """abs2rel set both thick columns to `ends_`, an absolute GTF end, until 2026-09-09.

        For this fixture that was 500 against a record ending at 402, so the coding range
        sat outside the feature it belonged to. abs2rel reads no CDS, so BED's "nothing
        marked" -- thick_start == thick_end == start -- is what it can honestly write.
        """
        record = self.bed_record()
        start, end = int(record[1]), int(record[2])
        thick_start, thick_end = int(record[6]), int(record[7])
        self.assertLessEqual(start, thick_start)
        self.assertLessEqual(thick_start, thick_end)
        self.assertLessEqual(thick_end, end)

    def test_no_coding_range_is_marked(self):
        """The only honest thick range for a command that never reads a CDS feature."""
        record = self.bed_record()
        self.assertEqual(int(record[6]), int(record[7]))

    def test_the_gene_id_is_the_reference_name_unchanged(self):
        """`chr.str.replace('^ENSG', 'IMMT')` sat at the end of abs2rel until 2026-09-09.

        pandas 2.0 made `regex` default to False, so '^ENSG' was matched literally and
        replaced nothing -- the line was dead. Had it fired it would have rewritten every
        gene_id in the output. The name written is the gene_id from the GTF.
        """
        record = self.bed_record()
        self.assertEqual(record[0], "ENSG00000000001.1")


# %%
# gtf2bed's thick columns, which are where a non-coding transcript used to kill the run

class Gtf2bedThickTests(unittest.TestCase):
    """The two failures a GTF holding one non-coding transcript exposes.

      **A KeyError on every real GENCODE file.** `cds_start_min_` is built from CDS rows
      alone, and `cds_start_min_[idx]` asked it for EVERY transcript. A label lookup with
      a missing label raises `KeyError: "[...] not in index"`, so one lncRNA, Mt_tRNA or
      retained_intron transcript aborted the command. The one-gene fixture above could not
      see it: its only transcript has a CDS.

      **A coding range that was never written.** The override was
      `bed_df["thick_start"].where(..., inplace=True)`, which under pandas copy-on-write
      writes to a temporary and leaves bed_df untouched. pandas 3.0 reports it as
      ChainedAssignmentError -- a warning -- so the run looked clean and every transcript
      came out with no coding range marked, coding ones included.

    Fixture arithmetic. The coding transcript spans 100..500, its CDS 100..200, its
    start_codon 100..102 and its stop_codon 498..500. gtf2bed takes a base off each start
    and then the min and max over both columns, so the CDS bounds are 99..200 and the
    codon bounds 99..500, and the codons override the CDS. The non-coding transcript spans
    1000..1500 and has neither, so it keeps thick_start == thick_end == start == 999.
    """

    CODING_THICK = (99, 500)
    NONCODING_START = 999

    def records(self):
        """gtf2bed output for the mixed fixture, as {name field -> split record}."""
        directory = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(directory, ignore_errors=True))
        gtf_path = write_mixed(directory)
        result = run_subcommand("gtf2bed", [gtf_path], directory)
        self.assertEqual(
            result.returncode,
            0,
            f"gtf2bed failed on a GTF holding a non-coding transcript:\n{result.stderr}",
        )
        with io.open(os.path.join(directory, f"{FIXTURE_STEM}.bed")) as f:
            rows = [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]
        return {row[3]: row for row in rows}

    def transcript_record(self, transcript_id):
        records = self.records()
        matching = [r for name, r in records.items() if name.startswith(f"ID={transcript_id};")]
        self.assertEqual(len(matching), 1, f"expected one record for {transcript_id}")
        return matching[0]

    def test_a_non_coding_transcript_does_not_abort_the_run(self):
        """The regression: one lncRNA transcript raised KeyError and no file was written."""
        self.assertIn(NONCODING_TRANSCRIPT, " ".join(self.records()))

    def test_the_coding_transcript_carries_its_codon_bounds(self):
        """The silent half: these were left at the transcript start by the chained write."""
        record = self.transcript_record(CODING_TRANSCRIPT)
        self.assertEqual((int(record[6]), int(record[7])), self.CODING_THICK)

    def test_the_non_coding_transcript_marks_no_coding_range(self):
        record = self.transcript_record(NONCODING_TRANSCRIPT)
        self.assertEqual(int(record[6]), self.NONCODING_START)
        self.assertEqual(int(record[7]), self.NONCODING_START)

    def test_every_record_satisfies_the_bed12_invariants(self):
        """start <= thickStart <= thickEnd <= end, and the last block ends at end - start.

        Run over the gene records too, which the thick override never touches.
        """
        for name, record in sorted(self.records().items()):
            with self.subTest(name=name):
                start, end = int(record[1]), int(record[2])
                thick_start, thick_end = int(record[6]), int(record[7])
                self.assertLessEqual(start, thick_start)
                self.assertLessEqual(thick_start, thick_end)
                self.assertLessEqual(thick_end, end)

                sizes = [int(v) for v in record[10].rstrip(",").split(",")]
                starts = [int(v) for v in record[11].rstrip(",").split(",")]
                self.assertEqual(int(record[9]), len(sizes))
                self.assertEqual(len(sizes), len(starts))
                self.assertEqual(starts[0], 0)
                self.assertEqual(starts[-1] + sizes[-1], end - start)

    def test_a_duplicated_transcript_id_is_an_error(self):
        """The uniqueness gtf2bed's reindex depends on, and which would mis-pair silently.

        A second transcript row under an id that already has one would otherwise align one
        transcript's coding range onto another.
        """
        directory = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(directory, ignore_errors=True))
        gtf_path = write_mixed(directory)
        with io.open(gtf_path, "a") as f:
            f.write(NONCODING_LINES[1] + "\n")
        result = run_subcommand("gtf2bed", [gtf_path], directory)
        self.assertNotEqual(result.returncode, 0, "a duplicated transcript_id was accepted")
        self.assertIn("more than one", result.stderr)


# %%
# The exon_id fallback, for a GTF that has none

class Gtf2bedWithoutExonIdTests(unittest.TestCase):
    """gtf2bed supplies synthetic exon ids when the GTF carries none, and that branch
    called `np.unicode` -- an alias numpy 2 removed. Any GTF without exon_id therefore
    raised AttributeError instead of taking the fallback. GENCODE always carries exon_id,
    so nothing in the repository reached this line.
    """

    def test_a_gtf_without_exon_id_still_converts(self):
        directory = tempfile.mkdtemp()
        self.addCleanup(lambda: __import__("shutil").rmtree(directory, ignore_errors=True))
        gtf_path = write_without_exon_id(directory)
        result = run_subcommand("gtf2bed", [gtf_path], directory)
        self.assertEqual(
            result.returncode,
            0,
            f"gtf2bed failed on a GTF without exon_id:\n{result.stderr}",
        )
        with io.open(os.path.join(directory, f"{FIXTURE_STEM}.bed")) as f:
            rows = [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]

        transcript = [r for r in rows if r[3].startswith("ID=TX1;")]
        self.assertEqual(len(transcript), 1, f"expected one transcript record, got {rows}")
        record = transcript[0]
        # The two exons must still be counted, which is what the synthetic ids are for.
        self.assertEqual(record[9], "2")
        self.assertEqual(record[10], "101,101,")
        self.assertEqual(record[11], "0,300,")


if __name__ == "__main__":
    unittest.main(verbosity=2)
