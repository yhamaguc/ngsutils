#!/usr/bin/env python3
"""Coordinate-sort a STAR unsorted BAM with samtools, verifying that no reads are lost.

Why this is a separate step, and not part of align_star
-------------------------------------------------------
STAR's own sorter is in-memory and bounded by --limitBAMsortRAM. The requirement scales
with read count, so ANY fixed ceiling is one large sample away from failing -- and it
fails only AFTER the whole mapping phase has been paid for:

    EXITING because of fatal ERROR: not enough memory for BAM sorting:
    SOLUTION: re-run STAR with at least --limitBAMsortRAM 39842852549

Observed history on real data: the default (= index size, ~30 GB) failed on samples
needing ~40 GB; a 64 GB ceiling failed on a sample that asked for 132,227,075,545 bytes
after 83 min of mapping. samtools sort spills to disk, so the ceiling becomes disk space
rather than RAM and a large sample costs time instead of a crash.

Keeping it a separate step means the memory-hungry half can be retried, given different
resources, or run on a different host without repeating the alignment. Pair it with
`ngsutils align_star --sort-with none`, which stops after writing Aligned.out.bam.

Output naming
-------------
The sample is identified by the DIRECTORY, not the filename: align_star is given an
out-prefix ending in "/" and STAR writes its fixed basenames inside it, so a run looks
like

    results/align_star/TCGA-UVM/<file_id>/<aliquot>/Aligned.out.bam

The only thing the filename has to carry is that the file is coordinate-sorted. So
<sorted-bam> is optional and derived from the input, keeping its directory and its prefix:

    <anything>Aligned.out.bam  ->  <anything>Aligned.sortedByCoord.out.bam
    <anything>.bam             ->  <anything>.sortedByCoord.bam

The first form is what STAR's own sorter would have written, so nothing downstream has to
know which sorter ran. Pass <sorted-bam> explicitly only to override that.

Usage:
  sort_star_bam [options] <unsorted-bam> [<sorted-bam>]
  sort_star_bam (-h | --help)

Arguments:
  <unsorted-bam>  STAR's Aligned.out.bam.
  <sorted-bam>    Where to write the coordinate-sorted BAM. Derived from the input when
                  omitted; see "Output naming" above.

NOTE: threads is an OPTION here, not a trailing positional as in align_star. Two optional
  positionals in a row are ambiguous and docopt fills them left to right, so
  `sort_star_bam in.bam 8` silently took 8 as the OUTPUT NAME and wrote a file called "8".
  Measured 2026-08-19 on real data, in the first run after the derivation was added. One
  optional positional is the most this command can safely carry.

Options:
  --threads=<threads>           Threads for samtools sort and index. [default: 8]
  --memory-per-thread=<memory>  samtools sort -m. PER THREAD, not total. [default: 4G]
  --keep-unsorted               Keep the input after a verified sort. Off by default:
                                removing it roughly halves the peak footprint per sample.
  --no-index                    Skip samtools index on the result.
  -h --help                     Show this message.
"""

from __future__ import annotations

import shlex
import subprocess
import sys
from pathlib import Path

from docopt import docopt

# %%
# Constants

DEFAULT_THREADS = 8

# 8 x 4G = 32 GB total, against the 132 GB that STAR's in-memory sorter demanded for the
# largest sample seen here -- the point being that samtools does not need to hold the
# whole thing at once.
DEFAULT_MEMORY_PER_THREAD = "4G"

BAM_SUFFIX = ".bam"
TMP_INFIX = ".sorttmp"

# STAR's own names for the two files, so a derived name is indistinguishable from what its
# internal sorter would have produced.
STAR_UNSORTED = "Aligned.out.bam"
STAR_SORTED = "Aligned.sortedByCoord.out.bam"
SORTED_INFIX = ".sortedByCoord"


class SortError(RuntimeError):
    """A refusal or a failure that main() turns into a message and a non-zero exit."""


# %%
# Sorting


def count_records(bam: str) -> int:
    result = subprocess.run(
        ["samtools", "view", "-c", bam], capture_output=True, text=True, check=False
    )
    if result.returncode != 0 or not result.stdout.strip():
        raise SortError(f"could not count records in {bam}: {result.stderr.strip()}")
    return int(result.stdout.strip())


def derive_sorted_name(unsorted: str) -> str:
    """Where the sorted BAM goes when the caller does not say.

    Keeps the directory and the prefix -- both of which carry the sample identity -- and
    changes only the part that says how the file is ordered.
    """
    if unsorted.endswith(STAR_UNSORTED):
        return unsorted[: -len(STAR_UNSORTED)] + STAR_SORTED
    if unsorted.endswith(BAM_SUFFIX):
        return unsorted[: -len(BAM_SUFFIX)] + SORTED_INFIX + BAM_SUFFIX
    return unsorted + SORTED_INFIX + BAM_SUFFIX


def temp_prefix(sorted_bam: str) -> str:
    """Where samtools spills.

    NOTE: -T must sit on a filesystem with room for roughly the size of the BAM, so it
      goes next to the output rather than in /tmp, which is small on these machines.
    """
    return f"{sorted_bam.removesuffix(BAM_SUFFIX)}{TMP_INFIX}"


def clear_temporaries(prefix: str) -> None:
    """Leftover fragments mean a sort died partway. Clear them either way, so a retry
    does not read another run's temporaries."""
    directory = Path(prefix).parent
    for leftover in directory.glob(f"{Path(prefix).name}.*{BAM_SUFFIX}"):
        leftover.unlink(missing_ok=True)


def sort_bam(
    unsorted: str,
    sorted_bam: str | None = None,
    threads: int = DEFAULT_THREADS,
    memory_per_thread: str = DEFAULT_MEMORY_PER_THREAD,
    keep_unsorted: bool = False,
    index: bool = True,
) -> int:
    """Sort, verify, index, and drop the input. Returns the verified record count.

    NOTE: -m is PER THREAD. samtools uses roughly threads x memory and exceeds it
      somewhat, which is the classic way to OOM a machine while believing the limit was
      set. The total is printed so the number is visible in the log.
    NOTE: a sort must not lose reads. samtools exits 0 on truncated input, so the record
      count is compared before and after and a mismatch is fatal. File existence and exit
      status do not distinguish a complete sort from a partial one.
    """
    if not Path(unsorted).exists() or Path(unsorted).stat().st_size == 0:
        raise SortError(f"missing or empty: {unsorted}")

    sorted_bam = sorted_bam or derive_sorted_name(unsorted)
    if sorted_bam == unsorted:
        raise SortError(f"the sorted output would overwrite the input: {unsorted}")

    # NOTE: this guard exists because the failure it catches is destructive. <sorted-bam>
    #   is an optional positional, so a caller reaching for align_star's trailing
    #   [<threads>] and typing `sort_star_bam in.bam 8` hands "8" over as a perfectly valid
    #   output name -- and this command deletes its input once the sort verifies, so the
    #   unsorted BAM is gone and the result is a file called "8". Observed twice on
    #   2026-08-19. A BAM must be named like one.
    if not sorted_bam.endswith(BAM_SUFFIX):
        raise SortError(
            f"the output name does not end in {BAM_SUFFIX}: {sorted_bam!r}. If that was "
            f"meant to be a thread count, it is --threads=<n> here, not a positional."
        )

    Path(sorted_bam).parent.mkdir(parents=True, exist_ok=True)
    prefix = temp_prefix(sorted_bam)

    print(f"Input   : {unsorted}")
    print(f"Output  : {sorted_bam}")
    print(f"Threads : {threads} x {memory_per_thread} per thread")
    print(f"Temp    : {prefix}.*{BAM_SUFFIX}")

    before = count_records(unsorted)
    print(f"Records : {before} in")

    argv = [
        "samtools", "sort",
        "-@", str(threads),
        "-m", memory_per_thread,
        "-T", prefix,
        "-o", sorted_bam,
        unsorted,
    ]
    print(f"CMD: {shlex.join(argv)}")
    status = subprocess.run(argv, check=False).returncode
    clear_temporaries(prefix)

    if status != 0:
        Path(sorted_bam).unlink(missing_ok=True)
        raise SortError(f"samtools sort failed (exit {status})")

    if subprocess.run(["samtools", "quickcheck", sorted_bam], check=False).returncode != 0:
        Path(sorted_bam).unlink(missing_ok=True)
        raise SortError(f"samtools quickcheck rejected {sorted_bam}")

    after = count_records(sorted_bam)
    if before != after:
        Path(sorted_bam).unlink(missing_ok=True)
        raise SortError(
            f"record count changed during sorting: {before} -> {after}. The sorted BAM "
            f"was incomplete and has been removed."
        )

    if index:
        subprocess.run(["samtools", "index", "-@", str(threads), sorted_bam], check=False)
    if not keep_unsorted:
        Path(unsorted).unlink(missing_ok=True)
        print("Removed the unsorted BAM.")

    print(f"Sorted  : {sorted_bam} ({after} records, verified)")
    return after


# %%
# Main


def main() -> int:
    opts = docopt(__doc__)
    try:
        sort_bam(
            opts["<unsorted-bam>"],
            opts["<sorted-bam>"],
            threads=int(opts["--threads"]),
            memory_per_thread=opts["--memory-per-thread"],
            keep_unsorted=opts["--keep-unsorted"],
            index=not opts["--no-index"],
        )
    except SortError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
