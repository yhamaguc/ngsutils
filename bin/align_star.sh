#!/usr/bin/env bash

#
# Usage:
#   align_star.sh --index=star_index --output-prefix=out/sample. r1.fq.gz r2.fq.gz
#
#   find fastq -name '*_r1.fastq.gz' \
#     | xargs -P2 -n1 -I{} align_star.sh --index=star_index --output-prefix=out/{}. {}
#

DOC="Align reads to a genome with STAR

Usage:
  align_star.sh --index=<PATH> --output-prefix=<PATH> [--threads=<n>] <r1> [<r2>]
  align_star.sh (-h | --help)

Arguments:
  <r1>  Read 1 FASTQ, gzip compressed
  <r2>  Read 2 FASTQ, gzip compressed; omit for single-end

Options:
  --index=<PATH>          STAR genome index directory
  --output-prefix=<PATH>  Prefix prepended to every STAR output file
  --threads=<n>           Threads [default: 16]
  -h --help               Show this message
"

#
# Subs
#
align_star() {
  local cmd_=(
    STAR
    --readFilesCommand zcat
    --outSAMtype BAM SortedByCoordinate
    --outSAMstrandField intronMotif
    --outSAMattributes NH HI AS nM NM ch
    --outSAMunmapped Within
    --outFilterType BySJout
    --outFilterMultimapNmax 20
    --alignSJoverhangMin 8
    --alignSJDBoverhangMin 1
    --outFilterMismatchNmax 999
    --outFilterMismatchNoverReadLmax 0.04
    --alignIntronMin 20
    --alignIntronMax 1000000
    --alignMatesGapMax 1000000
    --quantMode TranscriptomeSAM
    --twopassMode Basic
    --runThreadN "${args[--threads]}"
    --genomeDir "${args[--index]}"
    --outFileNamePrefix "${args[--output-prefix]}"
    --readFilesIn "${reads[@]}"
  )

  echo "CMD: ${cmd_[*]}"
  "${cmd_[@]}"
}


#
# Main
#
if ! command -v docopts > /dev/null; then
  echo "Error: docopts not found on PATH. bin/*.sh parse their arguments with" >&2
  echo "       it; install it from https://github.com/docopt/docopts" >&2
  exit 127
fi

# NOTE: docopts emits shell code to eval -- the parsed values on success, an
#   'echo usage >&2; exit 64' on a bad command line, 'echo help; exit 0' for
#   --help. Empty output means it did none of those (a malformed DOC makes it
#   panic), so nothing may be assumed parsed.
parsed_=$(docopts -A args -h "${DOC}" : "$@")
if [ -z "${parsed_}" ]; then
  echo "Error: docopts could not parse the usage of $(basename "$0");" >&2
  echo "       its own message is above." >&2
  exit 1
fi
eval "${parsed_}"

reads=("${args[<r1>]}")
if [ -n "${args[<r2>]}" ]; then
  reads+=("${args[<r2>]}")
fi

output_dir=$(dirname "${args[--output-prefix]}")
if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

STAR --version
align_star
