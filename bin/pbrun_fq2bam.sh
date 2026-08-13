#!/usr/bin/env bash

#
# Usage:
#   pbrun_fq2bam.sh --index=GRCh38.fa --output-dir=out r1.fq.gz r2.fq.gz
#
#   find fastq -name '*_r1.fastq.gz' \
#     | xargs -P1 -n1 -I{} pbrun_fq2bam.sh --index=GRCh38.fa --output-dir=out/{} {}
#

DOC="Align reads and recalibrate base qualities with Parabricks fq2bam

Usage:
  pbrun_fq2bam.sh --index=<PATH> [--output-dir=<PATH>] <r1> [<r2>]
  pbrun_fq2bam.sh (-h | --help)

Arguments:
  <r1>  Read 1 FASTQ
  <r2>  Read 2 FASTQ; omit for single-end

Options:
  --index=<PATH>       Reference FASTA passed to --ref
  --output-dir=<PATH>  Output directory [default: .]
  -h --help            Show this message
"

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

output_dir=${args[--output-dir]}

# NOTE: --in-fq takes one or two files. This used to pass "$3" "$4"
#   unconditionally, so a single-end run handed pbrun an empty second argument.
reads=("${args[<r1>]}")
if [ -n "${args[<r2>]}" ]; then
  reads+=("${args[<r2>]}")
fi

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

cmd_=(
  pbrun fq2bam
  --ref "${args[--index]}"
  --in-fq "${reads[@]}"
  --out-bam "${output_dir}/aligned.sorted.bam"
  --out-recal-file "${output_dir}/recal_data.table"
)

echo "CMD: ${cmd_[*]}"
"${cmd_[@]}"
