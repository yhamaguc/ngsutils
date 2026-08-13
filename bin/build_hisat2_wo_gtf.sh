#!/usr/bin/env bash

#
# Usage:
#   build_hisat2_wo_gtf.sh transcripts.fa
#
#   find ref -name '*.fa' \
#     | xargs -P4 -n1 build_hisat2_wo_gtf.sh --output-dir=index
#

DOC="Build a HISAT2 index without a gene annotation

Usage:
  build_hisat2_wo_gtf.sh [--output-dir=<PATH>] [--threads=<n>] <fasta>
  build_hisat2_wo_gtf.sh (-h | --help)

Arguments:
  <fasta>  Genome or transcript FASTA

Options:
  --output-dir=<PATH>  Output directory [default: .]
  --threads=<n>        Threads [default: 4]
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

fasta=${args[<fasta>]}
output_dir=${args[--output-dir]}

output_base=$(basename "$(basename "${fasta}" .fasta)" .fa)

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

cmd_=(
  hisat2-build
  -p "${args[--threads]}"
  "${fasta}"
  "${output_dir}/${output_base}"
)

echo "CMD: ${cmd_[*]}"
"${cmd_[@]}"
