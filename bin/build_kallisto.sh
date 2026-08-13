#!/usr/bin/env bash

#
# Usage:
#   build_kallisto.sh gencode.v50.transcripts.fa
#
#   find ref -name '*.transcripts.fa' \
#     | xargs -P4 -n1 build_kallisto.sh --output-dir=index
#

DOC="Build a transcriptome index for kallisto

Usage:
  build_kallisto.sh [--output-dir=<PATH>] <fasta>
  build_kallisto.sh (-h | --help)

Arguments:
  <fasta>  Transcript FASTA

Options:
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

fasta=${args[<fasta>]}
output_dir=${args[--output-dir]}

output_base=$(basename "$(basename "$(basename "${fasta}" .fasta)" .fa)" .fna)

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

cmd_=(
  kallisto index
  -i "${output_dir}/${output_base}"
  "${fasta}"
)

echo "CMD: ${cmd_[*]}"
"${cmd_[@]}"
