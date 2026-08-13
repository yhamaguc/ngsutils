#!/usr/bin/env bash

#
# Usage:
#   build_star_wo_gtf.sh GRCh38.fa
#
#   find ref -name '*.fa' \
#     | xargs -P2 -n1 build_star_wo_gtf.sh --output-dir=index
#

DOC="Build a STAR genome index without a gene annotation

Usage:
  build_star_wo_gtf.sh [--output-dir=<PATH>] [--threads=<n>] <fasta>
  build_star_wo_gtf.sh (-h | --help)

Arguments:
  <fasta>  Genome FASTA

Options:
  --output-dir=<PATH>  Parent of the generated genome directory [default: .]
  --threads=<n>        Threads [default: 8]
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

output_base=$(basename "$(basename "${fasta}" .fasta)" .fa)
genome_dir=${args[--output-dir]}/${output_base}

mkdir -p "${genome_dir}"

STAR --version

cmd_=(
  STAR
  --runThreadN "${args[--threads]}"
  --runMode genomeGenerate
  --genomeDir "${genome_dir}"
  --genomeFastaFiles "${fasta}"
)

echo "CMD: ${cmd_[*]}"
"${cmd_[@]}"
