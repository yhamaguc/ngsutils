#!/usr/bin/env bash

#
# Usage:
#   build_hisat2.sh --gtf=gencode.v50.annotation.gtf GRCh38.fa
#
#   find ref -name '*.fa' \
#     | xargs -P2 -n1 build_hisat2.sh --gtf=gencode.v50.annotation.gtf --output-dir=index
#

DOC="Build a HISAT2 index with splice sites and exons from a GTF

Usage:
  build_hisat2.sh --gtf=<PATH> [--output-dir=<PATH>] [--threads=<n>] <fasta>
  build_hisat2.sh (-h | --help)

Arguments:
  <fasta>  Genome FASTA

Options:
  --gtf=<PATH>         Gene annotation GTF
  --output-dir=<PATH>  Output directory [default: .]
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
gtf=${args[--gtf]}
output_dir=${args[--output-dir]}

output_base=$(basename "$(basename "${fasta}" .fasta)" .fa).$(basename "${gtf}" .gtf)
output_splicesites=${output_dir}/${output_base}.splicesites.txt
output_exons=${output_dir}/${output_base}.exons.txt

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

# NOTE: hisat2-build below reads both of these files. An empty splice-site or
#   exon table is silently accepted there and yields an index with no annotation
#   in it, which is indistinguishable from build_hisat2_wo_gtf.sh output. The
#   extraction status is therefore checked here rather than left to the caller.
if ! hisat2_extract_splice_sites.py "${gtf}" > "${output_splicesites}"; then
  echo "Error: hisat2_extract_splice_sites.py failed on ${gtf}" >&2
  exit 1
fi

if ! hisat2_extract_exons.py "${gtf}" > "${output_exons}"; then
  echo "Error: hisat2_extract_exons.py failed on ${gtf}" >&2
  exit 1
fi

cmd_=(
  hisat2-build
  -p "${args[--threads]}"
  --ss "${output_splicesites}"
  --exon "${output_exons}"
  "${fasta}"
  "${output_dir}/${output_base}"
)

echo "CMD: ${cmd_[*]}"
"${cmd_[@]}"
