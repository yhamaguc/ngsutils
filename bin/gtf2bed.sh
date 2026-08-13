#!/usr/bin/env bash

#
# Usage:
#   gtf2bed.sh gencode.v50.annotation.gtf
#
#   find ref -name '*.gtf' \
#     | xargs -P4 -n1 gtf2bed.sh --output-dir=bed
#

DOC="Convert a GTF to a coordinate-sorted BED12 using the UCSC utilities

Usage:
  gtf2bed.sh [--output-dir=<PATH>] <gtf>
  gtf2bed.sh (-h | --help)

Arguments:
  <gtf>  Gene annotation GTF

Options:
  --output-dir=<PATH>  Output directory [default: .]
  -h --help            Show this message
"

#
# Subs
#
run_() {
  echo "CMD: $*"
  "$@"
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

gtf=${args[<gtf>]}
output_dir=${args[--output-dir]}

input_root=$(basename "${gtf%.*}")

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

output_genepred=${output_dir}/${input_root}.genePred
output_bed12=${output_dir}/${input_root}.bed12
output_bed=${output_dir}/${input_root}.bed

# NOTE: the two intermediates are removed at the end, so every step is checked
#   here -- a failure that reached the rm would delete the only evidence of how
#   far the conversion got and leave a truncated .bed looking finished.
if ! run_ gtfToGenePred "${gtf}" "${output_genepred}"; then
  echo "Error: gtfToGenePred failed on ${gtf}" >&2
  exit 1
fi

if ! run_ genePredToBed "${output_genepred}" "${output_bed12}"; then
  echo "Error: genePredToBed failed on ${output_genepred}" >&2
  exit 1
fi

# NOTE: gtfToGenePred and genePredToBed both emit records in GTF order; the BED
#   consumers here expect coordinate order, so the sort is not optional.
echo "CMD: sort -k1,1 -k2,2n ${output_bed12} > ${output_bed}"
if ! sort -k1,1 -k2,2n "${output_bed12}" > "${output_bed}"; then
  echo "Error: sort failed on ${output_bed12}" >&2
  exit 1
fi

rm "${output_genepred}" "${output_bed12}"
