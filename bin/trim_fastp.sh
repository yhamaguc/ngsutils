#!/usr/bin/env bash

#
# Usage:
#   trim_fastp.sh --output-dir=trimmed sample_r1.fastq.gz sample_r2.fastq.gz
#
#   find fastq -name '*_r1.fastq.gz' \
#     | xargs -P4 -n1 trim_fastp.sh --output-dir=trimmed
#

DOC="Trim reads with fastp

Usage:
  trim_fastp.sh [--output-dir=<PATH>] [--threads=<n>] [--length-required=<n>]
                <r1> [<r2>]
  trim_fastp.sh (-h | --help)

Arguments:
  <r1>  Read 1 FASTQ, gzip compressed
  <r2>  Read 2 FASTQ, gzip compressed; omit for single-end

Options:
  --output-dir=<PATH>    Output directory [default: .]
  --threads=<n>          Threads [default: 8]
  --length-required=<n>  Discard reads shorter than this after trimming
                         [default: 31]
  -h --help              Show this message
"

#
# Subs
#
output_prefix_from_() {
  local name
  name=$(basename "$1")

  # NOTE: strip the read-1 marker first so both mates of a pair report under one
  #   prefix. The second substitution is what gives the single-end path a prefix
  #   at all -- it was missing until 2026-08-14, and every single-end run wrote
  #   its reports to the same bare '.fastp.json' and '.fastp.html'.
  name=$(printf '%s' "${name}" | sed -E 's/_r1\.fastq\.gz$//I')
  name=$(printf '%s' "${name}" | sed -E 's/\.f(ast)?q\.gz$//I')

  printf '%s' "${name}"
}


trim_pe_() {
  local cmd_=("${base_command_[@]}")

  cmd_+=(
    --in1 "${r1}"
    --in2 "${r2}"
    --out1 "${output_dir}/$(basename "${r1}")"
    --out2 "${output_dir}/$(basename "${r2}")"
    --json "${output_dir}/${output_prefix}.fastp.json"
    --html "${output_dir}/${output_prefix}.fastp.html"
    --correction
    --detect_adapter_for_pe
  )

  echo "CMD: ${cmd_[*]}"
  "${cmd_[@]}"
}


trim_sr_() {
  local cmd_=("${base_command_[@]}")

  cmd_+=(
    --in1 "${r1}"
    --out1 "${output_dir}/$(basename "${r1}")"
    --json "${output_dir}/${output_prefix}.fastp.json"
    --html "${output_dir}/${output_prefix}.fastp.html"
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

r1=${args[<r1>]}
r2=${args[<r2>]}
output_dir=${args[--output-dir]}

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

base_command_=(
  fastp
  --thread "${args[--threads]}"
  --trim_tail1 1
  --trim_tail2 1
  --length_required "${args[--length-required]}"
  --overrepresentation_analysis
)

output_prefix=$(output_prefix_from_ "${r1}")

if [ -n "${r2}" ]; then
  trim_pe_
else
  trim_sr_
fi
