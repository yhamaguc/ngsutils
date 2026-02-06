#!/bin/bash

#
# Quantify transcript abundance using Salmon
#
# Usage:
#   quant_salmon.sh <index> <libtype> <output_dir> <input1> [<input2>]
#
# Arguments:
#   <index>        Salmon index directory
#   <libtype>      Library type (e.g., ISR, IU, A)
#   <output_dir>   Output directory
#   <input1>       Input file (BAM or FASTQ R1). File type determined by extension
#   <input2>       Read 2 FASTQ file (optional, for paired-end FASTQ)
#
# Examples:
#   # BAM input (auto-detected)
#   quant_salmon.sh index.idx ISR output input.bam
#
#   # FASTQ paired-end
#   quant_salmon.sh index.idx ISR output r1.fq.gz r2.fq.gz
#
#   # FASTQ single-end
#   quant_salmon.sh index.idx U output reads.fq.gz
#

set -euo pipefail

#
# Subs
#
show_help() {
  sed -n '2,/^$/p' "$0" | sed 's/^# \?//'
  exit 0
}


#
# Main
#
if [[ $# -lt 4 || "$1" == "-h" || "$1" == "--help" ]]; then
  show_help
fi

THREADS=8
index="$1"
libtype="$2"
output_dir="$3"
input1="$4"
input2="${5:-}"

mkdir -p "$output_dir"

salmon_base=(
  salmon quant
  -p "$THREADS"
  -i "$index"
  -l "$libtype"
)

if [[ "$input1" =~ \.bam$ ]]; then
  filebase=$(basename "$input1" .bam)
  tmp_dir="${output_dir}/tmp_fastq"
  mkdir -p "$tmp_dir"

  bamtofastq_cmd=(
    bamtofastq
    gz=1
    fastq=0
    disablevalidation=1
    "S=${tmp_dir}/${filebase}.sr.fastq.gz"
    "F=${tmp_dir}/${filebase}.r1.fastq.gz"
    "F2=${tmp_dir}/${filebase}.r2.fastq.gz"
    O=/dev/null
    O2=/dev/null
    "filename=$input1"
  )

  printf '%q ' "${bamtofastq_cmd[@]}" && echo
  "${bamtofastq_cmd[@]}"

  if [[ "${libtype:0:1}" == "I" ]]; then
    r1="${tmp_dir}/${filebase}.r1.fastq.gz"
    r2="${tmp_dir}/${filebase}.r2.fastq.gz"
  else
    sr="${tmp_dir}/${filebase}.sr.fastq.gz"
  fi
else
  r1="$input1"
  r2="$input2"
fi

salmon_cmd=("${salmon_base[@]}")

if [[ -n "${r2:-}" ]]; then
  salmon_cmd+=(-1 "$r1" -2 "$r2")
elif [[ -n "${sr:-}" ]]; then
  salmon_cmd+=(-r "$sr")
else
  salmon_cmd+=(-1 "$r1")
fi

salmon_cmd+=(-o "$output_dir")

printf '%q ' "${salmon_cmd[@]}" && echo
"${salmon_cmd[@]}"

if [[ "$input1" =~ \.bam$ ]]; then
  rm -rf "$tmp_dir"
fi
