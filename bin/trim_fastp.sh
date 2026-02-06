#! /usr/bin/bash

# Usage:
#   trim_fastp.sh <output_dir> <r1> [<r2>]
#


#
# Subs
#
show_help() {
  sed -n '2,/^$/p' "$0"  | sed 's/^# \?//'
  exit 0
}

THREADS=8
LENGTH_REQUIRED=31

base_command=(
  fastp
    --thread $THREADS
    --trim_tail1 1
    --trim_tail2 1
    --length_required $LENGTH_REQUIRED
    --overrepresentation_analysis
)


trim_pe() {
  local cmd=("${base_command[@]}")
  output_prefix=$(basename "$in1" | sed -E 's/_r1\.fastq\.gz$//I')

  cmd+=(
    --in1 "$in1"
    --in2 "$in2"
    --out1 "${output_dir}/$(basename "$in1")"
    --out2 "${output_dir}/$(basename "$in2")"
    --json "${output_dir}/$(basename "${output_prefix}.fastp.json")"
    --html "${output_dir}/$(basename "${output_prefix}.fastp.html")"
    --correction
    --detect_adapter_for_pe
  )

  printf '%q ' "${cmd[@]}" && echo
  "${cmd[@]}"
}


trim_sr() {
  local cmd=("${base_command[@]}")

  cmd+=(
    --in1 "$in1"
    --out1 "${output_dir}/$(basename "$in1")"
    --json "${output_dir}/$(basename "${output_prefix}.fastp.json")"
    --html "${output_dir}/$(basename "${output_prefix}.fastp.html")"
  )

  printf '%q ' "${cmd[@]}" && echo
  "${cmd[@]}"
}


#
# Main
#
if [[ $# -eq 0 || "$1" == "-h" ]]; then
  show_help
fi

if [ $# -lt 2 ]; then
  echo "Usage (PE): $0 <r1> <r2> <output_dir>"
  echo "Usage (SR): $0 <r1> <output_dir>"
  exit 1
fi

output_dir=$1
mkdir -p "$output_dir"

if [ $# -eq 3 ]; then
  in1=$2
  in2=$3
elif [ $# -eq 2 ]; then
  in1=$2
  in2=""
fi

if [ -z "$in2" ]; then
  trim_sr
else
  trim_pe
fi
