#! /usr/bin/bash

#
# Usage:
#   this.sh <index> <output_dir> <r1> [<r2>]
#

#
# Subs
#
show_help() {
  sed -n '2,/^$/p' "$0"  | sed 's/^# \?//'
  exit 0
}

index=$1
output_dir=$2

if [ ! -e "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

cmd=(
  pbrun fq2bam
  --ref "$index"
  --in-fq "$3" "$4"
  --out-bam "$output_dir/aligned.sorted.bam"
  --out-recal-file "$output_dir/recal_data.table"
  )

printf '%q ' "${cmd[@]}" && echo
"${cmd[@]}"
