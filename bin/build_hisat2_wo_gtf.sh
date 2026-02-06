#! /usr/bin/bash

#
# Build index for HISAT2 without GTF
#
# Usage:
#   this.sh <transcript.fasta>
#

#
# Subs
#
show_help() {
  sed -n '2,/^$/p' "$0"  | sed 's/^# \?//'
  exit 0
}

#
# Main
#
if [[ $# -eq 0 || "$1" == "-h" ]]; then
  show_help
fi

fasta=$1

output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa)

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

hisat2-build -p 4 ${fasta} ${output_dir}/${output_base}
