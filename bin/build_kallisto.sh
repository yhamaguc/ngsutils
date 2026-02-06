#! /usr/bin/bash

#
# Build transcriptome index for kallisto
#
# Usage:
#   qsub this.sh <transcripts.fasta>
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
output=$(basename $(basename $(basename ${fasta} .fasta) .fa) .fna)

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

kallisto index -i ${output_dir}/${output} ${fasta}