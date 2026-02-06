#! /usr/bin/bash

#
# Prepare transcriptome index for salmon
#
# Usage:
#   qsub this.sh <fasta>
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
output=$(basename $(basename ${fasta} .fasta) .fa)
output_dir=.

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

salmon index -t ${fasta} -i ${output_dir}/${output}
