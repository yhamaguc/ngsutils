#! /usr/bin/bash

#
# Prepare index for STAR without annotation
#
# Usage:
#   qsub this.sh <genome.fasta>
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

mkdir -p ${output_dir}/${output_base}

STAR --version
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${output_dir}/${output_base} --genomeFastaFiles ${fasta}
