#! /usr/bin/bash

#
# Prepare index for STAR
#
# Usage:
#   qsub this.sh <genome.fasta> <gene-annotation.gtf>
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
gtf=$2
output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa).$(basename ${gtf} .gtf)

mkdir -p ${output_dir}/${output_base}

STAR --version
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${output_dir}/${output_basel} --genomeFastaFiles "${fasta}" --sjdbGTFfile ${gtf}