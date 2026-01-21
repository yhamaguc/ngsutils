#! /usr/bin/bash
#
# Prepare index for HISAT2
#
# Usage:
#   qsub this.sh <transcript.fasta>
#
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -l s_vmem=48G -l mem_req=48G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

fasta=$1

output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa)

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

hisat2-build -p 4 ${fasta} ${output_dir}/${output_base}
