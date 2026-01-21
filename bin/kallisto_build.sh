#! /usr/bin/bash
#
# Prepare transcriptome index for kallisto
#
# Usage:
#   qsub this.sh <transcripts.fasta>
#
#$ -S /bin/bash
#$ -l s_vmem=32G -l mem_req=32G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

fasta=$1
output_dir=.
output=$(basename $(basename $(basename ${fasta} .fasta) .fa) .fna)

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

cmd="kallisto index -i ${output_dir}/${output} ${fasta}"

echo $cmd
$cmd

