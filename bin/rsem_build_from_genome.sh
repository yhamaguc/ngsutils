#! /usr/bin/bash
#
# Prepare transcriptome index for RSEM
#
# Usage:
#   qsub this.sh <genome.fasta> <gene-annotation.gtf>
#
#$ -S /bin/bash
#$ -l s_vmem=16G -l mem_req=16G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

genome=$1
annotation=$2
output_dir=.
output=$(basename $(basename ${genome} .fasta) .fa).$(basename ${annotation} .gtf)

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

cmd="rsem-prepare-reference \
    --gtf ${annotation} \
    ${genome} \
    ${output_dir}/${output}"

echo $cmd
$cmd
