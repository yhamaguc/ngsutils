#! /usr/bin/bash
#
# Prepare index for STAR
#
# Usage:
#   qsub this.sh <genome.fasta> <gene-annotation.gtf>
#
#$ -S /bin/bash
#$ -pe def_slot 8
#$ -l s_vmem=8G -l mem_req=8G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

fasta=$1
gtf=$2
output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa).$(basename ${gtf} .gtf)

mkdir -p ${output_dir}/${output_base}

STAR --version
cmd="STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${output_dir}/${output_base} --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf}"
echo $cmd && eval $cmd
