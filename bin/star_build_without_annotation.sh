#! /usr/bin/bash
#
# Prepare index for STAR without annotation
#
# Usage:
#   qsub this.sh <genome.fasta>
#
#$ -S /bin/bash
#$ -pe def_slot 8
#$ -l s_vmem=8G -l mem_req=8G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

fasta=$1

output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa)

mkdir -p ${output_dir}/${output_base}

STAR --version
cmd="STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${output_dir}/${output_base} --genomeFastaFiles ${fasta}"
echo $cmd
$cmd
