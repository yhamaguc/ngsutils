#! /usr/bin/bash
#
# Prepare transcriptome index for RSEM
#   and generate Ng vector for EBSeq
#
# Usage:
#   qsub this.sh <gene2tx> <transcripts.fasta>
#
#$ -S /bin/bash
#$ -l s_vmem=32G -l mem_req=32G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

gene2tx=$1
fasta=$2
output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa)

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

# Build index
cmd="rsem-prepare-reference \
  --transcript-to-gene-map  ${gene2tx} \
  ${fasta} \
  ${output_dir}/${output_base}"
echo ${cmd}
eval ${cmd}
