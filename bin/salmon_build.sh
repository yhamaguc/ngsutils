#! /usr/bin/bash
#
# Prepare transcriptome index for salmon
#
# Usage:
#   qsub this.sh <fasta>
#
#$ -S /bin/bash
#$ -l s_vmem=16G -l mem_req=16G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

fasta=$1
output=$(basename $(basename ${fasta} .fasta) .fa)
output_dir=.

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

cmd="salmon index -t ${fasta} -i ${output_dir}/${output}"
echo $cmd
$cmd
