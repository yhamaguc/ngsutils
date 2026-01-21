#! /usr/bin/bash
#
# Prepare index for HISAT2
#
# Usage:
#   qsub this.sh <genome.fasta> <gene-annotation.gtf>
#
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -l s_vmem=48G -l mem_req=48G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

fasta=$1
gtf=$2

output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa).$(basename ${gtf} .gtf)
output_splicesites=${output_dir}/${output_base}.splicesites.txt
output_exons=${output_dir}/${output_base}.exons.txt

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

hisat2_extract_splice_sites.py ${gtf} > ${output_splicesites}
hisat2_extract_exons.py ${gtf} > ${output_exons}
hisat2-build -p 8 --ss ${output_splicesites} --exon ${output_exons} ${fasta} ${output_dir}/${output_base}
