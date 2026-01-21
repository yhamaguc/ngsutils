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
  rsem/${output_base}"
echo ${cmd}
eval ${cmd}

# Generate ngvector for EBSeq
# NOTE: Rscript execution in the container cannot find EBSeq library　included in RSEM

if [ ! -e ebseq ]; then
  mkdir -p ebseq
fi

cmd="rsem-generate-ngvector \
  ${fasta} \
  ebseq/${output_base}.tmp"
echo ${cmd}
eval ${cmd}

# Combine transcript name and ngvector
cmd="grep \">\" ${fasta} > ebseq/${output_base}.tmp.names"
echo ${cmd}
eval ${cmd}

cmd="sed -e 's/^>//' ebseq/${output_base}.tmp.names > ebseq/${output_base}.tmp.names.trimmed"
echo ${cmd}
eval ${cmd}

cmd="paste ebseq/${output_base}.tmp.names.trimmed ebseq/${output_base}.tmp.ngvec > ebseq/${output_base}.ngvec"
echo ${cmd}
eval ${cmd}

rm ebseq/${output_base}.tmp.names
rm ebseq/${output_base}.tmp.names.trimmed
rm ebseq/${output_base}.tmp.ngvec
