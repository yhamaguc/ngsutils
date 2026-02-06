#! /usr/bin/bash

#
# Build transcriptome index for RSEM
#   and generate Ng vector for EBSeq
#
# Usage:
#   qsub this.sh <gene2tx> <transcripts.fasta>
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

gene2tx=$1
fasta=$2
output_dir=.
output_base=$(basename $(basename ${fasta} .fasta) .fa)

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

# Build index
rsem-prepare-reference --transcript-to-gene-map ${gene2tx} ${fasta} rsem/${output_base}

# Generate ngvector for EBSeq
# NOTE: Rscript execution in the container cannot find EBSeq library included in RSEM

if [ ! -e ebseq ]; then
  mkdir -p ebseq
fi

rsem-generate-ngvector ${fasta} ebseq/${output_base}.tmp

# Combine transcript name and ngvector
grep ">" "${fasta}" > "ebseq/${output_base}.tmp.names"

sed -e 's/^>//' "ebseq/${output_base}.tmp.names" > "ebseq/${output_base}.tmp.names.trimmed"

paste "ebseq/${output_base}.tmp.names.trimmed" "ebseq/${output_base}.tmp.ngvec" > "ebseq/${output_base}.ngvec"

rm "ebseq/${output_base}.tmp.names" "ebseq/${output_base}.tmp.names.trimmed" "ebseq/${output_base}.tmp.ngvec"
