#!/usr/bin/bash

# Usage:
#   align_star.sh <index> <output_prefix> <r1> [<r2>]
#


#
# Subs
#
THREADS=16

show_help() {
  sed -n '2,/^$/p' "$0"  | sed 's/^# \?//'
  exit 0
}

align_star () {
  cmd=(
    STAR
    --readFilesCommand zcat
    --outSAMtype BAM SortedByCoordinate
    --outSAMstrandField intronMotif
    --outSAMattributes NH HI AS nM NM ch
    --outSAMunmapped Within
    --outFilterType BySJout
    --outFilterMultimapNmax 20
    --alignSJoverhangMin 8
    --alignSJDBoverhangMin 1
    --outFilterMismatchNmax 999
    --outFilterMismatchNoverReadLmax 0.04
    --alignIntronMin 20
    --alignIntronMax 1000000
    --alignMatesGapMax 1000000
    --quantMode TranscriptomeSAM
    --twopassMode Basic
    --runThreadN "$THREADS"
    --genomeDir "$index"
    --outFileNamePrefix "$output_prefix"
    --readFilesIn "${reads[@]}"
  )

  printf '%q ' "${cmd[@]}" && echo
  "${cmd[@]}"
}


#
# Main
#
if [[ $# -eq 0 || "$1" == "-h" ]]; then
  show_help
fi

index=$1
output_prefix=$2

shift 2
reads=( "$@" )

output_dir=$(dirname "$output_prefix")
if [ ! -e $output_dir ]; then
  mkdir -p "$output_dir"
fi

STAR --version
align_star
