#! /usr/bin/bash
#
# Usage:
#   gtf2bed.sh <gtf>
#

input=${1}
input_root=$(basename ${input%.*})

cmd="gtfToGenePred ${input} ${input_root}.genePred"
echo $cmd
eval $cmd

cmd="genePredToBed ${input_root}.genePred ${input_root}.bed12"
echo $cmd
eval $cmd

cmd="sort -k1,1 -k2,2n ${input_root}.bed12 > ${input_root}.bed"
echo $cmd
eval $cmd

rm ${input_root}.genePred ${input_root}.bed12
