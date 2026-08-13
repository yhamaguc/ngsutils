#!/usr/bin/env bash

#
# Usage:
#   build_rsem_ebseq.sh --gene2tx=gene2tx.tsv transcripts.fa
#
#   find ref -name '*.transcripts.fa' \
#     | xargs -P2 -n1 build_rsem_ebseq.sh --gene2tx=gene2tx.tsv --output-dir=index
#

DOC="Build a transcriptome index for RSEM and an Ng vector for EBSeq

Usage:
  build_rsem_ebseq.sh --gene2tx=<PATH> [--output-dir=<PATH>] <fasta>
  build_rsem_ebseq.sh (-h | --help)

Arguments:
  <fasta>  Transcript FASTA

Options:
  --gene2tx=<PATH>     Transcript-to-gene map passed to rsem-prepare-reference
  --output-dir=<PATH>  Output directory; 'rsem' and 'ebseq' are created inside
                       it [default: .]
  -h --help            Show this message
"

#
# Main
#
if ! command -v docopts > /dev/null; then
  echo "Error: docopts not found on PATH. bin/*.sh parse their arguments with" >&2
  echo "       it; install it from https://github.com/docopt/docopts" >&2
  exit 127
fi

# NOTE: docopts emits shell code to eval -- the parsed values on success, an
#   'echo usage >&2; exit 64' on a bad command line, 'echo help; exit 0' for
#   --help. Empty output means it did none of those (a malformed DOC makes it
#   panic), so nothing may be assumed parsed.
parsed_=$(docopts -A args -h "${DOC}" : "$@")
if [ -z "${parsed_}" ]; then
  echo "Error: docopts could not parse the usage of $(basename "$0");" >&2
  echo "       its own message is above." >&2
  exit 1
fi
eval "${parsed_}"

fasta=${args[<fasta>]}
output_dir=${args[--output-dir]}

output_base=$(basename "$(basename "${fasta}" .fasta)" .fa)
rsem_dir=${output_dir}/rsem
ebseq_dir=${output_dir}/ebseq

if [ ! -e "${rsem_dir}" ]; then
  mkdir -p "${rsem_dir}"
fi

if [ ! -e "${ebseq_dir}" ]; then
  mkdir -p "${ebseq_dir}"
fi

#
# Build index
#
rsem_=(
  rsem-prepare-reference
  --transcript-to-gene-map "${args[--gene2tx]}"
  "${fasta}"
  "${rsem_dir}/${output_base}"
)

echo "CMD: ${rsem_[*]}"
"${rsem_[@]}"

#
# Generate ngvector for EBSeq
#
# NOTE: Rscript execution in the container cannot find EBSeq library included in RSEM
ngvec_=(
  rsem-generate-ngvector
  "${fasta}"
  "${ebseq_dir}/${output_base}.tmp"
)

echo "CMD: ${ngvec_[*]}"
if ! "${ngvec_[@]}"; then
  echo "Error: rsem-generate-ngvector failed on ${fasta}" >&2
  exit 1
fi

#
# Combine transcript name and ngvector
#
# NOTE: the paste below is a positional join -- row i of the name list is
#   assumed to describe row i of the ngvector. Nothing downstream can detect a
#   length mismatch, so it is checked here, at the one step that pairs them.
grep ">" "${fasta}" \
  | sed -e 's/^>//' \
  > "${ebseq_dir}/${output_base}.tmp.names"

n_names=$(wc -l < "${ebseq_dir}/${output_base}.tmp.names")
n_ngvec=$(wc -l < "${ebseq_dir}/${output_base}.tmp.ngvec")

if [ "${n_names}" -ne "${n_ngvec}" ]; then
  echo "Error: ${n_names} transcript names but ${n_ngvec} ngvector rows for ${fasta};" >&2
  echo "       these are joined by position, so the mismatch is not recoverable." >&2
  exit 1
fi

paste \
  "${ebseq_dir}/${output_base}.tmp.names" \
  "${ebseq_dir}/${output_base}.tmp.ngvec" \
  > "${ebseq_dir}/${output_base}.ngvec"

rm \
  "${ebseq_dir}/${output_base}.tmp.names" \
  "${ebseq_dir}/${output_base}.tmp.ngvec"
