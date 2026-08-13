#!/usr/bin/env bash

#
# Usage:
#   build_rsem_ebseq.sh --gene2tx=gene2tx.tsv transcripts.fa.gz
#
#   find ref -name '*.transcripts.fa.gz' \
#     | xargs -P2 -n1 build_rsem_ebseq.sh --gene2tx=gene2tx.tsv --output-dir=index
#

DOC="Build a transcriptome index for RSEM and an Ng vector for EBSeq

Usage:
  build_rsem_ebseq.sh --gene2tx=<PATH> [--output-dir=<PATH>] <fasta>
  build_rsem_ebseq.sh (-h | --help)

Arguments:
  <fasta>  Transcript FASTA, optionally gzipped

Options:
  --gene2tx=<PATH>     Transcript-to-gene map passed to rsem-prepare-reference
  --output-dir=<PATH>  Output directory; 'rsem' and 'ebseq' are created inside
                       it [default: .]
  -h --help            Show this message
"

#
# Subs
#
# NOTE: STAR rejects a gzipped reference outright ('Make sure the file is
#   uncompressed'), and rsem-prepare-reference fails with 'Number of transcripts
#   in the reference is less than 1!', which does not point at the compression
#   at all. hisat2-build does not document .gz support either. salmon and
#   kallisto do read .gz directly, but every build_*.sh expands it the same way
#   regardless -- one behaviour to document beats a per-tool table to keep
#   correct against five tools that each change independently.
tmp_dir=

ungzip_() {
  local src=$1
  local work_dir=$2

  ungzipped=${src}

  case "${src}" in
    *.gz) ;;
    *) return 0 ;;
  esac

  if [ -z "${tmp_dir}" ]; then
    # NOTE: expanded beside the output, not under ${TMPDIR} -- an uncompressed
    #   primary assembly runs to a few GB and /tmp on a compute node rarely
    #   holds one. A private directory also means an existing file of the same
    #   name in the output directory is never overwritten, then deleted by the
    #   cleanup below.
    tmp_dir=$(mktemp -d "${work_dir}/.ngsutils_XXXXXX") || return 1
  fi

  ungzipped=${tmp_dir}/$(basename "${src}" .gz)
  echo "Decompressing ${src} -> ${ungzipped}"

  if command -v unpigz > /dev/null; then
    unpigz -c "${src}" > "${ungzipped}"
  else
    gzip -dc "${src}" > "${ungzipped}"
  fi
}


cleanup_() {
  if [ -n "${tmp_dir}" ]; then
    rm -rf "${tmp_dir}"
  fi
}


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

output_dir=${args[--output-dir]}
rsem_dir=${output_dir}/rsem
ebseq_dir=${output_dir}/ebseq

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

trap cleanup_ EXIT

if ! ungzip_ "${args[<fasta>]}" "${output_dir}"; then
  echo "Error: failed to decompress ${args[<fasta>]}" >&2
  exit 1
fi
fasta=${ungzipped}

# NOTE: the reference is named after the decompressed file, so a .gz input and
#   its expanded form produce the same reference name.
output_base=$(basename "$(basename "${fasta}" .fasta)" .fa)

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
