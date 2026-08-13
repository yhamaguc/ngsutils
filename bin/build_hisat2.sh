#!/usr/bin/env bash

#
# Usage:
#   build_hisat2.sh --gtf=gencode.v50.annotation.gtf.gz GRCh38.fa.gz
#
#   find ref -name '*.fa.gz' \
#     | xargs -P2 -n1 build_hisat2.sh --gtf=gencode.v50.annotation.gtf.gz --output-dir=index
#

DOC="Build a HISAT2 index with splice sites and exons from a GTF

Usage:
  build_hisat2.sh --gtf=<PATH> [--output-dir=<PATH>] [--threads=<n>] <fasta>
  build_hisat2.sh (-h | --help)

Arguments:
  <fasta>  Genome FASTA, optionally gzipped

Options:
  --gtf=<PATH>         Gene annotation GTF, optionally gzipped
  --output-dir=<PATH>  Output directory [default: .]
  --threads=<n>        Threads [default: 8]
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

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

trap cleanup_ EXIT

if ! ungzip_ "${args[<fasta>]}" "${output_dir}"; then
  echo "Error: failed to decompress ${args[<fasta>]}" >&2
  exit 1
fi
fasta=${ungzipped}

if ! ungzip_ "${args[--gtf]}" "${output_dir}"; then
  echo "Error: failed to decompress ${args[--gtf]}" >&2
  exit 1
fi
gtf=${ungzipped}

# NOTE: the index is named after the decompressed files, so a .gz input and its
#   expanded form produce the same index prefix.
output_base=$(basename "$(basename "${fasta}" .fasta)" .fa).$(basename "${gtf}" .gtf)
output_splicesites=${output_dir}/${output_base}.splicesites.txt
output_exons=${output_dir}/${output_base}.exons.txt

# NOTE: hisat2-build below reads both of these files. An empty splice-site or
#   exon table is silently accepted there and yields an index with no annotation
#   in it, which is indistinguishable from build_hisat2_wo_gtf.sh output. The
#   extraction status is therefore checked here rather than left to the caller.
if ! hisat2_extract_splice_sites.py "${gtf}" > "${output_splicesites}"; then
  echo "Error: hisat2_extract_splice_sites.py failed on ${gtf}" >&2
  exit 1
fi

if ! hisat2_extract_exons.py "${gtf}" > "${output_exons}"; then
  echo "Error: hisat2_extract_exons.py failed on ${gtf}" >&2
  exit 1
fi

cmd_=(
  hisat2-build
  -p "${args[--threads]}"
  --ss "${output_splicesites}"
  --exon "${output_exons}"
  "${fasta}"
  "${output_dir}/${output_base}"
)

echo "CMD: ${cmd_[*]}"
"${cmd_[@]}"
