#!/usr/bin/env bash

#
# Usage:
#   build_star.sh --gtf=gencode.v50.annotation.gtf.gz GRCh38.fa.gz
#   build_star.sh GRCh38.fa.gz
#
#   find ref -name '*.fa.gz' \
#     | xargs -P2 -n1 build_star.sh --gtf=gencode.v50.annotation.gtf.gz --output-dir=index
#
#   # 2-pass: re-build with the junctions found by a 1st-pass alignment
#   find pass1 -name '*SJ.out.tab' | sed 's|^|--sjdb-file=|' \
#     | xargs build_star.sh --gtf=gencode.v50.annotation.gtf.gz --output-dir=index GRCh38.fa.gz
#

DOC="Build a STAR genome index, optionally with a gene annotation and 1st-pass junctions

Usage:
  build_star.sh [--gtf=<PATH>] [--sjdb-file=<PATH>...] [--output-dir=<PATH>] [--threads=<n>] <fasta>
  build_star.sh (-h | --help)

Arguments:
  <fasta>  Genome FASTA, optionally gzipped

Options:
  --gtf=<PATH>         Gene annotation GTF, optionally gzipped; omit to build
                       an index without an annotation
  --sjdb-file=<PATH>   Splice junctions from a 1st-pass alignment (SJ.out.tab),
                       optionally gzipped; repeat the option to insert the
                       junctions of several 1st-pass runs
  --output-dir=<PATH>  Parent of the generated genome directory [default: .]
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

# NOTE: ${out_name} is the one place this copy of ungzip_ differs from the other
#   build_*.sh -- several SJ.out.tab arrive as per-sample files that all share
#   that basename, so the caller has to name the expanded copies apart.
ungzip_() {
  local src=$1
  local work_dir=$2
  local out_name=$3

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

  if [ -z "${out_name}" ]; then
    out_name=$(basename "${src}" .gz)
  fi

  ungzipped=${tmp_dir}/${out_name}
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

# NOTE: a repeatable option lands in the associative array as one entry per
#   occurrence -- args[--sjdb-file,0], [,1], ... -- with the count in
#   args[--sjdb-file,#]. Nothing repeated leaves the count at 0.
n_sjdb_file=${args[--sjdb-file,#]}

# NOTE: this is the choke point for input existence. Everything below assumes
#   the paths resolve; STAR is only reached once they do.
inputs=("${args[<fasta>]}")
if [ -n "${args[--gtf]}" ]; then
  inputs+=("${args[--gtf]}")
fi
for i in $(seq 0 $(( n_sjdb_file - 1 ))); do
  inputs+=("${args[--sjdb-file,${i}]}")
done

n_missing=0
for input in "${inputs[@]}"; do
  if [ ! -f "${input}" ]; then
    echo "Error: no such file: ${input}" >&2
    n_missing=$(( n_missing + 1 ))
  fi
done
if [ "${n_missing}" -gt 0 ]; then
  exit 1
fi

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

trap cleanup_ EXIT

if ! ungzip_ "${args[<fasta>]}" "${output_dir}" ""; then
  echo "Error: failed to decompress ${args[<fasta>]}" >&2
  exit 1
fi
fasta=${ungzipped}

gtf=
if [ -n "${args[--gtf]}" ]; then
  if ! ungzip_ "${args[--gtf]}" "${output_dir}" ""; then
    echo "Error: failed to decompress ${args[--gtf]}" >&2
    exit 1
  fi
  gtf=${ungzipped}
fi

sjdb_files=()
for i in $(seq 0 $(( n_sjdb_file - 1 ))); do
  sjdb_file=${args[--sjdb-file,${i}]}
  if ! ungzip_ "${sjdb_file}" "${output_dir}" "${i}.$(basename "${sjdb_file}" .gz)"; then
    echo "Error: failed to decompress ${sjdb_file}" >&2
    exit 1
  fi
  sjdb_files+=("${ungzipped}")
done

# NOTE: the index is named after the decompressed files, so a .gz input and its
#   expanded form produce the same genome directory.
output_base=$(basename "$(basename "${fasta}" .fasta)" .fa)
if [ -n "${gtf}" ]; then
  output_base=${output_base}.$(basename "${gtf}" .gtf)
fi
# NOTE: the suffix records how many junction files went in, not which ones --
#   two different 1st-pass sample sets of the same size land in the same genome
#   directory and the second build overwrites the first. Give them separate
#   --output-dir when that matters.
if [ "${n_sjdb_file}" -gt 0 ]; then
  output_base=${output_base}.sj${n_sjdb_file}
fi

# NOTE: --genomeDir and mkdir must name the same directory. They did not until
#   2026-08-14: --genomeDir expanded ${output_basel}, a typo for ${output_base},
#   so STAR wrote the index into the output directory itself and left the
#   directory mkdir had just created empty.
genome_dir=${output_dir}/${output_base}

mkdir -p "${genome_dir}"

STAR --version

cmd_=(
  STAR
  --runThreadN "${args[--threads]}"
  --runMode genomeGenerate
  --genomeDir "${genome_dir}"
  --genomeFastaFiles "${fasta}"
)

if [ -n "${gtf}" ]; then
  cmd_+=(--sjdbGTFfile "${gtf}")
fi

# NOTE: STAR reads SJ.out.tab as-is here -- it takes the first four columns
#   (chr, intron start, intron end, strand) and ignores the rest. The junctions
#   are inserted with --sjdbOverhang, left at STAR's default of 100; reads much
#   longer than 101 bp want it set to read length - 1 at build time.
if [ "${n_sjdb_file}" -gt 0 ]; then
  cmd_+=(--sjdbFileChrStartEnd "${sjdb_files[@]}")
fi

echo "CMD: ${cmd_[*]}"
"${cmd_[@]}"
