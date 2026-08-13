#!/usr/bin/env bash

#
# Usage:
#   download_gencode_hs.sh --output-dir=ref/gencode_v50 50
#
# NOTE: not written for the find | xargs form the other scripts use -- one
#   release is one invocation, and the GENCODE mirror is the bottleneck.
#

DOC="Download a GENCODE human release, then index the annotation GTFs

Usage:
  download_gencode_hs.sh [--output-dir=<PATH>] [--threads=<n>] <release>
  download_gencode_hs.sh (-h | --help)

Arguments:
  <release>  GENCODE release version, e.g. 50

Options:
  --output-dir=<PATH>  Directory the files are downloaded into [default: .]
  --threads=<n>        Threads passed to bgzip [default: 8]
  -h --help            Show this message
"

#
# Constants
#
BASE_URL_TEMPLATE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_%s"

SUFFIXES=(
  2wayconspseudos.gff3.gz
  2wayconspseudos.gtf.gz
  annotation.gff3.gz
  annotation.gtf.gz
  basic.annotation.gff3.gz
  basic.annotation.gtf.gz
  chr_patch_hapl_scaff.annotation.gff3.gz
  chr_patch_hapl_scaff.annotation.gtf.gz
  chr_patch_hapl_scaff.basic.annotation.gff3.gz
  chr_patch_hapl_scaff.basic.annotation.gtf.gz
  lncRNA_transcripts.fa.gz
  long_noncoding_RNAs.gff3.gz
  long_noncoding_RNAs.gtf.gz
  metadata.Annotation_remark.gz
  metadata.EntrezGene.gz
  metadata.Exon_supporting_feature.gz
  metadata.Gene_source.gz
  metadata.HGNC.gz
  metadata.PDB.gz
  metadata.PolyA_feature.gz
  metadata.Pubmed_id.gz
  metadata.RefSeq.gz
  metadata.Selenocysteine.gz
  metadata.SwissProt.gz
  metadata.TrEMBL.gz
  metadata.Transcript_source.gz
  metadata.Transcript_supporting_feature.gz
  pc_transcripts.fa.gz
  pc_translations.fa.gz
  polyAs.gff3.gz
  polyAs.gtf.gz
  primary_assembly.annotation.gff3.gz
  primary_assembly.annotation.gtf.gz
  primary_assembly.basic.annotation.gff3.gz
  primary_assembly.basic.annotation.gtf.gz
  promoter_windows.gff3.gz
  tRNAs.gff3.gz
  tRNAs.gtf.gz
  transcript_rankings.txt.gz
  transcripts.fa.gz
)

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

release=${args[<release>]}
output_dir=${args[--output-dir]}

if [ ! -e "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

# NOTE: everything below names files relative to the download directory --
#   unpigz, gtf2bed.sh, gtf2sqlite and tabix all write beside their input -- so
#   the directory is entered once here rather than threaded through each call.
if ! cd "${output_dir}"; then
  echo "Error: cannot enter ${output_dir}" >&2
  exit 1
fi

# shellcheck disable=SC2059
base_url=$(printf "${BASE_URL_TEMPLATE}" "${release}")

echo "Downloading GENCODE v${release} files into $(pwd) ..."

# NOTE: not every suffix exists in every release, so a failed download is
#   tolerated per file. The count is reported below and again before indexing --
#   a partial release must not read as a complete one.
n_failed=0
for suffix in "${SUFFIXES[@]}"; do
  if ! wget "${base_url}/gencode.v${release}.${suffix}"; then
    echo "Warning: failed to download ${suffix}" >&2
    n_failed=$((n_failed + 1))
  fi
done

if [ ${n_failed} -gt 0 ]; then
  echo "Warning: ${n_failed} of ${#SUFFIXES[@]} files failed to download" >&2
fi

echo "Decompressing files..."
unpigz gencode.v"${release}".*.gz

indexing_targets=(
  gencode.v"${release}".basic.annotation.gtf
  gencode.v"${release}".annotation.gtf
)

n_indexed=0
n_skipped=0
for target in "${indexing_targets[@]}"; do
  if [ ! -f "${target}" ]; then
    echo "Skip: ${target} not found." >&2
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo "Converting ${target} to BED and SQLITE format and indexing..."

  if ! gtf2bed.sh "${target}"; then
    echo "Warning: gtf2bed.sh failed on ${target}; not indexing it" >&2
    n_skipped=$((n_skipped + 1))
    continue
  fi

  if ! gtf2sqlite --gtf "${target}"; then
    echo "Warning: gtf2sqlite failed on ${target}" >&2
  fi

  # NOTE: tabix needs the bgzip to have finished, so an unindexed .bed.gz is
  #   never left looking queryable.
  if ! bgzip --threads "${args[--threads]}" "${target%.gtf}.bed"; then
    echo "Warning: bgzip failed on ${target%.gtf}.bed; not indexing it" >&2
    n_skipped=$((n_skipped + 1))
    continue
  fi

  if ! tabix -p bed "${target%.gtf}.bed.gz"; then
    echo "Warning: tabix failed on ${target%.gtf}.bed.gz" >&2
    n_skipped=$((n_skipped + 1))
    continue
  fi

  n_indexed=$((n_indexed + 1))
done

echo "Downloaded ${#SUFFIXES[@]} targets with ${n_failed} failure(s);" \
     "indexed ${n_indexed} of ${#indexing_targets[@]} annotation GTF(s)," \
     "${n_skipped} skipped."

if [ ${n_failed} -gt 0 ] || [ ${n_skipped} -gt 0 ]; then
  exit 1
fi
