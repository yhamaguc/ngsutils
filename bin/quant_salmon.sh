#!/usr/bin/env bash

#
# Usage:
#   quant_salmon.sh --index=salmon_index --libtype=ISR --output-dir=out input.bam
#
#   find bam -name '*.bam' \
#     | xargs -P4 -n1 -I{} quant_salmon.sh --index=salmon_index --libtype=ISR \
#         --output-dir=out/{} {}
#

DOC="Quantify transcript abundance with salmon, from BAM or FASTQ

Usage:
  quant_salmon.sh --index=<PATH> --libtype=<TYPE> [--output-dir=<PATH>]
                  [--threads=<n>] <input1> [<input2>]
  quant_salmon.sh (-h | --help)

Arguments:
  <input1>  BAM, or read 1 FASTQ. The file type is taken from the extension:
            a name ending in .bam is converted to FASTQ first
  <input2>  Read 2 FASTQ; omit for single-end, and for BAM input

Options:
  --index=<PATH>       Salmon index directory
  --libtype=<TYPE>     Library type, e.g. ISR, IU, A
  --output-dir=<PATH>  Output directory [default: .]
  --threads=<n>        Threads [default: 8]
  -h --help            Show this message

Examples:
  # BAM input, paired-end library
  quant_salmon.sh --index=index.idx --libtype=ISR --output-dir=out input.bam

  # FASTQ paired-end
  quant_salmon.sh --index=index.idx --libtype=ISR --output-dir=out r1.fq.gz r2.fq.gz

  # FASTQ single-end
  quant_salmon.sh --index=index.idx --libtype=U --output-dir=out reads.fq.gz
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

input1=${args[<input1>]}
libtype=${args[--libtype]}
output_dir=${args[--output-dir]}

mkdir -p "${output_dir}"

salmon_=(
  salmon quant
  -p "${args[--threads]}"
  -i "${args[--index]}"
  -l "${libtype}"
)

r1=
r2=
sr=

if [[ "${input1}" =~ \.bam$ ]]; then
  filebase=$(basename "${input1}" .bam)
  tmp_dir=${output_dir}/tmp_fastq
  mkdir -p "${tmp_dir}"

  bamtofastq_=(
    bamtofastq
    gz=1
    fastq=0
    disablevalidation=1
    "S=${tmp_dir}/${filebase}.sr.fastq.gz"
    "F=${tmp_dir}/${filebase}.r1.fastq.gz"
    "F2=${tmp_dir}/${filebase}.r2.fastq.gz"
    O=/dev/null
    O2=/dev/null
    "filename=${input1}"
  )

  echo "CMD: ${bamtofastq_[*]}"

  # NOTE: salmon would happily quantify a truncated FASTQ and report a plausible
  #   TPM table, so the conversion is checked before it is consumed.
  if ! "${bamtofastq_[@]}"; then
    echo "Error: bamtofastq failed on ${input1}; leaving ${tmp_dir} in place" >&2
    exit 1
  fi

  # NOTE: an 'I' prefix on the library type is salmon's marker for an inward
  #   paired-end library, so it is what decides which of bamtofastq's outputs
  #   to hand back to salmon.
  if [[ "${libtype:0:1}" == "I" ]]; then
    r1=${tmp_dir}/${filebase}.r1.fastq.gz
    r2=${tmp_dir}/${filebase}.r2.fastq.gz
  else
    sr=${tmp_dir}/${filebase}.sr.fastq.gz
  fi
else
  r1=${input1}
  r2=${args[<input2>]}
fi

if [ -n "${r2}" ]; then
  salmon_+=(-1 "${r1}" -2 "${r2}")
elif [ -n "${sr}" ]; then
  salmon_+=(-r "${sr}")
else
  salmon_+=(-1 "${r1}")
fi

salmon_+=(-o "${output_dir}")

echo "CMD: ${salmon_[*]}"
"${salmon_[@]}"
salmon_status=$?

# NOTE: the temporary FASTQs are only removed on success. A failed run keeps
#   them so the salmon invocation can be repeated against the same input.
if [ -n "${tmp_dir:-}" ]; then
  if [ ${salmon_status} -eq 0 ]; then
    rm -rf "${tmp_dir}"
  else
    echo "Error: salmon exited ${salmon_status}; leaving ${tmp_dir} in place" >&2
  fi
fi

exit ${salmon_status}
