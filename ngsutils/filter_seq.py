#! /bin/env python

"""
Filter transcript FASTA with transcript ids in GTF

Usage:
  filter_seq <fasta> <gtf>

"""

import sys

from Bio import SeqIO
from gtfparse import read_gtf


def main():
    fastq_path = sys.argv[1]
    gtf_path = sys.argv[2]

    gtf_df = read_gtf(gtf_path)
    transcript_ids = set(gtf_df.transcript_id)

    for r in SeqIO.parse(fastq_path, 'fasta'):
        if r.id in transcript_ids:
            print(">{}".format(r.id))
            print(r.seq)


if __name__ == '__main__':
    main()
