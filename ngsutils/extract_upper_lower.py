#!/usr/bin/env python

"""
Extract upper/lower region from feature

Usage:
  extract_upper_lower.py [options] <gtf> <fasta> <feature_ids>...

Options:
  --region <INT>    : Region; The value of plus indicates the upper region
                      from target, and the value of minus indicates the lower
                      region from target, respectively [default: 1000]
  <gtf>             : GTF file
  <fasta>           : Genomic FASTA file
  <feature_ids>...  : Target feature (gene/transcript) ID(s)

"""


import sys

from docopt import docopt
import pandas as pd

from gtfparse import read_gtf
from Bio import SeqIO


def extract_seq(x, seq: dict, region: int):
    if x.strand not in ["+", "-"]:
        sys.stderr.write("Could not determin strand.\n")
        sys.exit(1)

    # Note: Outputs adjusted slice start/end
    def _range(x, region):
        if region > 0:
            if x.strand == "+":
                return (x.start - region - 1, x.start - 1)
            if x.strand == "-":
                return (x.end, x.end + region)
        else:
            if x.strand == "+":
                return (x.end, x.end - region)
            if x.strand == "-":
                return (x.start + region - 1, x.start - 1)

    # def _slice(x, region):
    #     if x.strand == "+":
    #         return slice(*region)
    #     elif x.strand == "-":
    #         return slice(region[0] - 1, region[1] - 1)
    #     else:
    #         sys.stderr.write("Could not determin strand.\n")

    r = _range(x, region)

    seq = seq[x.seqname].seq[slice(*r)]

    if x.strand == "-":
        seq = seq.reverse_complement()

    print(
        "\n".join(
            [
                f">{x.feature_id}{region:+d}|{x.seqname}:{r[0]}-{r[1]}|{x.strand}".format(
                    x.feature_id,
                    region,
                ),
                str(seq),
            ]
        )
    )


def main():
    options = docopt(__doc__)

    fasta_path = options["<fasta>"]
    gtf_path = options["<gtf>"]
    feature_ids = options["<feature_ids>"]
    region = int(options["--region"])

    gtf_df = read_gtf(gtf_path).query("feature == 'exon'")

    features = []

    features.append(
        gtf_df[["gene_id", "seqname", "start", "end", "strand"]]
        .groupby("gene_id")
        .agg({"seqname": min, "start": min, "end": max, "strand": min})
        .rename({"gene_id": "feature_id"}, axis=1)
    )

    features.append(
        gtf_df[["transcript_id", "seqname", "start", "end", "strand"]]
        .groupby("transcript_id")
        .agg({"seqname": min, "start": min, "end": max, "strand": min})
        .rename({"transcript_id": "feature_id"}, axis=1)
    )

    feature_df = pd.concat(features)

    feature_df = (
        feature_df.loc[feature_ids]
        .reset_index()
        .rename({"index": "feature_id"}, axis=1)
    )

    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    _ = feature_df.apply(lambda x: extract_seq(x, genome, region), axis=1)


if __name__ == "__main__":
    main()
