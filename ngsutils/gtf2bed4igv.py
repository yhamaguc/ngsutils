#! /usr/bin/env python3

"""
Convert gene annotation GTF to BED

Usage:
  gtf2bed4igv [options] <gtf>

Options:
  --tx-only   : Output transcript records only [default: False]
  --simplify  : Exclude feature name and feature type from ID [default: False]
  <gtf>       : GTF file

"""


import os
from functools import partial

from docopt import docopt
import numpy as np
import pandas as pd
from gtfparse import read_gtf


def pack_name(id: pd.Series, feature_name: pd.Series, biotype: pd.Series):
    name_packed = "ID=" + id + ";Name=" + feature_name + ";BioType=" + biotype
    return name_packed


colors_gene = {"+": "128,0,0", "-": "0,0,128", ".": "0,128,0"}

colors_transcript = {"+": "205,0,0", "-": "0,0,205", ".": "0,205,0"}

ITEM_COLORS = {"gene": colors_gene, "transcript": colors_transcript}


def assign_color(strand, feature):
    return ITEM_COLORS[feature][strand]


def main():
    options = docopt(__doc__)
    gtf_path = options["<gtf>"]
    tx_only = options["--tx-only"]
    simplify = options["--simplify"]

    _pack_name = pack_name

    if simplify:
        _pack_name = lambda id, _y, _z: id

    gtf_df = read_gtf(gtf_path).to_pandas()

    gtf_df.strand = gtf_df.strand.replace("nan", ".")
    # NOTE: gtf_df.columns;
    # Index(['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
    #        'frame', 'havana_gene', 'gene_id', 'tag', 'index', 'gene_type', 'level',
    #        'gene_name', 'transcript_id', 'transcript_name', 'havana_transcript',
    #        'ont', 'transcript_type', 'exon_id', 'exon_number',
    #        'transcript_support_level', 'protein_id', 'ccdsid'],
    #       dtype='object')

    optional_columns = ["transcript_name", "gene_type", "transcript_type"]
    for c in optional_columns:
        if c not in gtf_df.columns:
            gtf_df[c] = "unknown"

    features = ["gene", "transcript", "exon", "CDS", "start_codon", "stop_codon"]

    gtf_dfs = {}
    for f in features:
        gtf_dfs[f] = gtf_df.query("feature == '{}'".format(f))

    bed_dfs = {}
    # NOTE: bed_columns;
    # 'chr', 'start', 'end', 'name', 'score', 'strand',
    # 'thick_start', 'thick_end', 'block_count', 'block_sizes', 'block_starts'

    if not tx_only:
        # Gene record
        feature = "gene"
        assign_color_ = partial(assign_color, feature=feature)
        gtf_df = gtf_dfs[feature]
        zeros_ = np.zeros(len(gtf_df.index)).astype(int)
        starts_ = gtf_df.start - 1
        ends_ = gtf_df.end
        bed_df = pd.DataFrame(
            {
                "chr": gtf_df.seqname,
                "start": starts_,
                "end": ends_,
                "name": _pack_name(gtf_df.gene_id, gtf_df.gene_name, gtf_df.gene_type),
                "score": zeros_,
                "strand": gtf_df.strand,
                "thick_start": starts_,
                "thick_end": ends_,
                "item_rgb": gtf_df.strand.apply(assign_color_),
                "block_count": np.ones(len(gtf_df.index)),
                "block_sizes": (gtf_df.end - gtf_df.start + 1).astype(str) + ",",
                "block_starts": pd.Series(zeros_).astype(str).values + ",",
            }
        )

        bed_dfs["gene"] = bed_df

    # NOTE: Exon record
    gtf_df = gtf_dfs["exon"].copy()

    optional_columns = ["exon_id"]
    for c in optional_columns:
        if c not in gtf_df.columns:
            gtf_df[c] = np.array(range(len(gtf_df)), dtype=np.unicode)

    gtf_df["size"] = gtf_df.end - gtf_df.start + 1
    columns_ = ["transcript_id", "start", "end"]
    outer_df = gtf_df
    inner_df = gtf_dfs["transcript"].filter(columns_)
    gtf_df = pd.merge(outer_df, inner_df, on="transcript_id")
    gtf_df["rel_start"] = gtf_df.start_x - gtf_df.start_y
    gtf_df = gtf_df.sort_values(["transcript_id", "start_x"])
    gtf_df["size"] = gtf_df["size"].astype(str)
    gtf_df["rel_start"] = gtf_df["rel_start"].astype(str)

    sizes_ = gtf_df.groupby("transcript_id")["size"].apply(lambda x: ",".join(x))
    rel_starts_ = gtf_df.groupby("transcript_id")["rel_start"].apply(
        lambda x: ",".join(x)
    )
    counts_ = (
        gtf_df.groupby("transcript_id")["exon_id"]
        .count()
        .rename({"exon_id": "count"})
    )

    # NOTE: CDS record
    gtf_df = gtf_dfs["CDS"].copy()
    gtf_df["start"] = gtf_df.start - 1

    cds_start_min_ = gtf_df[["transcript_id", "start", "end"]].melt(id_vars="transcript_id").groupby("transcript_id")["value"].min()
    cds_end_max_ = gtf_df[["transcript_id", "start", "end"]].melt(id_vars="transcript_id").groupby("transcript_id")["value"].max()

    # NOTE: start/stop codon records; if exists override thick columns
    gtf_df = pd.concat([gtf_dfs["start_codon"].copy(), gtf_dfs["stop_codon"].copy()])
    gtf_df["start"] = gtf_df.start - 1

    ss_codons_start_min_ = gtf_df[["transcript_id", "start", "end"]].melt(id_vars="transcript_id").groupby("transcript_id")["value"].min()
    ss_codons_end_max_ = gtf_df[["transcript_id", "start", "end"]].melt(id_vars="transcript_id").groupby("transcript_id")["value"].max()

    # NOTE: Transcript record
    feature = "transcript"
    assign_color_ = partial(assign_color, feature=feature)
    gtf_df = gtf_dfs[feature].copy()
    zeros_ = np.zeros(len(gtf_df.index)).astype(int)
    starts_ = gtf_df.start - 1
    ends_ = gtf_df.end

    bed_df = pd.DataFrame(
        {
            "chr": gtf_df.seqname,
            "start": starts_,
            "end": ends_,
            "name": _pack_name(
                gtf_df.transcript_id, gtf_df.transcript_name, gtf_df.transcript_type
            ),
            "score": zeros_,
            "strand": gtf_df.strand,
            "thick_start": starts_,
            "thick_end": starts_,
            "item_rgb": gtf_df.strand.apply(assign_color_),
        }
    ).reset_index(drop=True)

    idx = gtf_df.transcript_id

    bed_df = pd.concat(
        [
            bed_df,
            pd.DataFrame(
                {
                    "block_count": counts_[idx].astype(str),
                    "block_sizes": (sizes_[idx].astype(str) + ","),
                    "block_starts": (rel_starts_[idx].astype(str) + ","),
                }
            ).reset_index(drop=True),
        ],
        sort=False,
        axis="columns",
    )

    bed_df.set_index(idx, inplace=True, drop=True)

    bed_df["thick_start"].where(
        cds_start_min_[idx].isna(), cds_start_min_[idx], inplace=True
    )
    bed_df["thick_start"].where(
        ss_codons_start_min_[idx].isna(), ss_codons_start_min_[idx], inplace=True
    )

    bed_df["thick_end"].where(cds_end_max_[idx].isna(), cds_end_max_[idx], inplace=True)
    bed_df["thick_end"].where(
        ss_codons_end_max_[idx].isna(), ss_codons_end_max_[idx], inplace=True
    )

    bed_dfs["transcript"] = bed_df

    bed_df_merged = pd.concat(bed_dfs.values())
    bed_df_merged = bed_df_merged.sort_values(["chr", "start", "name", "end"])

    output_dir = os.path.dirname(gtf_path)
    if output_dir == "":
        output_dir = "."

    gtf_root, _ = os.path.splitext(os.path.basename(gtf_path))
    output_path = os.path.join(output_dir, "{}.bed".format(gtf_root))

    for c in ["score", "thick_start", "thick_end", "block_count"]:
        bed_df_merged[c] = bed_df_merged[c].astype(int)

    with open(output_path, "w") as f:
        f.write("#gffTags\n")
        bed_df_merged.to_csv(f, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
