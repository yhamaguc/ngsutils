#! /usr/bin/env python3
import sys
import os
from functools import partial

import numpy as np
import pandas as pd
from gtfparse import read_gtf

"""
Convert positions from genome corrdinates to transcript coordinates

Usage:
  abspos2relpos <gtf>

"""


def pack_name(id: pd.Series, feature_name: pd.Series, biotype: pd.Series):
    name_packed = "ID=" + id + ";Name=" + feature_name + ";BioType=" + biotype
    return name_packed


colors_gene = {'+': '128,0,0',
               '-': '0,0,128',
               '.': '0,128,0'}

colors_transcript = {'+': '205,0,0',
                     '-': '0,0,205',
                     '.': '0,205,0'}

ITEM_COLORS = {'gene': colors_gene, 'transcript': colors_transcript}


def assign_color(strand, feature):
    return ITEM_COLORS[feature][strand]


def main():
    gtf_path = sys.argv[1]
    gtf_path = '/Volumes/External/Assets/annotations/human/GRCh38/gencode/gencode.v29.annotation.gtf'

    gtf_df = read_gtf(gtf_path)
    # gtf_df.columns
    # Index(['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
    #        'frame', 'havana_gene', 'gene_id', 'tag', 'index', 'gene_type', 'level',
    #        'gene_name', 'transcript_id', 'transcript_name', 'havana_transcript',
    #        'ont', 'transcript_type', 'exon_id', 'exon_number',
    #        'transcript_support_level', 'protein_id', 'ccdsid'],
    #       dtype='object')

    features = ['gene', 'transcript', 'exon']

    gtf_dfs = {}
    for f in features:
        gtf_dfs[f] = gtf_df.query("feature == '{}'".format(f))

    bed_dfs = {}
    # NOTE: bed_columns
    # 'chr', 'start', 'end', 'name', 'score', 'strand',
    # 'thick_start', 'thick_end', 'block_count', 'block_sizes', 'block_starts'

    # Gene record
    feature = 'gene'
    assign_color_ = partial(assign_color, feature=feature)
    gtf_df = gtf_dfs[feature]
    zeros_ = np.zeros(len(gtf_df.index)).astype(int)
    starts_ = gtf_df.start - 1
    ends_ = gtf_df.end
    bed_df = pd.DataFrame({
        'chr': gtf_df.gene_id,
        'start': zeros_,
        'end': ends_ - starts_ + 1,
        'name': pack_name(gtf_df.gene_id, gtf_df.gene_name, gtf_df.gene_type),
        'score': zeros_,
        'strand': gtf_df.strand,
        'thick_start': ends_,
        'thick_end': ends_,
        'item_rgb': gtf_df.strand.apply(assign_color_),
        'block_count': np.ones(len(gtf_df.index)),
        'block_sizes': (gtf_df.end - gtf_df.start + 1).astype(str) + ',',
        'block_starts': pd.Series(zeros_).astype(str).values + ','})

    bed_dfs['gene'] = bed_df

    # Exon record
    gtf_df = gtf_dfs['exon'].copy()
    gtf_df['size'] = gtf_df.end - gtf_df.start + 1
    gtf_df.columns
    columns_ = ['gene_id', 'start', 'end']
    outer_df = gtf_df
    inner_df = gtf_dfs['gene'].filter(columns_)
    gtf_df = pd.merge(outer_df, inner_df, on='gene_id', left_index=True)
    gtf_df['rel_start'] = gtf_df.start_x - gtf_df.start_y
    gtf_df = gtf_df.sort_values(['gene_id', 'start_x'])
    gtf_df['size'] = gtf_df['size'].astype(str)
    gtf_df['rel_start'] = gtf_df['rel_start'].astype(str)

    sizes_ = gtf_df.groupby('gene_id')['size'].apply(lambda x: ','.join(x))
    rel_starts_ = gtf_df.groupby('gene_id')['rel_start'].apply(lambda x: ','.join(x))
    counts_ = gtf_df.groupby('gene_id')['exon_id'].count().rename(columns={'exon_id': 'count'})

    feature = 'gene'
    gtf_df = gtf_dfs[feature]
    bed_df['block_count'] = counts_[gtf_df.gene_id].values
    bed_df['block_sizes'] = sizes_[gtf_df.gene_id].values + ','
    bed_df['block_starts'] = rel_starts_[gtf_df.gene_id].values + ','

#    # Transcript record
#    feature = 'transcript'
#    assign_color_ = partial(assign_color, feature=feature)
#    gtf_df = gtf_dfs[feature].copy()
#    zeros_ = np.zeros(len(gtf_df.index)).astype(int)
#    starts_ = gtf_df.start - 1
#    ends_ = gtf_df.end
#    bed_df = pd.DataFrame({
#        'chr': gtf_df.seqname,
#        'start': starts_,
#        'end': ends_,
#        'name': pack_name(gtf_df.transcript_id,
#                          gtf_df.transcript_name,
#                          gtf_df.transcript_type),
#        'score': zeros_,
#        'strand': gtf_df.strand,
#        'thick_start': ends_,
#        'thick_end': ends_,
#        'item_rgb': gtf_df.strand.apply(assign_color_)})
#
#    bed_df['block_count'] = counts_[gtf_df.transcript_id].values
#    bed_df['block_sizes'] = sizes_[gtf_df.transcript_id].values + ','
#    bed_df['block_starts'] = rel_starts_[gtf_df.transcript_id].values + ','
#
#    bed_dfs['transcript'] = bed_df

#    bed_df_merged = pd.concat(bed_dfs.values())
#    bed_df_merged = bed_df_merged.sort_values(['chr', 'start', 'name', 'end'])
    bed_df_merged = bed_df

    output_dir = os.path.dirname(gtf_path)
    if output_dir == '':
        output_dir = '.'

    gtf_root, _ = os.path.splitext(os.path.basename(gtf_path))
    output_path = os.path.join(output_dir, "{}_2.bed".format(gtf_root))
    
    for c in ['score', 'thick_start', 'thick_end', 'block_count']:
        bed_df_merged[c] = bed_df_merged[c].astype(int)


    bed_df_merged['chr'] = bed_df_merged.chr.str.replace('^ENSG', 'IMMT')

    with open(output_path, 'w') as f:
        f.write("#gffTags\n")
        bed_df_merged.to_csv(f,  sep="\t", header=False, index=False)


if __name__ == '__main__':
    main()
