#! /usr/bin/env python3
"""
Convert positions from genome coordinates to transcript coordinates

Usage:
  abs2rel <gtf>

Arguments:
  <gtf>  GTF file, .gz accepted

"""

import sys
import os
from functools import partial

from docopt import docopt
import numpy as np
import pandas as pd
from ngsutils.gtf import gtf_stem, read_gtf


def pack_name(id: pd.Series, feature_name: pd.Series, biotype: pd.Series):
    name_packed = "ID=" + id + ";Name=" + feature_name + ";BioType=" + biotype
    return name_packed


colors_gene = {
    '+': '128,0,0',
    '-': '0,0,128',
    '.': '0,128,0'
}

colors_transcript = {
    '+': '205,0,0',
    '-': '0,0,205',
    '.': '0,205,0'
}


ITEM_COLORS = {'gene': colors_gene, 'transcript': colors_transcript}


def assign_color(strand, feature):
    return ITEM_COLORS[feature][strand]


def main():
    options = docopt(__doc__)
    gtf_path = options['<gtf>']

    gtf_df = read_gtf(gtf_path, result_type='pandas')

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
    # NOTE: the gene becomes its own reference sequence here, so every coordinate below is
    #   an offset from the gene start and the record spans [0, length) in the 0-based
    #   half-open BED frame. length is the 1-based inclusive GTF span, end - start + 1,
    #   which is ends_ - starts_ because starts_ has already taken the base off. Adding
    #   one again -- as this line did until 2026-09-09 -- counts one base twice and left
    #   the final block ending one short of the record end, which is not valid BED12.
    lengths_ = ends_ - starts_
    bed_df = pd.DataFrame({
        'chr': gtf_df.gene_id,
        'start': zeros_,
        'end': lengths_,
        'name': pack_name(gtf_df.gene_id, gtf_df.gene_name, gtf_df.gene_type),
        'score': zeros_,
        'strand': gtf_df.strand,
        # NOTE: thick_start == thick_end == 0 is BED's "no coding range marked". abs2rel
        #   reads gene, transcript and exon only -- no CDS and no codon features -- so it
        #   has nothing to mark. Both columns held ends_, an absolute GTF coordinate,
        #   which fell outside the record it belonged to.
        'thick_start': zeros_,
        'thick_end': zeros_,
        'item_rgb': gtf_df.strand.apply(assign_color_),
        'block_count': np.ones(len(gtf_df.index)),
        'block_sizes': lengths_.astype(str) + ',',
        'block_starts': pd.Series(zeros_).astype(str).values + ','})

    bed_dfs['gene'] = bed_df

    # Exon record
    gtf_df = gtf_dfs['exon'].copy()
    gtf_df['size'] = gtf_df.end - gtf_df.start + 1
    gtf_df.columns
    columns_ = ['gene_id', 'start', 'end']
    outer_df = gtf_df
    inner_df = gtf_dfs['gene'].filter(columns_)
    # NOTE: on='gene_id' alone, as gtf2bed's equivalent merge does. The
    #   left_index=True this line carried until 2026-09-09 is rejected by pandas
    #   in combination with on=, so the line had most likely never run.
    gtf_df = pd.merge(outer_df, inner_df, on='gene_id')
    gtf_df['rel_start'] = gtf_df.start_x - gtf_df.start_y
    gtf_df = gtf_df.sort_values(['gene_id', 'start_x'])
    gtf_df['size'] = gtf_df['size'].astype(str)
    gtf_df['rel_start'] = gtf_df['rel_start'].astype(str)

    sizes_ = gtf_df.groupby('gene_id')['size'].apply(lambda x: ','.join(x))
    rel_starts_ = gtf_df.groupby(
        'gene_id')['rel_start'].apply(lambda x: ','.join(x))
    # NOTE: no rename. The .rename(columns={'exon_id': 'count'}) this line carried is a
    #   TypeError on a Series, and gtf2bed's .rename({'exon_id': 'count'}) form renames
    #   index LABELS, none of which is 'exon_id' -- so it renamed nothing there either.
    #   Nothing below reads this by name; block_count indexes it by gene_id.
    counts_ = gtf_df.groupby('gene_id')['exon_id'].count()

    feature = 'gene'
    gtf_df = gtf_dfs[feature]
    bed_df['block_count'] = counts_[gtf_df.gene_id].values
    bed_df['block_sizes'] = sizes_[gtf_df.gene_id].values + ','
    bed_df['block_starts'] = rel_starts_[gtf_df.gene_id].values + ','

    bed_df_merged = bed_df

    output_dir = os.path.dirname(gtf_path)
    if output_dir == '':
        output_dir = '.'

    output_path = os.path.join(
        output_dir, "{}_2.bed".format(gtf_stem(gtf_path)))

    for c in ['score', 'thick_start', 'thick_end', 'block_count']:
        bed_df_merged[c] = bed_df_merged[c].astype(int)

    with open(output_path, 'w') as f:
        f.write("#gffTags\n")
        bed_df_merged.to_csv(f,  sep="\t", header=False, index=False)


if __name__ == '__main__':
    main()
