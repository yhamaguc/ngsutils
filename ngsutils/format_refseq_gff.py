"""
Format RefSeq GFF to GTF and Parsed SQLite

Usage:
   format_refseq_gff.py <gff> <feature_table> <assembly_conversion_table> <fasta>

"""


import sys
import os
import re
import csv

import numpy as np
import pandas as pd
import sqlite3
from Bio import SeqIO


REQUIRED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]


def output_sqlite(df, path):
    path = "{}.sqlite".format(path)
    conn = sqlite3.connect(path)
    target_table = 'annotations'
    df.to_sql(target_table, conn, if_exists='replace')


def output_gtf(df, path):
    path = "{}.gtf".format(path)
    df.to_csv(path, index=False, header=False, sep='\t', quoting=csv.QUOTE_NONE)


def output_fasta(path, keep_ids, dup_ids, postfix='_chrY'):
    root, ext = os.path.splitext(path)
    output_path = "{}_filtered{}".format(root, ext)

    with open(output_path, 'w') as f:
        for r in SeqIO.parse(path, 'fasta'):
            if r.id in keep_ids:
                f.write(">{}\n".format(r.id))
                f.write("{}\n".format(r.seq))
                if r.id in dup_ids:
                    f.write(">{}{}\n".format(r.id, postfix))
                    f.write("{}\n".format(r.seq))


def load_gff(path):
    # path = gff_path
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=REQUIRED_COLUMNS)

    df = df.query("feature == 'exon'")
    attributes = df.attribute.apply(lambda x: re.split('[=;]', x))

    attribute_cols = list()
    attributes.apply(lambda x: attribute_cols.extend(x[0::2]))

    def extract_subcolumn(x, col):
        try:
            val = x[x.index(col) + 1]
        except Exception:
            val = ''
        return val

    REQUIRED_COLUMNS_ATTRIBUTE = ['ID', 'gene', 'gbkey', 'transcript_id', 'Dbxref']
    REQUIRED_COLUMNS_ATTRIBUTE_RENAMED = ['exon_id', 'gene_name', 'gbkey', 'transcript_id', 'dbxref']
    attribute_serieses = {}
    for col, key in zip(REQUIRED_COLUMNS_ATTRIBUTE, REQUIRED_COLUMNS_ATTRIBUTE_RENAMED):
        attribute_serieses[key] = attributes.apply(lambda x: extract_subcolumn(x, col))

    attributes_dbxref = attribute_serieses['dbxref'].apply(lambda x: re.split('[,:]', x))

    REQUIRED_COLUMNS_DBXREF = ['GeneID']
    REQUIRED_COLUMNS_DBXREF_RENAMED = ['gene_id']
    for col, key in zip(REQUIRED_COLUMNS_DBXREF, REQUIRED_COLUMNS_DBXREF_RENAMED):
        attribute_serieses[key] = attributes_dbxref.apply(lambda x: extract_subcolumn(x, col))

    df = df.drop('attribute', axis=1)
    ORDERED_ATTRIBUTE_KEYS = ['gene_id', 'gene_name', 'transcript_id', 'exon_id', 'gbkey', 'dbxref']
    for k in ORDERED_ATTRIBUTE_KEYS:
        df[k] = attribute_serieses[k]

    output_sqlite(df, path + '.tmp')
    df_omitted = df.query("gene_id == '' | transcript_id == ''")
    df = df.query("gene_id != '' & transcript_id != ''")
    output_sqlite(df_omitted, path + '.omit')

    return df


def load_feature(path):
    df = pd.read_csv(
        path,
        sep="\t")

    cond = df['class'].isnull()
    df.loc[cond, 'class'] = df.loc[cond, 'feature']
    df = df.query("product_accession == product_accession").filter(['class', 'product_accession', 'name']).drop_duplicates(keep='first')

    return df


def load_assembly(path):
    df = pd.read_csv(
        path,
        sep="\t")

    df = df[df['Sequence-Role'] == 'assembled-molecule'].filter(['RefSeq-Accn', 'UCSC-style-name'])

    return df


def attributes_to_str(df, attribute_keys, replacement=''):
    def formatted(df, k):
        return df[k].fillna(replacement).apply(lambda x: '{} "{}";'.format(k, x))

    attributes = None

    for k in attribute_keys:
        if attributes is None:
            attributes = formatted(df, k)
            continue

        attributes = attributes + formatted(df, k)

    return attributes


if __name__ == '__main__':
    if __debug__:
        gff_path = sys.argv[1]
        feature_path = sys.argv[2]
        assembly_path = sys.argv[3]
        fasta_path = sys.argv[4]
    else:
        print('Run as debug mode')
        base_dir = os.path.join(os.path.expanduser('~'), '/home/yamada/Shared/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13/')
        gff_path = os.path.join(base_dir, 'GCF_000001405.39_GRCh38.p13_genomic.gff')
        feature_path = os.path.join(base_dir, 'GCF_000001405.39_GRCh38.p13_feature_table.txt')
        assembly_path = os.path.join(base_dir, 'GCF_000001405.39_GRCh38.p13_assembly_table.txt')
        fasta_path = os.path.join(base_dir, 'GCF_000001405.39_GRCh38.p13_rna.fna')

    gff_df = load_gff(gff_path)
    feature_df = load_feature(feature_path)
    assembly_df = load_assembly(assembly_path)

    gff_df_merged = gff_df.merge(assembly_df, how='inner', left_on='seqname', right_on='RefSeq-Accn').merge(feature_df, how='left', left_on='transcript_id', right_on='product_accession')
    set(gff_df_merged['seqname'])

    gff_df_merged = gff_df_merged.drop(['seqname', 'product_accession', 'RefSeq-Accn'], axis=1)
    gff_df_merged = gff_df_merged.rename(columns={'UCSC-style-name': 'seqname', 'name': 'description'})
    gff_df_merged['transcript_name'] = gff_df_merged['gene_name'] + '-' + gff_df_merged['transcript_id']
    ORDERED_ATTRIBUTE_KEYS = ['gene_id', 'gene_name', 'transcript_id', 'transcript_name', 'exon_id', 'class', 'gbkey', 'description', 'dbxref']

    def count_distinct(x):
        return len(set(list(x)))

    def group_concat(x):
        return ','.join(set(list(x)))

    funcs_agg = {
            'seqname': [count_distinct, group_concat],
            }

    dfg = gff_df_merged.groupby('transcript_id').agg(funcs_agg)
    dfg.columns = ['_'.join(c) for c in dfg.columns.tolist()]

    dfg_dup = dfg.query("seqname_count_distinct > 1")
    if np.any('chrX,chrY' != dfg_dup['seqname_group_concat']) and np.any('chrY,chrX' != dfg_dup['seqname_group_concat']):
        raise Exception("Duplication records are NOT in PAR. {}".format('|'.join(set(dfg_dup['seqname_group_concat']))))

    dup_transcript_ids = dfg_dup.index
    cond = (gff_df_merged['transcript_id'].isin(dup_transcript_ids)) & (gff_df_merged['seqname'] == 'chrY')

    gff_df_merged.loc[cond, 'transcript_id'] = gff_df_merged.loc[cond, 'transcript_id'] + '_' + gff_df_merged.loc[cond, 'seqname']
    gff_df_merged.loc[cond, 'gene_id'] = gff_df_merged.loc[cond, 'gene_id'] + '_' + gff_df_merged.loc[cond, 'seqname']

    output_sqlite(gff_df_merged.filter(REQUIRED_COLUMNS + ORDERED_ATTRIBUTE_KEYS), gff_path)

    attributes = attributes_to_str(gff_df_merged, ORDERED_ATTRIBUTE_KEYS)
    gtf_df = gff_df_merged.filter(REQUIRED_COLUMNS)
    gtf_df['attribute'] = attributes
    output_gtf(gtf_df, gff_path)

    keep_transcript_ids = set(gff_df_merged['transcript_id'])
    output_fasta(fasta_path, keep_transcript_ids, dup_transcript_ids)
