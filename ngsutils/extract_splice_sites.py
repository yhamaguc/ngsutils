#!/usr/bin/env python

"""
Extract splice site positions from GTF

Usage:
  extract_splice_sites.py [options] <gtf>

Options:
  --width <INT>  : donor / acceptor width [default: 2]
  <gtf>          : GTF file

"""

import sys

from docopt import docopt
import pandas as pd
import numpy as np

from gtfparse import read_gtf


def sort_(X):
    return X.sort_values(by=['seqname', 'start', 'end', 'strand'])


def to0base(X):
    X.start = X.start - 1
    X.end = X.end
    return X


def extract_splice_sites(exons, strand, class_, width=2):
    if class_ == IS_SINGLE:
        return None

    exons['name'] = exons['transcript_id'] + \
        '.' + exons['exon_number_x'].astype(str)

    def donor(exon, strand):
        if exon['strand'] == '+':
            return (exon['seqname'], exon['end'] + 1, exon['end'] + width, exon['name'] + 'd', 0, strand)
        elif exon['strand'] == '-':
            return (exon['seqname'], exon['start'] - width, exon['start'] - 1, exon['name'] + 'd', 0, strand)

    def acceptor(exon, strand):
        if exon['strand'] == '+':
            return (exon['seqname'], exon['start'] - width, exon['start'] - 1, exon['name'] + 'a', 0, strand)
        elif exon['strand'] == '-':
            return (exon['seqname'], exon['end'] + 1, exon['end'] + width, exon['name'] + 'a', 0, strand)

    if class_ == IS_FIRST:
        return exons.apply(lambda r: donor(r, strand), axis=1)

    if class_ == IS_LAST:
        return exons.apply(lambda r: acceptor(r, strand), axis=1)

    if class_ == IS_MIDDLE:
        return (
            pd.concat([exons.apply(lambda r: donor(r, strand), axis=1),
                       exons.apply(lambda r: acceptor(r, strand), axis=1)],
                      axis=0,
                      ignore_index=True
                      )
        )


def collapse(X):
    X = X.groupby(['seqname', 'start', 'end', 'score', 'strand']).agg({'name' : lambda x: ';'.join(x.sort_values())})
    X = X.drop_duplicates().reset_index().reindex(columns=['seqname', 'start', 'end', 'name', 'score', 'strand'])

    return X


IS_SINGLE = 0
IS_FIRST = 1
IS_MIDDLE = 2
IS_LAST = 3


def encode(x, y):
    d = {
        (True, True): IS_SINGLE,
        (True, False): IS_FIRST,
        (False, True): IS_LAST,
        (False, False): IS_MIDDLE
    }

    return d[(x, y)]


def main():
    options = docopt(__doc__)

    gtf_path = options["<gtf>"]
    # gtf_path = "/home/yhamaguc/share/assets/references/annotations/gencode.v45.annotation.gtf"

    exons = read_gtf(gtf_path).to_pandas().query("feature == 'exon'")
    exons.exon_number = exons.exon_number.astype(int)

    exons = exons.merge(
        exons.groupby('transcript_id')['exon_number'].max(),
        on='transcript_id',
        how='left'
    ).reset_index().pipe(sort_)

    exons['is_first'] = np.where((exons['exon_number_x'] == 1), True, False)
    exons['is_last'] = np.where(
        (exons['exon_number_x'] == exons['exon_number_y']), True, False)

    exons['class'] = exons.apply(lambda r: encode(
        r['is_first'], r['is_last']), axis=1)

    byclass_splice_sites = {k: extract_splice_sites(v, *k)
                            for k, v in exons.groupby(['strand', 'class'])}

    bed = pd.concat(
        [pd.DataFrame(v.tolist(), columns=['seqname', 'start', 'end', 'name', 'score', 'strand'])
         for v in byclass_splice_sites.values() if v is not None],
        axis=0,
        ignore_index=True
    ).pipe(sort_).pipe(to0base).pipe(collapse)

    bed.to_csv(sys.stdout, sep='\t', index=False, header=False)


if __name__ == "__main__":
    main()
