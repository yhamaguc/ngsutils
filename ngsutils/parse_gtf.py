#! /usr/bin/env python3

"""
Parse GTF to pandas dataframe

Usage:
  parse_gtf <gtf>

"""

import sys
import csv

from docopt import docopt
import pandas as pd


COLUMN_NAMES = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute',
]


def expand_attributes(df, quote_chr='\"', replacement=''):
    def to_key_value(attribute):
        if not attribute:
            return {}

        return dict(kv.strip().split(' ') for kv in attribute.split(';') if kv)

    attributes = df['attribute'].replace('; ', ';', regex=True)
    expanded = pd.DataFrame(
        [to_key_value(a) for a in attributes],
        ).fillna(replacement).applymap(lambda x: x.strip(quote_chr))

    # FIXME: Save order
    df_expanded = pd.concat([
        df.iloc[:, 0:8],
        expanded
    ],
        axis=1,
        sort=False
    )

    del expanded

    return df_expanded


def read_gtf(path, names=COLUMN_NAMES):
    return expand_attributes(
        pd.read_csv(
            path,
            sep="\t",
            comment='#',
            names=names,
            skipinitialspace=True,
            skip_blank_lines=True,
            error_bad_lines=False,
            warn_bad_lines=True,
            dtype=str
        )
    )


def main():
    opt = docopt(__doc__)
    gtf_path = opt['<gtf>']

    gtf = read_gtf(gtf_path)

    gtf.to_csv(
        sys.stdout,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE
        )


if __name__ == '__main__':
    main()
