#! /usr/bin/env python3

"""
Extract gene-transcript relationship from GTF

Usage:
  gene2tx [options] <gtf>

Options:
  -o --output-dir <PATH>  : Output directory [default: .]
  <gtf>                : GTF formatted gene annotation file

"""

import os

from docopt import docopt
from gtfparse import read_gtf


def main():
    options = docopt(__doc__)
    gtf_path = options['<gtf>']

    root, _ = os.path.splitext(os.path.basename(gtf_path))

    output_path = os.path.join(
        options['--output-dir'],
        f"{root}.gene2tx.txt"
    )

    cols = ['gene_id', 'transcript_id']
    annotations = read_gtf(gtf_path)
    annotations = annotations.query(
        "feature == 'exon'"
    ).filter(cols).sort_values(by=cols)

    annotations = annotations.drop_duplicates(cols)
    annotations.to_csv(output_path, header=False, index=False, sep='\t')


if __name__ == '__main__':
    main()
