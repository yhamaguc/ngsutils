#! /usr/bin/env python3

"""
Extract gene-transcript relationship from GTF

Usage:
  gene2tx [options] <gtf>

Options:
  -o --output-dir <PATH>  : Output directory [default: .]
  <gtf>                : GTF formatted gene annotation file, .gz accepted

"""

import os

from docopt import docopt
from ngsutils.gtf import gtf_stem, read_gtf


def main():
    options = docopt(__doc__)
    gtf_path = options['<gtf>']

    output_path = os.path.join(
        options['--output-dir'],
        f"{gtf_stem(gtf_path)}.gene2tx.txt"
    )

    cols = ['gene_id', 'transcript_id']
    annotations = read_gtf(gtf_path, result_type='pandas')
    annotations = annotations.query(
        "feature == 'exon'"
    ).filter(cols).sort_values(by=cols)

    annotations = annotations.drop_duplicates(cols)
    annotations.to_csv(output_path, header=False, index=False, sep='\t')


if __name__ == '__main__':
    main()
