#!/usr/bin/env python3

"""
Convert gene annotation GTF to SQLite

Usage:
  gtf2sqlite <gtf> [<target_table_name>]

Arguments:
  <gtf>                    GTF file, .gz accepted
  <target_table_name>      Target table name [default: annotations]

"""

import sys
import os

from docopt import docopt
import ngsutils.gtf as gtfparse


def main():
    options = docopt(__doc__)
    gtf_path = options['<gtf>']
    target_table = options['<target_table_name>'] or 'annotations'

    print(f"gtfparse.version: {gtfparse.__version__}", file=sys.stderr)

    gtf_df = gtfparse.read_gtf(gtf_path)

    output_dir = os.path.dirname(os.path.abspath(gtf_path))
    output_path = os.path.join(
        output_dir, gtfparse.gtf_stem(gtf_path) + '.sqlite')

    conn = 'sqlite:///' + output_path
    gtf_df.write_database(target_table, conn, if_table_exists='replace', engine='sqlalchemy')
    print("File wrote: {}".format(output_path))


if __name__ == '__main__':
    main()
