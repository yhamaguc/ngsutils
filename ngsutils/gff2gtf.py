#! /usr/bin/env python3
import sys
import os

import sqlite3
from gffparse import read_gtf

"""
Convert gene annotation GFF to SQLite

Usage:
  gff2sqlite <gtf> [target_table_name]

"""


def main():
    gtf_path = sys.argv[1]

    try:
        target_table = sys.argv[2]
    except Exception:
        target_table = 'annotations'


    gtf_df = read_gtf(gtf_path)

    gtf_abspath = os.path.abspath(gtf_path)
    root, _ = os.path.splitext(gtf_abspath)
    output_path = root + '.sqlite'

    conn = sqlite3.connect(output_path)
    gtf_df.to_sql(target_table, conn, if_exists='replace')
    conn.close()
    print("File wrote: {}".format(output_path))


if __name__ == '__main__':
    main()
