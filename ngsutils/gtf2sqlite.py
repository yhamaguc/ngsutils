#! /usr/bin/env python3
import sys
import os

import sqlite3
from gtfparse import read_gtf

"""
Convert gene annotation GTF to SQLite

Usage:
  gtf2sqlite <gtf> [target_table_name]

"""


def main():
    gtf_file = sys.argv[1]
    try:
        target_table = sys.argv[2]
    except Exception:
        target_table = 'annotations'

    gtf_df = read_gtf(gtf_file)

    output_dir = os.path.dirname(gtf_file)
    output_file, _ = os.path.splitext(os.path.basename(gtf_file))
    if output_dir == '':
        output_dir = '.'

    conn = sqlite3.connect(output_dir + '/' + output_file + '.sqlite')
    gtf_df.to_sql(target_table, conn, if_exists='replace')
    conn.close()
    print('File wrote: ' + output_dir + '/' + output_file + '.sqlite')


if __name__ == '__main__':
    main()
