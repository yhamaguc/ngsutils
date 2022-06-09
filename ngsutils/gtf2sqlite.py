#! /usr/bin/env python3

"""
Convert gene annotation GTF to SQLite

Usage:
  gtf2sqlite <gtf> [target_table_name]

"""

import sys
import os

import sqlite3
import gtfparse


def main():
    gtf_path = sys.argv[1]
    try:
        target_table = sys.argv[2]
    except Exception:
        target_table = 'annotations'

    print(f"gtfparse.version: {gtfparse.__version__}", file=sys.stderr)

    gtf_df = gtfparse.read_gtf(gtf_path)

    gtf_abspath = os.path.abspath(gtf_path)
    root, _ = os.path.splitext(gtf_abspath)
    output_path = root + '.sqlite'

    conn = sqlite3.connect(output_path)
    gtf_df.to_sql(target_table, conn, if_exists='replace')
    conn.close()
    print("File wrote: {}".format(output_path))


if __name__ == '__main__':
    main()
