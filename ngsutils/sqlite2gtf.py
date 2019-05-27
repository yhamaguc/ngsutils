#! /usr/bin/env python3
import sys

import sqlite3


"""
Convert gene annotation SQLite to GTF

Usage:
  sqlite2gtf <sqlite> [target_table_name]

"""


def ifnull(values: list, replacement=''):
    list_ = []

    for v in values:
        value = str(v) if (v is not None and v != '') else str(replacement).strip()
        list_.append(value)

    return list_


def attributes_to_str(attribute_values: list, attribute_keys: list, replacement=''):
    values = attribute_values
    keys = attribute_keys
    list_ = []

    values = ifnull(values, replacement)
    values[-1] = values[-1] + ';'

    for k, v in zip(keys, values):
        if v is not None and v != '':
            list_.append(' '.join((k, v)))

    return list_


def main():
    sqlite_path = sys.argv[1]
    try:
        target_table = sys.argv[2]
    except Exception:
        target_table = 'annotations'

    conn = sqlite3.connect(sqlite_path)
    conn.row_factory = sqlite3.Row

    cur = conn.cursor()

    cur.execute("select * from {} order by gene_id, transcript_id, exon_number, exon_id;".format(target_table))
    keys = []

    for i, r in enumerate(cur):
        keys = keys or r.keys()
        print("{}\t{}".format(
            "\t".join(ifnull(r[:8])),
            '; '.join(attributes_to_str(r[8:], keys[8:])))
        )


if __name__ == '__main__':
    main()
