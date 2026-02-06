#! /usr/bin/env python3
"""
Convert gene annotation SQLite to GTF

Usage:
  sqlite2gtf <sqlite> [<target_table_name>]

Arguments:
  <sqlite>                 SQLite database file
  <target_table_name>      Target table name [default: annotations]

"""

from itertools import compress

import sqlite3

from docopt import docopt


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
    values = ['"{}"'.format(v) for v in values]
    values[-1] = values[-1] + ';'

    for k, v in zip(keys, values):
        if v is not None and v != '':
            list_.append(' '.join((k, v)))

    return list_


def main():
    options = docopt(__doc__)
    sqlite_path = options['<sqlite>']
    target_table = options['<target_table_name>'] or 'annotations'

    conn = sqlite3.connect(sqlite_path)
    conn.row_factory = sqlite3.Row

    cur = conn.cursor()

    # Dynamically determine column order based on table schema
    cur.execute("PRAGMA table_info({})".format(target_table))
    columns = [row[1] for row in cur.fetchall()]

    # Determine order column: use 'index' if exists, otherwise use 'gene_id, transcript_id'
    order_col = '[index]' if 'index' in columns else 'gene_id, transcript_id'

    sql = "SELECT * FROM {} ORDER BY {};".format(target_table, order_col)
    cur.execute(sql)

    keys = []

    for i, r in enumerate(cur):
        keys = keys or r.keys()
        offset = 0
        if keys[0] == 'index':
            offset = 1

        is_not_nulls = [True if v is not None and v != '' else False for v in r[8 + offset:]]
        _values = list(compress(r[8 + offset:], is_not_nulls))
        _keys = list(compress(keys[8 + offset:], is_not_nulls))

        print("{}\t{}".format(
            "\t".join(ifnull(r[offset:8 + offset])),
            '; '.join(attributes_to_str(_values, _keys)))
        )


if __name__ == '__main__':
    main()
