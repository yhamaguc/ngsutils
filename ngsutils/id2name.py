#! /usr/bin/env python3
import sys
import os

import pandas as pd
from sqlalchemy import create_engine
from sqlite3 import Row

import ngsutils.utils as utils


"""

Convert Ensembl ID to feature (gene/transcript) name

Usage:
  cat ids.txt | id2name

"""


def main():
    outer_sqlite_path = utils.from_root('assets/gencode.annotation.gtf.sqlite')

    if not os.path.exists(outer_sqlite_path):
        raise FileNotFoundError
        sys.exit(1)

    db_engine = create_engine('sqlite://', echo=False)
    db_engine.row_factory = Row

    df = pd.read_csv(sys.stdin, sep='\t', header=None)
    df.columns = ['c' + str(c) for c in range(len(df.columns))]
    df.to_sql('stdin', con=db_engine, if_exists='replace')

    db_engine.execute(f"attach database '{outer_sqlite_path}' as __ext__;")
    v = "(select distinct gene_id as id, gene_name as name from __ext__.annotations union all select distinct transcript_id as id, transcript_name as name from __ext__.annotations)"
    query = "select {}, ext.* from stdin left join {} as ext "\
            "on substr(stdin.c0, 1, 15) = substr(ext.id, 1, 15) order by stdin.[index];".format(', '.join(df.columns), v)

    results = db_engine.execute(query)

    for row in results:
        print('\t'.join([str(v) if v is not None else '' for v in row]))

    results.close()


if __name__ == '__main__':
    main()
