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
    # TODO: If sys.stdin is empty: exit
    outer_sqlite_path = utils.from_root('assets/gencode.v29.annotation.sqlite')
    if not os.path.exists(outer_sqlite_path):
        raise FileNotFoundError
        sys.exit(1)

    db_engine = create_engine('sqlite://', echo=False)
    db_engine.row_factory = Row

    USECOLS = [0]
    df = pd.read_csv(sys.stdin, sep='\t', header=None, usecols=USECOLS)
    df.columns = ['c' + str(c) for c in USECOLS]
    df.to_sql('stdin', con=db_engine, if_exists='replace')

    db_engine.execute("attach database '{}' as __outer__;".format(outer_sqlite_path))
    query = "select outer.* from stdin left join __outer__.features as outer "\
            "on stdin.c0 = outer.feature_id order by stdin.[index];"
    results = db_engine.execute(query)

    for row in results:
        print('\t'.join(row.values()))

    results.close()


if __name__ == '__main__':
    main()
