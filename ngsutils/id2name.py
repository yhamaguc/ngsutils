#! /usr/bin/env python3

"""
Convert Ensembl ID to feature (gene/transcript) name

Usage:
  id2name [options]

Options:
  --db <PATH>    : Annotation data base (SQLite)
  --file <PATH>  : Input file [default: stdin]
  --col <INT>    : ID column position [default: 1]

"""

import sys
import os

from docopt import docopt
import pandas as pd
from sqlalchemy import create_engine
from sqlite3 import Row

import ngsutils.utils as utils


def main():
    options = docopt(__doc__)

    outer_sqlite_path = (
        options["--db"]
        if options["--db"]
        else utils.from_root("assets/gencode.annotation.gtf.sqlite")
    )

    col_position = int(options["--col"]) - 1

    if not os.path.exists(outer_sqlite_path):
        raise FileNotFoundError
        sys.exit(1)

    db_engine = create_engine("sqlite://", echo=False)
    db_engine.row_factory = Row

    input = sys.stdin if options["--file"] == 'stdin' else open(options["--file"])

    df = pd.read_csv(input, sep="\t", header=None)
    df.columns = ["c" + str(c) for c in range(len(df.columns))]
    df.to_sql("stdin", con=db_engine, if_exists="replace")

    db_engine.execute(f"attach database '{outer_sqlite_path}' as __ext__;")
    v = "(select distinct gene_id as id, gene_name as name from __ext__.annotations union all select distinct transcript_id as id, transcript_name as name from __ext__.annotations)"
    query = (
        "select {}, ext.name from stdin left join {} as ext "
        "on substr(stdin.c{}, 1, 15) = substr(ext.id, 1, 15) order by stdin.[index];".format(
            ", ".join(df.columns), v, str(col_position)
        )
    )

    results = db_engine.execute(query)

    for row in results:
        print("\t".join([str(v) if v is not None else "" for v in row]))

    results.close()


if __name__ == "__main__":
    main()
