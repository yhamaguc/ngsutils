#! /usr/bin/env python3

"""
Convert Ensembl ID to feature (gene/transcript) name

Usage:
  id2name [options]

Options:
  --gtf <PATH>   : Annotation file in GTF format
  --file <PATH>  : Input file [default: stdin]
  --col <INT>    : ID column position [default: 1]

"""

import sys
import os
import pickle
import hashlib
import time
import logging

# Suppress all INFO and DEBUG logging from all modules - set EARLY
logging.basicConfig(level=logging.ERROR, force=True)
for logger_name in list(logging.root.manager.loggerDict.keys()):
    logging.getLogger(logger_name).setLevel(logging.ERROR)

import gtfparse
import polars as pl

from docopt import docopt


def main():
    options = docopt(__doc__)

    gtf_path = options["--gtf"]

    col_position = int(options["--col"]) - 1

    input_file = sys.stdin if options["--file"] == 'stdin' else open(
        options["--file"])

    features = None

    if gtf_path:
        if not os.path.exists(gtf_path):
            print("Please specify a valid GTF file.", file=sys.stderr)
            raise FileNotFoundError

        gtf_stat = os.stat(gtf_path)
        gtf_hash = hashlib.md5(
            f"{gtf_path}:{gtf_stat.st_mtime}".encode()).hexdigest()
        cache_path = os.path.expanduser(
            f"~/.cache/ngsutils_id2name_{gtf_hash}.pkl")
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)

        if os.path.exists(cache_path):
            try:
                with open(cache_path, 'rb') as f:
                    features = pickle.load(f)
            except:
                features = None

        if features is None:
            print("Generating cache file based on GTF...", file = sys.stderr)
            gtf = gtfparse.read_gtf(gtf_path)

            genes = gtf.select([
                pl.col("gene_id").str.slice(0, 15),
                "gene_name"
            ]).unique()

            transcripts = gtf.select([
                pl.col("transcript_id").str.slice(0, 15),
                "transcript_name"
            ]).unique()

            features = dict(genes.iter_rows()) | dict(transcripts.iter_rows())

            with open(cache_path, 'wb') as f:
                pickle.dump(features, f)
    else:
        default_pkl = os.path.join(os.path.dirname(__file__), 'data', 'id2name_gencode.v45.pkl')

        if os.path.exists(default_pkl):
            try:
                with open(default_pkl, 'rb') as f:
                    features = pickle.load(f)
            except:
                features = None

        if features is None:
            print("Error: No GTF file specified and no default annotation data found.", file=sys.stderr)
            print("Please specify a valid GTF file as --gtf option or GTF environment variable.", file=sys.stderr)
            raise FileNotFoundError

    for line in input_file:
        fields = line.rstrip('\n').split('\t')
        if len(fields) > col_position:
            query = fields[col_position][:15]
            name = features.get(query, '')
        else:
            name = ''

        output = fields + [name]
        print("\t".join(output))


if __name__ == "__main__":
    main()
