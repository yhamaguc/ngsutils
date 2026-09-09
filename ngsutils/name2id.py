#! /usr/bin/env python3

"""
Convert feature (gene/transcript) name to Ensembl ID

Names are matched exactly. With no annotation file given, the map named by
ngsutils/data/maps.json is used: a vendored copy if the package has one, else the
local cache, else the GitHub release asset the manifest points at. A name GENCODE
reuses across identifiers resolves to one of them.

Usage:
  name2id [options]

Options:
  --gtf <PATH>   : Annotation file in GTF format, .gz accepted
  -f --file <PATH>  : Input file [default: stdin]
  -c --col <INT>    : Name column position [default: 1]

"""

import sys
import os
import pickle
import hashlib
import logging

# Suppress all INFO and DEBUG logging from all modules - set EARLY
logging.basicConfig(level=logging.ERROR, force=True)
for logger_name in list(logging.root.manager.loggerDict.keys()):
    logging.getLogger(logger_name).setLevel(logging.ERROR)

from docopt import docopt

import ngsutils.maps as maps

# NOTE: polars and ngsutils.gtf are imported inside main() instead of here. They are
#   reached only when an annotation file is given; on the prebuilt-map path they cost
#   484 ms of the command's 520 ms and nothing uses them. Measured 2026-08-21 by
#   scripts/bench_map_load.py.


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
            f"~/.cache/ngsutils_name2id_{gtf_hash}.pkl")
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)

        if os.path.exists(cache_path):
            try:
                with open(cache_path, 'rb') as f:
                    features = pickle.load(f)
            except:
                features = None

        if features is None:
            print("Generating cache file based on GTF...", file = sys.stderr)
            import ngsutils.gtf as gtfparse
            import polars as pl

            gtf = gtfparse.read_gtf(gtf_path)

            genes = gtf.select([
                "gene_name",
                pl.col("gene_id").str.slice(0, 15)
            ]).unique()

            transcripts = gtf.select([
                "transcript_name",
                pl.col("transcript_id").str.slice(0, 15)
            ]).unique()

            features = dict(genes.iter_rows()) | dict(transcripts.iter_rows())

            with open(cache_path, 'wb') as f:
                pickle.dump(features, f)
    else:
        # NOTE: ngsutils.maps owns which map this is and what its bytes must hash to.
        #   It raises rather than returning a partial map, on purpose: a lookup that
        #   answers nothing for every row is worse than a command that fails.
        try:
            features = maps.load("name2id")
        except maps.MapUnavailable as error:
            raise SystemExit(str(error))

    for l in input_file:
        fields = l.rstrip('\n').split('\t')
        if len(fields) > col_position:
            query = fields[col_position]
            ensembl_id = features.get(query, '')
        else:
            ensembl_id = ''

        output = fields + [ensembl_id]
        print("\t".join(output))


if __name__ == "__main__":
    main()
