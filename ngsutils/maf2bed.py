#!/usr/bin/env python3
"""
Convert MAF to BED

Usage:
  maf2bed [options]

Options:
  -s --slop <INT>  : Slop size [default: 0]

"""

# NOTE: This script was taken from https://github.com/edawson/maf_to_bed and modified

import sys

from docopt import docopt


def main():
    options = docopt(__doc__)
    slop = int(options['--slop'])
    for line in sys.stdin:
        if line .startswith("#") or line.startswith("individual") or line.startswith("Hugo_Symbol"):
            continue

        tokens = line.strip().split("\t")

        chr_1, start, end = tokens.pop(4), int(
            tokens.pop(4)), int(tokens.pop(4))

        tokens = [f'"{x}"' if x != "NA" else '' for x in tokens]

        line_one = "\t".join(
            [chr_1, str(start - slop - 1), str(end + slop), ";".join(tokens)])
        print(line_one)

if __name__ == '__main__':
    main()