#! /usr/bin/env python

# NOTE: This script was taken from https://github.com/edawson/maf_to_bed and modified

import sys

if __name__ == "__main__":
    slop = 0
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
