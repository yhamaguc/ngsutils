#! /bin/env python
#
# Remove redundant `chr` prefix from sequence name
#
# Usage:
#   this.py <fasta>
#

import sys


def main():
    fasta = sys.argv[1]
    with open(fasta) as f:
        for l in f:
            # NOTE: Eliminate blank line
            if not l:
                continue

            # NOTE: Remove additional information from sequence name
            if l.startswith(">chr"):
                print(l.replace(">chr", '>'))
            else:
                print(l, end='')


if __name__ == '__main__':
    main()
