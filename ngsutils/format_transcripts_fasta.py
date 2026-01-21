#! /bin/env python
# Format sequence name for kallisto indexing
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
            if l.startswith('>'):
                print(l.split()[0].split('|')[0].split(' ')[0])
            else:
                print(l, end='')


if __name__ == '__main__':
    main()
