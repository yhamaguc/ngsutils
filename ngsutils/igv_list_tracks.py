#! /usr/bin/env python3

"""
List track names and ids in an IGV session

Usage:
  igv_list_tracks [options] <xml>

Options:
  -a --all  : List all tracks including hidden ones
  <xml>     : IGV session XML file

"""

import xml.etree.ElementTree as ET

from docopt import docopt


def main():
    options = docopt(__doc__)

    root = ET.parse(options['<xml>']).getroot()
    query = ".//Track" if options['--all'] else ".//Track[@visible='true']"

    for t in root.findall(query):
        print(f"{t.get('name')}\t{t.get('id')}")


if __name__ == '__main__':
    main()
