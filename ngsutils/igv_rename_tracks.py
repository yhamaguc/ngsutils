#! /usr/bin/env python3

"""
Rename visible tracks in an IGV session

Usage:
  igv_rename_tracks [options] <xml> <names>

Options:
  <xml>    : IGV session XML file
  <names>  : Text file with one track name per line, ordered as `igv_list_tracks`

"""

import sys
import xml.etree.ElementTree as ET

from docopt import docopt


def main():
    options = docopt(__doc__)

    tree = ET.parse(options['<xml>'])
    tracks = tree.getroot().findall(".//Track[@visible='true']")

    with open(options['<names>']) as f:
        names = [l.rstrip() for l in f]

    if len(tracks) != len(names):
        print(f"Number of names ({len(names)}) does not match "
              f"number of visible tracks ({len(tracks)})", file=sys.stderr)
        return 1

    for (t, n) in zip(tracks, names):
        t.set("name", n)

    tree.write(sys.stdout, encoding="unicode", xml_declaration=True)


if __name__ == '__main__':
    main()
