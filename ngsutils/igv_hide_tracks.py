#! /usr/bin/env python3

"""
Hide tracks of the given classes in an IGV session

Usage:
  igv_hide_tracks [options] <xml>

Options:
  -c --track-class <STR>  : Track classes to hide, comma separated [default: org.broad.igv.sam.CoverageTrack]
  <xml>                   : IGV session XML file

"""

import sys
import xml.etree.ElementTree as ET

from docopt import docopt


def main():
    options = docopt(__doc__)

    tree = ET.parse(options['<xml>'])
    root = tree.getroot()

    for c in options['--track-class'].split(','):
        for t in root.findall(f".//Track[@clazz='{c}']"):
            t.set("visible", "false")

    tree.write(sys.stdout, encoding="unicode", xml_declaration=True)


if __name__ == '__main__':
    main()
