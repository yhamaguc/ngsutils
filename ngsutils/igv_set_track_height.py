#! /usr/bin/env python3

"""
Set panel and track heights in an IGV session to fit the display

Usage:
  igv_set_track_height [options] <xml>

Options:
  -d --display-height <INT>  : Total display height in pixels [default: 900]
  -m --margin <INT>          : Height reserved for the feature panel [default: 300]
  -c --track-class <STR>     : Track classes to resize, comma separated [default: org.broad.igv.sam.AlignmentTrack]
  -n --count <INT>           : Number of panels to divide by [default: auto]
  <xml>                      : IGV session XML file

"""

# NOTE: Intended to run after `igv_hide_tracks`, so that the resized track
#       occupies the whole panel

import sys
import xml.etree.ElementTree as ET

from docopt import docopt

FEATURE_PANEL = "FeaturePanel"


def main():
    options = docopt(__doc__)

    tree = ET.parse(options['<xml>'])
    root = tree.getroot()

    panels = [p for p in root.findall(".//Panel")
              if p.get("name") != FEATURE_PANEL]

    if options['--count'] == 'auto':
        count = len(panels)
    else:
        count = int(options['--count'])

    if count < 1:
        print("No data panel to resize", file=sys.stderr)
        return 1

    display_height = int(options['--display-height'])
    margin = int(options['--margin'])
    height = (display_height - margin) // count

    if height < 1:
        print(f"Display height ({display_height}) minus margin ({margin}) "
              f"leaves no room for {count} panels", file=sys.stderr)
        return 1

    for p in panels:
        p.set("height", str(height))

    for p in root.findall(f".//Panel[@name='{FEATURE_PANEL}']"):
        p.set("height", str(margin // 2))

    for c in options['--track-class'].split(','):
        for t in root.findall(f".//Track[@clazz='{c}']"):
            t.set("height", str(height))

    tree.write(sys.stdout, encoding="unicode", xml_declaration=True)


if __name__ == '__main__':
    main()
