#!/usr/bin/env python3

"""
Compare two id2name or name2id maps, key by key

What a new GENCODE release does to a lookup table is three numbers -- keys added,
keys that stop answering, and keys whose answer changed -- and only the second and
third can break a caller. This prints them, with examples, and writes the changed
keys out in full so a claim about them can be checked later.

Either side may be a pickle or a gzipped pickle. An earlier release that is no
longer in the working tree can be read straight out of git:

  git show <revision>:ngsutils/data/id2name_gencode.v45.pkl > /tmp/before.pkl

Usage:
  compare_id_name_maps.py --before <PATH> --after <PATH> [--examples <INT>] [--output <PATH>]

Options:
  --before <PATH>   : The older map (required)
  --after <PATH>    : The newer map (required)
  --examples <INT>  : How many examples of each kind to print [default: 5]
  --output <PATH>   : Where a TSV of every changed key is written [default: none]

"""

from __future__ import annotations

import gzip
import os
import pickle

from docopt import docopt


# %%
# Helper functions

def read_map(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return pickle.load(f)
    with open(path, "rb") as f:
        return pickle.load(f)


def show(label, keys, examples, render):
    print(f"{label}: {len(keys)}")
    for key in sorted(keys)[:examples]:
        print(f"  {render(key)}")
    if len(keys) > examples:
        print(f"  ... and {len(keys) - examples} more")


# %%
# Main logic

def main():
    options = docopt(__doc__)
    examples = int(options["--examples"])

    before = read_map(options["--before"])
    after = read_map(options["--after"])

    added = set(after) - set(before)
    removed = set(before) - set(after)
    changed = {k for k in set(before) & set(after) if before[k] != after[k]}

    print(f"before {options['--before']}: {len(before)} keys")
    print(f"after  {options['--after']}: {len(after)} keys\n")
    show("added", added, examples, lambda k: f"{k} -> {after[k]}")
    show("removed", removed, examples, lambda k: f"{k} was {before[k]}")
    show("changed", changed, examples, lambda k: f"{k}: {before[k]} -> {after[k]}")

    if options["--output"] == "none":
        return
    output_path = os.path.abspath(options["--output"])
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write("kind\tkey\tbefore\tafter\n")
        for key in sorted(removed):
            f.write(f"removed\t{key}\t{before[key]}\t\n")
        for key in sorted(changed):
            f.write(f"changed\t{key}\t{before[key]}\t{after[key]}\n")
        for key in sorted(added):
            f.write(f"added\t{key}\t\t{after[key]}\n")
    print(f"\nwrote {output_path}")


if __name__ == "__main__":
    main()
