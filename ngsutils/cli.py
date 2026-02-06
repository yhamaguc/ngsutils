"""#!/usr/bin/env python3Subcommand dispatcher for ngsutils.

Usage:
    ngsutils <subcommand> [<args>...]
    ngsutils (-h | --help)

Options:
    -h --help  Show this help message.
"""

import importlib
import sys

from docopt import docopt


_SUBCOMMANDS = {
    "abs2rel": "ngsutils.abs2rel:main",
    "extract_splice_sites": "ngsutils.extract_splice_sites:main",
    "gene2tx": "ngsutils.gene2tx:main",
    "get_biomart": "ngsutils.get_biomart:main",
    "gtf2bed": "ngsutils.gtf2bed:main",
    "gtf2sqlite": "ngsutils.gtf2sqlite:main",
    "gtf2tsv": "ngsutils.gtf2tsv:main",
    "id2name": "ngsutils.id2name:main",
    "maf2bed": "ngsutils.maf2bed:main",
    "name2id": "ngsutils.name2id:main",
    "sqlite2gtf": "ngsutils.sqlite2gtf:main",
}


def _print_help(prog):
    print(f"Usage: {prog} <subcommand> [args...]\n")
    print("Run `ngsutils <subcommand> --help` for command help.\n")
    print("Available subcommands:")
    for name in sorted(_SUBCOMMANDS):
        print(f"  {name}")


def _load_callable(target):
    module_name, func_name = target.split(":", 1)
    module = importlib.import_module(module_name)
    return getattr(module, func_name)


def main():
    opts = docopt(__doc__, options_first=True, help=False)
    if opts.get("--help"):
        _print_help("ngsutils")
        return 0

    subcommand_raw = opts.get("<subcommand>")
    if not subcommand_raw:
        _print_help("ngsutils")
        return 0

    target = _SUBCOMMANDS.get(subcommand_raw)
    if not target:
        print(f"Unknown subcommand: {subcommand_raw}\n", file=sys.stderr)
        _print_help("ngsutils")
        return 2

    func = _load_callable(target)
    sys.argv = [f"ngsutils {subcommand_raw}"] + (opts.get("<args>") or [])
    return func()


if __name__ == "__main__":
    raise SystemExit(main())
