#!/usr/bin/env python3
"""test_entry_points — the command table and the docstrings docopt reads it through.

Two silent failures, one per class below.

  **A command that stops existing.** Until 0.2.0 every utility was its own console script
  (`id2name`, `gene2tx`, …) and the move from setup.cfg to pyproject.toml kept only
  `ngsutils`. Nothing failed at build or install time; the commands simply were not there
  on the next `pip install .`, and a caller that had them on PATH got "command not found"
  from a package that installed cleanly. `ngsutils/cli.py` owns the subcommand table, so
  these tests assert `[project.scripts]` matches it — a subcommand added in one place only
  fails here rather than at the next install.

  **A `[default:]` that evaporates.** docopt 0.6.2 splits the *whole* docstring on
  dash-initial lines, so a prose line beginning with a flag, a dash-underlined heading, or
  one option's description naming another option redefines that option and drops its
  default. Nothing raises: the value arrives as `None` or `''` and is read downstream as a
  missing argument. The two tests the failure needs are the intended key/default table and
  the cross-reference scan, run over every command.

    python3 -m unittest discover tests
    pytest tests/test_entry_points.py

Python 3.11+ for `tomllib`. Needs `ngsutils` importable and the repository's pyproject.toml
readable; nothing is installed, run, or written.
"""

from __future__ import annotations

import ast
import importlib
import os
import re
import sys
import tomllib
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from docopt import parse_defaults, printable_usage  # noqa: E402

from ngsutils.cli import _SUBCOMMANDS  # noqa: E402

PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"

# The dispatcher itself, which is a console script but not a subcommand.
DISPATCHER = "ngsutils"
DISPATCHER_TARGET = "ngsutils.cli:main"

# NOTE: the intended command line of every subcommand: long option -> the default docopt
#   must deliver for it (False for a flag, None where no default is declared). Written out
#   rather than derived from the docstring, because the docstring is what is under test --
#   deriving the expectation from it would agree with any breakage.
EXPECTED_OPTIONS = {
    "abs2rel": {},
    "extract_splice_sites": {"--width": "2"},
    "gene2tx": {"--output-dir": "."},
    "get_biomart": {"--reference": "grch38"},
    "gtf2bed": {"--tx-only": False, "--simplify": False},
    "gtf2sqlite": {},
    "gtf2tsv": {},
    "id2name": {"--gtf": None, "--file": "stdin", "--col": "1"},
    "igv_hide_tracks": {"--track-class": "org.broad.igv.sam.CoverageTrack"},
    "igv_list_tracks": {"--all": False},
    "igv_rename_tracks": {},
    "igv_set_track_height": {
        "--display-height": "900",
        "--margin": "300",
        "--track-class": "org.broad.igv.sam.AlignmentTrack",
        "--count": "auto",
    },
    "maf2bed": {"--slop": "0", "--help": False},
    "name2id": {"--gtf": None, "--file": "stdin", "--col": "1"},
    "sqlite2gtf": {},
    "write_coverage_track": {
        "--output-dir": None,
        "--gtf": None,
        "--id": None,
        "--region": None,
        "--mode": "span",
        "--lower": "0",
        "--upper": "0",
        "--flank": None,
        "--label": None,
        "--label-from": "stem",
        "--no-track-line": False,
        "--min-mq": None,
        "--excl-flags": None,
        "--threads": None,
        "--quiet": False,
        "--help": False,
        "--version": False,
    },
}


def console_scripts():
    with open(PYPROJECT_PATH, "rb") as f:
        return tomllib.load(f)["project"]["scripts"]


def module_doc(subcommand):
    module_name = _SUBCOMMANDS[subcommand].split(":", 1)[0]
    return importlib.import_module(module_name).__doc__ or ""


def option_descriptions(doc):
    """(option definition, description) per option, split the way docopt 0.6.2 splits."""
    split = re.split(r"\n *(<\S+?>|-\S+?)", doc)[1:]
    chunks = [head + tail for head, tail in zip(split[::2], split[1::2])]
    pairs = []
    for chunk in chunks:
        if not chunk.startswith("-"):
            continue
        definition, _, description = chunk.partition("  ")
        pairs.append((definition.strip(), description.strip()))
    return pairs


def usage_patterns(doc):
    """The Usage section's pattern lines, continuation lines folded into the line above."""
    lines = printable_usage(doc).splitlines()[1:]
    patterns = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        indent = len(line) - len(line.lstrip())
        if patterns and indent > patterns[-1][0]:
            patterns[-1] = (patterns[-1][0], patterns[-1][1] + " " + stripped)
        else:
            patterns.append((indent, stripped))
    return [pattern for _, pattern in patterns]


# %%
# The command table, which pyproject.toml and cli.py have to agree on

class CommandTableTests(unittest.TestCase):
    def test_every_subcommand_is_installed_under_its_own_name(self):
        """The regression this file exists for: `id2name` and friends were console scripts
        before 0.2.0, and dropping them broke every caller that used the bare name."""
        scripts = console_scripts()
        missing = sorted(set(_SUBCOMMANDS) - set(scripts))
        self.assertEqual(missing, [], f"no console script for: {missing}")

    def test_no_console_script_is_absent_from_the_subcommand_table(self):
        """The other direction: a script pointing at a module the dispatcher does not list
        is a command reachable one way and not the other."""
        scripts = console_scripts()
        extra = sorted(set(scripts) - set(_SUBCOMMANDS) - {DISPATCHER})
        self.assertEqual(extra, [], f"console script with no subcommand: {extra}")

    def test_the_two_routes_reach_the_same_function(self):
        scripts = console_scripts()
        self.assertEqual(scripts[DISPATCHER], DISPATCHER_TARGET)
        for subcommand, target in sorted(_SUBCOMMANDS.items()):
            with self.subTest(subcommand=subcommand):
                self.assertEqual(scripts[subcommand], target)

    def test_every_target_imports_and_is_callable(self):
        for subcommand, target in sorted(_SUBCOMMANDS.items()):
            with self.subTest(subcommand=subcommand):
                module_name, function_name = target.split(":", 1)
                module = importlib.import_module(module_name)
                self.assertTrue(callable(getattr(module, function_name)))


# %%
# The docstrings, which are the command line docopt actually parses

class DocstringTests(unittest.TestCase):
    def test_the_option_keys_are_exactly_the_intended_set(self):
        """A dash-initial line anywhere in a docstring adds an option key; a description
        naming another option removes one. Either way this set changes."""
        for subcommand in sorted(_SUBCOMMANDS):
            with self.subTest(subcommand=subcommand):
                parsed = {o.long for o in parse_defaults(module_doc(subcommand)) if o.long}
                self.assertEqual(parsed, set(EXPECTED_OPTIONS[subcommand]))

    def test_the_defaults_survived_the_docstring(self):
        for subcommand in sorted(_SUBCOMMANDS):
            doc = module_doc(subcommand)
            parsed = {o.long: o.value for o in parse_defaults(doc) if o.long}
            for option, default in sorted(EXPECTED_OPTIONS[subcommand].items()):
                with self.subTest(subcommand=subcommand, option=option):
                    self.assertEqual(parsed.get(option), default)

    def test_no_option_description_names_another_declared_option(self):
        """`--flank  Shorthand for --lower BP` costs `--lower` its default, because docopt
        starts a new option definition at the dash it finds in the description."""
        for subcommand in sorted(_SUBCOMMANDS):
            doc = module_doc(subcommand)
            declared = {o.long for o in parse_defaults(doc) if o.long}
            for definition, description in option_descriptions(doc):
                named = set(re.findall(r"--[A-Za-z][\w-]*", description)) & declared
                named -= set(re.findall(r"--[A-Za-z][\w-]*", definition))
                with self.subTest(subcommand=subcommand, option=definition):
                    self.assertEqual(sorted(named), [], f"{definition} names {sorted(named)}")

    def test_every_usage_pattern_names_the_command(self):
        """The Usage section is the help text a caller reads, and it is written by hand:
        `abspos2relpos <gtf>` outlived the rename to `abs2rel` in it for two years."""
        for subcommand in sorted(_SUBCOMMANDS):
            for pattern in usage_patterns(module_doc(subcommand)):
                with self.subTest(subcommand=subcommand, pattern=pattern):
                    self.assertEqual(pattern.split()[0], subcommand)

    def test_no_option_is_listed_in_more_than_one_usage_pattern(self):
        """Registering an option twice is a DocoptLanguageError at import -- "is not a
        unique prefix" -- so this one is loud, but only for whoever imports the module."""
        for subcommand in sorted(_SUBCOMMANDS):
            seen = {}
            for pattern in usage_patterns(module_doc(subcommand)):
                for option in set(re.findall(r"--[A-Za-z][\w-]*", pattern)):
                    seen.setdefault(option, []).append(pattern)
            for option, patterns in sorted(seen.items()):
                with self.subTest(subcommand=subcommand, option=option):
                    self.assertEqual(len(patterns), 1, f"{option} in {len(patterns)} patterns")


# %%
# The standalone scripts, which are not subcommands but parse the same way

# NOTE: bin/ and scripts/ hold CLI scripts that the subcommand table does not know
#   about, so the loop above cannot reach them. The docopt traps bite them the same
#   way -- a description naming another option costs it its default -- so the same
#   two assertions are made here against the docstring read off disk.
STANDALONE_EXPECTED_OPTIONS = {
    "bin/build_id_name_maps.py": {
        "--gtf": None,
        "--release": "auto",
        "--tag": "auto",
        "--repository": "auto",
        "--output-dir": "auto",
        "--manifest": "auto",
        "--cross-check": False,
    },
    "scripts/bench_map_load.py": {
        "--map": None,
        "--rounds": "15",
        "--output": "none",
        "--skip-commands": False,
    },
    "scripts/bench_build_routes.py": {
        "--gtf": None,
        "--route": "both",
        "--output": "none",
    },
    "scripts/package_map_archives.py": {
        "--release": None,
        "--id2name": None,
        "--name2id": None,
        "--recovered-from": "unknown",
        "--tag": "auto",
        "--repository": "auto",
        "--output-dir": "auto",
    },
    "scripts/compare_id_name_maps.py": {
        "--before": None,
        "--after": None,
        "--examples": "5",
        "--output": "none",
    },
}


def standalone_doc(relative_path):
    """The module docstring, read without importing -- these scripts are not modules."""
    source = (REPO_ROOT / relative_path).read_text()
    return ast.get_docstring(ast.parse(source)) or ""


class StandaloneScriptTests(unittest.TestCase):
    def test_every_expected_script_is_there(self):
        for relative_path in sorted(STANDALONE_EXPECTED_OPTIONS):
            with self.subTest(script=relative_path):
                self.assertTrue((REPO_ROOT / relative_path).exists())

    def test_the_option_keys_are_exactly_the_intended_set(self):
        for relative_path, expected in sorted(STANDALONE_EXPECTED_OPTIONS.items()):
            parsed = {o.long for o in parse_defaults(standalone_doc(relative_path))
                      if o.long}
            with self.subTest(script=relative_path):
                self.assertEqual(parsed, set(expected))

    def test_the_defaults_survived_the_docstring(self):
        for relative_path, expected in sorted(STANDALONE_EXPECTED_OPTIONS.items()):
            doc = standalone_doc(relative_path)
            parsed = {o.long: o.value for o in parse_defaults(doc) if o.long}
            for option, default in sorted(expected.items()):
                with self.subTest(script=relative_path, option=option):
                    self.assertEqual(parsed.get(option), default)

    def test_no_option_description_names_another_declared_option(self):
        for relative_path in sorted(STANDALONE_EXPECTED_OPTIONS):
            doc = standalone_doc(relative_path)
            declared = {o.long for o in parse_defaults(doc) if o.long}
            for definition, description in option_descriptions(doc):
                named = set(re.findall(r"--[A-Za-z][\w-]*", description)) & declared
                named -= set(re.findall(r"--[A-Za-z][\w-]*", definition))
                with self.subTest(script=relative_path, option=definition):
                    self.assertEqual(sorted(named), [],
                                     f"{definition} names {sorted(named)}")

    def test_the_usage_section_is_one_pattern_naming_the_script(self):
        for relative_path in sorted(STANDALONE_EXPECTED_OPTIONS):
            patterns = usage_patterns(standalone_doc(relative_path))
            with self.subTest(script=relative_path):
                self.assertEqual(len(patterns), 1, patterns)
                self.assertEqual(patterns[0].split()[0],
                                 os.path.basename(relative_path))


if __name__ == "__main__":
    unittest.main(verbosity=2)
