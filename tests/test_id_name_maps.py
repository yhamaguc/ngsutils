#!/usr/bin/env python3
"""test_id_name_maps -- the manifest, and the map it resolves to.

The maps are release assets, not committed files, so what the repository holds is
`ngsutils/data/maps.json`: a release, a tag, two file names and their SHA-256.
Three silent failures follow from that arrangement, one per class below.

  **A manifest that points at nothing.** Every field in it is a string, and a wrong
  one fails at the moment a user runs a command with no annotation file -- not at
  install, not at import, not in CI. These tests assert the shape of the manifest
  and that the asset names follow from the release it declares.

  **A map committed by accident.** Committing one puts 23 MB per release into every
  clone forever, which is the cost the release-asset arrangement exists to avoid. A
  vendored pickle is a supported way to run offline, so the check is on what git
  tracks rather than on what the directory holds.

  **A map that resolves to the wrong bytes.** A pickle is opaque and executes on
  load, so ngsutils.maps checks it against the manifest before unpickling. These
  tests cover the check itself, and -- only when a map is actually resolvable
  without the network -- that the two commands answer with it.

    python3 -m unittest discover tests
    pytest tests/test_id_name_maps.py

Nothing here downloads anything. The tests that need a map skip unless one is
already vendored or cached.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ngsutils import maps  # noqa: E402

DATA_DIR = REPO_ROOT / "ngsutils" / "data"

# NOTE: id2name truncates every identifier it looks up to this many characters, so a
#   key longer than that, or one carrying a version suffix, is unreachable.
ID_PREFIX_LENGTH = 15

# Genes read out of gencode.v50.primary_assembly.annotation.gtf.gz on 2026-08-21,
# chosen for being unlikely to be renamed or retired.
STABLE_GENES = {
    "ENSG00000141510": "TP53",
    "ENSG00000012048": "BRCA1",
    "ENSG00000146648": "EGFR",
    "ENSG00000111640": "GAPDH",
}

# NOTE: the canary for a stale map under a fresh release label. GENCODE called this
#   gene C2orf83 in v45 and SLC19A4P in v50, so a v45 map fails this and a v50 map
#   passes it. Replace it, do not delete it, if a later release renames it again.
RENAMED_SINCE_V45 = ("ENSG00000042304", "SLC19A4P")

REQUIRED_ASSET_FIELDS = (
    "file", "archive_bytes", "archive_sha256", "pickle_bytes", "pickle_sha256", "keys")


def resolvable_without_network(kind):
    """Whether load() can answer from a vendored copy or the cache alone."""
    try:
        manifest = maps.manifest()
    except maps.MapUnavailable:
        return False
    release = manifest["release"]
    return (os.path.exists(maps.vendored_path(kind, release))
            or os.path.exists(maps.cache_path(kind, release)))


def run_command(module, lines):
    """The command as a caller runs it: no --gtf, input on stdin."""
    return subprocess.run(
        [sys.executable, "-m", module],
        input="".join(f"{line}\n" for line in lines),
        capture_output=True, text=True, cwd=REPO_ROOT)


# %%
# The manifest, which is all the repository holds

class ManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = maps.manifest()

    def test_the_manifest_declares_one_release_and_both_maps(self):
        self.assertIn("release", self.manifest)
        self.assertEqual(sorted(self.manifest["assets"]), sorted(maps.KINDS))

    def test_every_asset_entry_carries_what_the_check_needs(self):
        for kind in maps.KINDS:
            entry = self.manifest["assets"][kind]
            for field in REQUIRED_ASSET_FIELDS:
                with self.subTest(kind=kind, field=field):
                    self.assertIn(field, entry)
            with self.subTest(kind=kind):
                self.assertEqual(len(entry["pickle_sha256"]), 64)
                self.assertGreater(entry["keys"], 0)

    def test_the_asset_names_follow_from_the_release(self):
        """The regression this file exists for: a release bumped in one field only."""
        release = self.manifest["release"]
        for kind in maps.KINDS:
            with self.subTest(kind=kind):
                self.assertEqual(self.manifest["assets"][kind]["file"],
                                 f"{kind}_{release}.pkl.gz")

    def test_the_download_url_is_the_public_asset_form(self):
        for kind in maps.KINDS:
            url = maps.asset_url(self.manifest["assets"][kind],
                                 self.manifest["repository"], self.manifest["tag"])
            with self.subTest(kind=kind):
                self.assertEqual(
                    url,
                    f"https://github.com/{self.manifest['repository']}/releases/"
                    f"download/{self.manifest['tag']}/"
                    f"{self.manifest['assets'][kind]['file']}")

    def test_no_map_is_committed(self):
        """23 MB per release in every clone forever is what the assets avoid."""
        try:
            completed = subprocess.run(
                ["git", "-C", str(REPO_ROOT), "ls-files", "ngsutils/data"],
                capture_output=True, text=True, check=True)
        except (OSError, subprocess.CalledProcessError):
            self.skipTest("no git to ask what is tracked")
        tracked = [line for line in completed.stdout.split() if line.endswith(".pkl")]
        self.assertEqual(tracked, [], f"maps committed: {tracked}")


# %%
# The check that stands between an opaque file and pickle.loads

class VerificationTests(unittest.TestCase):
    def test_matching_bytes_pass(self):
        payload = b"the bytes the manifest describes"
        maps._verify(payload, maps._sha256(payload), "a fixture")

    def test_altered_bytes_are_refused_by_name(self):
        payload = b"the bytes the manifest describes"
        expected = maps._sha256(payload)
        with self.assertRaises(maps.MapUnavailable) as raised:
            maps._verify(payload + b"!", expected, "a fixture")
        self.assertIn("a fixture", str(raised.exception))
        self.assertIn(expected, str(raised.exception))

    def test_an_unknown_map_is_a_programming_error_not_a_download(self):
        with self.assertRaises(ValueError):
            maps.load("id2names")


# %%
# The map itself, when one is there to be read

@unittest.skipUnless(resolvable_without_network("id2name")
                     and resolvable_without_network("name2id"),
                     "no vendored or cached map; nothing to read without the network")
class ResolvedMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.id2name_map = maps.load("id2name")
        cls.name2id_map = maps.load("name2id")
        cls.manifest = maps.manifest()

    def test_the_map_is_the_size_the_manifest_says(self):
        self.assertEqual(len(self.id2name_map),
                         self.manifest["assets"]["id2name"]["keys"])
        self.assertEqual(len(self.name2id_map),
                         self.manifest["assets"]["name2id"]["keys"])

    def test_no_key_carries_a_version_suffix(self):
        """A key of ENSG00000141510.21 is dead weight: the lookup truncates first."""
        for name, mapping in (("id2name", self.id2name_map),
                              ("name2id", self.name2id_map)):
            longest = max(len(key) for key in mapping if key.startswith("ENS"))
            versioned = [k for k in mapping if k.startswith("ENS") and "." in k]
            with self.subTest(map=name):
                self.assertLessEqual(longest, ID_PREFIX_LENGTH)
                self.assertEqual(versioned[:5], [])

    def test_stable_genes_resolve_both_ways(self):
        for identifier, gene_name in sorted(STABLE_GENES.items()):
            with self.subTest(gene=gene_name):
                self.assertEqual(self.id2name_map.get(identifier), gene_name)
                self.assertEqual(self.name2id_map.get(gene_name), identifier)

    def test_the_map_is_not_the_release_before_it(self):
        """An opaque map relabelled rather than rebuilt passes every test above."""
        identifier, gene_name = RENAMED_SINCE_V45
        self.assertEqual(self.id2name_map.get(identifier), gene_name)

    def test_id2name_answers_from_the_resolved_map(self):
        completed = run_command("ngsutils.id2name", ["ENSG00000141510.21", "not_an_id"])
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stdout.splitlines(),
                         ["ENSG00000141510.21\tTP53", "not_an_id\t"])

    def test_name2id_answers_from_the_resolved_map(self):
        completed = run_command("ngsutils.name2id", ["TP53", "not_a_name"])
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stdout.splitlines(),
                         ["TP53\tENSG00000141510", "not_a_name\t"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
