"""
ngsutils.maps
~~~~~~~~~~~~~

The id2name / name2id lookup maps, and where they come from.

The maps are too large to commit -- 23 MB each for GENCODE v50, growing with every
release -- so they are GitHub release assets instead. `ngsutils/data/maps.json` is
committed in their place: it names the release, the tag, the asset file names, and
the SHA-256 of both the archive and the pickle inside it. bin/build_id_name_maps.py
writes that manifest and the archives together.

load() resolves a map in this order, and says on stderr which one it took:

  1. a plain pickle vendored into ngsutils/data, for an install with no network
  2. the local cache, under XDG_CACHE_HOME or ~/.cache
  3. the release asset, downloaded and cached

NOTE: this module is the only place that decides what a map is. The SHA-256 in the
  manifest is checked against the downloaded bytes *before* they are unpickled, and
  against the cached bytes every time they are read, so a corrupted or substituted
  file is an error rather than a wrong answer. Callers do not repeat the check.
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import sys
import urllib.error
import urllib.request

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
MANIFEST_PATH = os.path.join(DATA_DIR, "maps.json")

KINDS = ("id2name", "name2id")

DOWNLOAD_TIMEOUT_SECONDS = 300
READ_BLOCK_BYTES = 1 << 20


class MapUnavailable(Exception):
    """No route to the map: nothing vendored, nothing cached, nothing downloadable."""


def _log(message):
    print(message, file=sys.stderr, flush=True)


def manifest():
    """The committed manifest, or an error naming what is missing."""
    if not os.path.exists(MANIFEST_PATH):
        raise MapUnavailable(
            f"No manifest at {MANIFEST_PATH}. The package is incomplete; reinstall it "
            "or pass an annotation file on the command line.")
    with open(MANIFEST_PATH) as f:
        return json.load(f)


def cache_dir():
    root = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
    return os.path.join(root, "ngsutils", "maps")


def vendored_path(kind, release):
    """Where a plain pickle would sit if someone put one in the package by hand."""
    return os.path.join(DATA_DIR, f"{kind}_{release}.pkl")


def cache_path(kind, release):
    """The cache holds the unpacked pickle, so a warm run pays no decompression."""
    return os.path.join(cache_dir(), f"{kind}_{release}.pkl")


def asset_url(entry, repository, tag):
    """The public download URL of a release asset.

    A 302 to the asset host is expected here and urllib follows it. Public
    repositories need no credentials; verified against the API's browser_download_url
    on 2026-08-21.
    """
    return f"https://github.com/{repository}/releases/download/{tag}/{entry['file']}"


def _sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def _read(path):
    with open(path, "rb") as f:
        return f.read()


def _verify(payload, expected_sha256, what):
    found = _sha256(payload)
    if found != expected_sha256:
        raise MapUnavailable(
            f"{what} does not match the manifest: expected SHA-256 {expected_sha256}, "
            f"found {found}. Refusing to unpickle it.")


def _download(url, expected_bytes):
    _log(f"Downloading {expected_bytes / 1e6:.1f} MB: {url}")
    try:
        with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response:
            payload = response.read()
    except (urllib.error.URLError, OSError) as error:
        raise MapUnavailable(
            f"Cannot download {url}: {error}. Pass an annotation file on the command "
            "line instead, or fetch that file by hand into "
            f"{cache_dir()} after unpacking it.") from error

    if len(payload) != expected_bytes:
        raise MapUnavailable(
            f"{url} returned {len(payload)} bytes, the manifest says {expected_bytes}.")
    return payload


def _write_cache(path, payload):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    temporary_path = f"{path}.{os.getpid()}.tmp"
    with open(temporary_path, "wb") as f:
        f.write(payload)
    os.replace(temporary_path, path)


def load(kind):
    """The map named by `kind`, from whichever of the three routes has it.

    Raises MapUnavailable rather than returning a partial or unverified map.
    """
    if kind not in KINDS:
        raise ValueError(f"Unknown map {kind!r}, expected one of {KINDS}")

    entry_manifest = manifest()
    release = entry_manifest["release"]
    entry = entry_manifest["assets"][kind]

    vendored = vendored_path(kind, release)
    if os.path.exists(vendored):
        _log(f"Reading the vendored map {vendored}")
        payload = _read(vendored)
        _verify(payload, entry["pickle_sha256"], vendored)
        return pickle.loads(payload)

    cached = cache_path(kind, release)
    if os.path.exists(cached):
        payload = _read(cached)
        _verify(payload, entry["pickle_sha256"], cached)
        return pickle.loads(payload)

    import gzip

    url = asset_url(entry, entry_manifest["repository"], entry_manifest["tag"])
    archive = _download(url, entry["archive_bytes"])
    _verify(archive, entry["archive_sha256"], url)

    payload = gzip.decompress(archive)
    _verify(payload, entry["pickle_sha256"], f"the pickle inside {url}")

    _write_cache(cached, payload)
    _log(f"Cached {cached}")
    return pickle.loads(payload)
