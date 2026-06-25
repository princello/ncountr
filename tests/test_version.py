"""Tests that the package version is single-sourced and consistent.

The version is declared once in ``pyproject.toml`` and surfaced through
``importlib.metadata``; ``ncountr.__version__`` reads it back.  These tests
guard against the two copies drifting apart (the historical 0.1.0 / 0.2.0
mismatch).
"""

from __future__ import annotations

import re
from importlib.metadata import version as pkg_version
from pathlib import Path

import pytest

import ncountr


PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"


def _pyproject_version() -> str | None:
    """Extract [project] version from pyproject.toml via a simple regex.

    Avoids a hard dependency on tomllib (absent on Python 3.9/3.10).
    """
    if not PYPROJECT.exists():
        return None
    for line in PYPROJECT.read_text().splitlines():
        m = re.match(r"""\s*version\s*=\s*["']([^"']+)["']""", line)
        if m:
            return m.group(1)
    return None


def test_version_is_a_nonempty_string():
    assert isinstance(ncountr.__version__, str)
    assert ncountr.__version__


def test_version_matches_installed_metadata():
    assert ncountr.__version__ == pkg_version("ncountr")


def test_version_is_not_the_uninstalled_fallback():
    """If the package is importable from metadata, the fallback must not trigger."""
    assert ncountr.__version__ != "0.0.0+unknown"


@pytest.mark.skipif(_pyproject_version() is None, reason="pyproject.toml not available")
def test_version_matches_pyproject():
    assert ncountr.__version__ == _pyproject_version()
