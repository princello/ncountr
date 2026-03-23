"""Download nCounter RCC files from NCBI GEO."""

from __future__ import annotations

import gzip
import io
import os
import tarfile
import tempfile
from pathlib import Path
from typing import Union
from urllib.request import urlretrieve
from urllib.error import URLError


def _gse_to_ftp_dir(gse: str) -> str:
    """Convert a GSE accession to the GEO FTP supplement directory URL."""
    gse = gse.upper()
    # GEO organizes by GSEnnn (first digits, thousands)
    prefix = gse[:len(gse) - 3] + "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/"
    )


def _download_with_progress(url: str, dest: Path, quiet: bool = False):
    """Download a file with optional progress reporting."""
    if not quiet:
        import sys
        print(f"  Downloading {url.split('/')[-1]}...", end="", flush=True)

    def _reporthook(count, block_size, total_size):
        if not quiet and total_size > 0:
            pct = min(100, count * block_size * 100 // total_size)
            print(f"\r  Downloading {url.split('/')[-1]}... {pct}%", end="", flush=True)

    try:
        urlretrieve(url, dest, reporthook=_reporthook)
    except URLError as e:
        raise RuntimeError(f"Failed to download {url}: {e}") from e

    if not quiet:
        print(" done.")


def fetch_geo(
    accession: str,
    output_dir: Union[str, Path] = ".",
    *,
    quiet: bool = False,
) -> Path:
    """Download and extract RCC files from a GEO accession.

    Looks for the ``GSE*_RAW.tar`` supplement file, downloads it,
    and extracts any ``.RCC`` (or ``.RCC.gz``) files into *output_dir*.

    Parameters
    ----------
    accession : str
        GEO series accession (e.g. ``"GSE275334"``).
    output_dir : str or Path
        Directory to extract RCC files into.  A subdirectory named
        after the accession will be created.
    quiet : bool
        Suppress progress output.

    Returns
    -------
    Path
        Path to the directory containing extracted RCC files.
    """
    accession = accession.upper().strip()
    if not accession.startswith("GSE"):
        raise ValueError(f"Expected a GSE accession, got: {accession}")

    output_dir = Path(output_dir) / accession
    output_dir.mkdir(parents=True, exist_ok=True)

    # Try to download the _RAW.tar supplement
    base_url = _gse_to_ftp_dir(accession)
    tar_name = f"{accession}_RAW.tar"
    tar_url = base_url + tar_name

    with tempfile.TemporaryDirectory() as tmpdir:
        tar_path = Path(tmpdir) / tar_name
        _download_with_progress(tar_url, tar_path, quiet=quiet)

        if not quiet:
            print(f"  Extracting RCC files to {output_dir}/")

        n_extracted = 0
        with tarfile.open(tar_path) as tar:
            for member in tar.getmembers():
                name_lower = member.name.lower()
                if name_lower.endswith(".rcc") or name_lower.endswith(".rcc.gz"):
                    # Extract to a flat directory (no subdirs)
                    member_name = Path(member.name).name
                    f = tar.extractfile(member)
                    if f is None:
                        continue

                    if member_name.lower().endswith(".rcc.gz"):
                        # Decompress gzipped RCC
                        out_name = member_name[:-3]  # strip .gz
                        with gzip.open(f) as gz:
                            content = gz.read()
                        (output_dir / out_name).write_bytes(content)
                    else:
                        (output_dir / member_name).write_bytes(f.read())

                    n_extracted += 1

    if n_extracted == 0:
        raise RuntimeError(
            f"No RCC files found in {tar_name}. The archive may contain "
            "a different format — check the GEO page manually."
        )

    if not quiet:
        print(f"  Extracted {n_extracted} RCC files.")

    return output_dir
