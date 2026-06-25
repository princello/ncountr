"""Tests for ncountr.io.geo — GEO download (network mocked)."""

from __future__ import annotations

import gzip
import io
import shutil
import tarfile
from pathlib import Path

import pytest

from ncountr.io import geo
from ncountr.io.geo import _gse_to_ftp_dir, fetch_geo


_RCC_BYTES = b"<Header>\nFileVersion,1.7\n</Header>\n"


def _build_raw_tar(path: Path, *, include_rcc: bool = True) -> None:
    """Create a GSE*_RAW.tar with a plain RCC, a gzipped RCC, and junk."""
    with tarfile.open(path, "w") as tar:
        if include_rcc:
            info = tarfile.TarInfo(name="GSM1_lung.RCC")
            info.size = len(_RCC_BYTES)
            tar.addfile(info, io.BytesIO(_RCC_BYTES))

            gz = gzip.compress(_RCC_BYTES)
            info_gz = tarfile.TarInfo(name="nested/GSM2_lung.RCC.gz")
            info_gz.size = len(gz)
            tar.addfile(info_gz, io.BytesIO(gz))

        junk = b"not an rcc file"
        info_junk = tarfile.TarInfo(name="filelist.txt")
        info_junk.size = len(junk)
        tar.addfile(info_junk, io.BytesIO(junk))


@pytest.fixture()
def mock_download(monkeypatch, tmp_path: Path):
    """Patch urlretrieve to copy a prebuilt tar instead of hitting the network."""
    def _factory(include_rcc: bool = True):
        prebuilt = tmp_path / "prebuilt_RAW.tar"
        _build_raw_tar(prebuilt, include_rcc=include_rcc)

        def fake_urlretrieve(url, dest, reporthook=None):
            shutil.copy(prebuilt, dest)
            return dest, None

        monkeypatch.setattr(geo, "urlretrieve", fake_urlretrieve)

    return _factory


# ---------------------------------------------------------------------------
# _gse_to_ftp_dir
# ---------------------------------------------------------------------------


class TestGseToFtpDir:
    def test_standard_accession(self):
        url = _gse_to_ftp_dir("GSE275334")
        assert url == (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275334/suppl/"
        )

    def test_lowercased_input_uppercased(self):
        assert _gse_to_ftp_dir("gse275334") == _gse_to_ftp_dir("GSE275334")


# ---------------------------------------------------------------------------
# fetch_geo
# ---------------------------------------------------------------------------


class TestFetchGeo:
    def test_extracts_rcc_files(self, mock_download, tmp_path: Path):
        mock_download(include_rcc=True)
        out = fetch_geo("GSE999999", output_dir=tmp_path, quiet=True)
        assert out == tmp_path / "GSE999999"
        rcc_files = sorted(p.name for p in out.glob("*.RCC"))
        assert rcc_files == ["GSM1_lung.RCC", "GSM2_lung.RCC"]

    def test_gz_member_is_decompressed(self, mock_download, tmp_path: Path):
        mock_download(include_rcc=True)
        out = fetch_geo("GSE999999", output_dir=tmp_path, quiet=True)
        # The .RCC.gz extracted to .RCC with decompressed content
        assert (out / "GSM2_lung.RCC").read_bytes() == _RCC_BYTES

    def test_flattens_nested_paths(self, mock_download, tmp_path: Path):
        mock_download(include_rcc=True)
        out = fetch_geo("GSE999999", output_dir=tmp_path, quiet=True)
        assert not (out / "nested").exists()

    def test_non_gse_accession_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="Expected a GSE accession"):
            fetch_geo("PRJNA12345", output_dir=tmp_path, quiet=True)

    def test_no_rcc_in_archive_raises(self, mock_download, tmp_path: Path):
        mock_download(include_rcc=False)
        with pytest.raises(RuntimeError, match="No RCC files found"):
            fetch_geo("GSE999999", output_dir=tmp_path, quiet=True)
