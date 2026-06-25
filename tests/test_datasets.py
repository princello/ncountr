"""Tests for ncountr.datasets — built-in gene sets and cell markers."""

from __future__ import annotations

import pytest

from ncountr.datasets import get_gene_set, list_gene_sets, get_cell_markers


# ---------------------------------------------------------------------------
# get_gene_set
# ---------------------------------------------------------------------------


class TestGetGeneSet:
    def test_custom_set(self):
        genes = get_gene_set("IFN_JAKSTAT")
        assert "MX1" in genes
        assert isinstance(genes, list)

    def test_case_insensitive(self):
        assert get_gene_set("ifn_jakstat") == get_gene_set("IFN_JAKSTAT")

    def test_returned_list_is_a_copy(self):
        genes = get_gene_set("IFN_JAKSTAT")
        genes.append("FAKE")
        assert "FAKE" not in get_gene_set("IFN_JAKSTAT")

    def test_falls_through_to_hallmark(self):
        # Not a custom set, but a Hallmark set (no prefix)
        genes = get_gene_set("INTERFERON_ALPHA_RESPONSE")
        assert "MX1" in genes

    def test_unknown_raises_keyerror(self):
        with pytest.raises(KeyError, match="Unknown gene set"):
            get_gene_set("DEFINITELY_NOT_A_SET")


# ---------------------------------------------------------------------------
# list_gene_sets
# ---------------------------------------------------------------------------


class TestListGeneSets:
    def test_includes_custom_and_hallmark(self):
        names = list_gene_sets()
        assert "IFN_JAKSTAT" in names
        assert any(n.startswith("HALLMARK_") for n in names)


# ---------------------------------------------------------------------------
# get_cell_markers
# ---------------------------------------------------------------------------


class TestGetCellMarkers:
    def test_all_markers_dict(self):
        markers = get_cell_markers()
        assert isinstance(markers, dict)
        assert "T_cells" in markers
        assert markers["T_cells"] == ["CD3D", "CD3E", "CD3G"]

    def test_single_cell_type(self):
        markers = get_cell_markers("CD8_T")
        assert markers == ["CD8A", "CD8B"]

    def test_space_in_name(self):
        # "B cells" should map to the "B_cells" key
        markers = get_cell_markers("B cells")
        assert "MS4A1" in markers

    def test_case_insensitive(self):
        assert get_cell_markers("nk_cells") == get_cell_markers("NK_cells")

    def test_unknown_raises_keyerror(self):
        with pytest.raises(KeyError, match="Unknown cell type"):
            get_cell_markers("Dragons")
