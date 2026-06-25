"""Tests for ncountr.genesets — GMT parsing, Hallmark sets, filtering."""

from __future__ import annotations

from pathlib import Path

import pytest

from ncountr.genesets import (
    load_gmt,
    save_gmt,
    get_hallmark_set,
    list_hallmark_sets,
    get_all_gene_sets,
    filter_gene_sets,
)


# ---------------------------------------------------------------------------
# GMT parsing
# ---------------------------------------------------------------------------


class TestGmt:
    def test_load_gmt(self, tmp_path: Path):
        gmt = tmp_path / "sets.gmt"
        gmt.write_text(
            "SET_A\tdescription A\tGENE1\tGENE2\tGENE3\n"
            "SET_B\thttp://example.org\tGENE4\tGENE5\n"
        )
        sets = load_gmt(gmt)
        assert sets["SET_A"] == ["GENE1", "GENE2", "GENE3"]
        assert sets["SET_B"] == ["GENE4", "GENE5"]

    def test_load_gmt_skips_short_lines(self, tmp_path: Path):
        gmt = tmp_path / "sets.gmt"
        gmt.write_text("SET_A\tdesc\tGENE1\nBADLINE\n\n")
        sets = load_gmt(gmt)
        assert list(sets.keys()) == ["SET_A"]

    def test_save_then_load_roundtrip(self, tmp_path: Path):
        original = {"S1": ["A", "B", "C"], "S2": ["D", "E"]}
        gmt = tmp_path / "out.gmt"
        save_gmt(original, gmt)
        loaded = load_gmt(gmt)
        assert loaded == original

    def test_save_with_descriptions(self, tmp_path: Path):
        gmt = tmp_path / "out.gmt"
        save_gmt({"S1": ["A", "B"]}, gmt, descriptions={"S1": "my desc"})
        text = gmt.read_text()
        assert "my desc" in text


# ---------------------------------------------------------------------------
# Hallmark sets
# ---------------------------------------------------------------------------


class TestHallmark:
    def test_get_with_prefix(self):
        genes = get_hallmark_set("HALLMARK_INTERFERON_ALPHA_RESPONSE")
        assert "MX1" in genes

    def test_get_without_prefix(self):
        genes = get_hallmark_set("INTERFERON_ALPHA_RESPONSE")
        assert "MX1" in genes

    def test_case_insensitive(self):
        genes = get_hallmark_set("interferon_alpha_response")
        assert "MX1" in genes

    def test_unknown_raises_keyerror(self):
        with pytest.raises(KeyError, match="Unknown Hallmark set"):
            get_hallmark_set("NOT_A_REAL_SET")

    def test_list_returns_known_sets(self):
        names = list_hallmark_sets()
        assert "HALLMARK_INTERFERON_GAMMA_RESPONSE" in names
        assert all(n.startswith("HALLMARK_") for n in names)

    def test_returned_list_is_a_copy(self):
        genes = get_hallmark_set("COMPLEMENT")
        genes.append("FAKE")
        assert "FAKE" not in get_hallmark_set("COMPLEMENT")


# ---------------------------------------------------------------------------
# get_all_gene_sets / filter_gene_sets
# ---------------------------------------------------------------------------


class TestGetAllGeneSets:
    def test_merges_custom_and_hallmark(self):
        merged = get_all_gene_sets()
        assert "IFN_JAKSTAT" in merged          # custom
        assert "HALLMARK_COMPLEMENT" in merged   # hallmark

    def test_values_are_lists(self):
        merged = get_all_gene_sets()
        assert all(isinstance(v, list) for v in merged.values())


class TestFilterGeneSets:
    def test_keeps_sufficient_overlap(self):
        measured = ["MX1", "IFIT1", "ISG15", "OAS1", "STAT1", "CXCL10"]
        sets = {"IFN": ["MX1", "IFIT1", "ISG15", "OAS1", "STAT1", "ZZZ"]}
        filtered = filter_gene_sets(sets, measured, min_overlap=5)
        assert "IFN" in filtered
        # Only measured genes are retained
        assert "ZZZ" not in filtered["IFN"]
        assert set(filtered["IFN"]) == {"MX1", "IFIT1", "ISG15", "OAS1", "STAT1"}

    def test_drops_insufficient_overlap(self):
        measured = ["MX1", "IFIT1"]
        sets = {"IFN": ["MX1", "IFIT1", "ISG15", "OAS1", "STAT1"]}
        filtered = filter_gene_sets(sets, measured, min_overlap=5)
        assert filtered == {}

    def test_drops_oversized_overlap(self):
        measured = [f"G{i}" for i in range(20)]
        sets = {"big": [f"G{i}" for i in range(20)]}
        filtered = filter_gene_sets(sets, measured, min_overlap=5, max_overlap=10)
        assert filtered == {}
