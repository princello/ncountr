"""Built-in gene sets and cell type markers."""

from __future__ import annotations

_GENE_SETS: dict[str, list[str]] = {
    "IFN_JAKSTAT": [
        "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
        "ISG15", "ISG20", "IFI44", "IFI44L", "IFI35", "IFI6", "IFI16", "IFI27",
        "IFITM1", "IFITM2", "IFITM3",
        "STAT1", "STAT2", "STAT3", "STAT4",
        "IRF1", "IRF3", "IRF7", "IRF9",
        "OAS1", "OAS2", "OAS3", "OASL",
        "RSAD2", "BST2", "DDX58", "HERC5", "USP18",
        "CXCL10", "CXCL11", "CXCL9",
        "IFNG", "IFNB1", "IFNL1",
        "PARP9", "PARP14", "TRIM21", "TRIM25",
        "GBP1", "GBP2", "GBP4", "GBP5",
        "JAK1", "JAK2", "TYK2",
        "PSMB8", "PSMB9", "PSMB10",
    ],
}

_CELL_MARKERS: dict[str, list[str]] = {
    "T_cells": ["CD3D", "CD3E", "CD3G"],
    "CD4_T": ["CD4"],
    "CD8_T": ["CD8A", "CD8B"],
    "Monocytes": ["CD14", "FCGR3A", "CD68"],
    "B_cells": ["MS4A1", "CD19", "CD79A"],
    "NK_cells": ["GNLY", "NKG7", "KLRD1"],
}


def get_gene_set(name: str) -> list[str]:
    """Return a built-in gene set by name.

    Parameters
    ----------
    name : str
        Case-insensitive gene set name.  Use :func:`list_gene_sets` to
        see available names.

    Raises
    ------
    KeyError
        If the name is not found.
    """
    key = name.upper()
    for k, v in _GENE_SETS.items():
        if k.upper() == key:
            return list(v)
    raise KeyError(
        f"Unknown gene set {name!r}. Available: {list(_GENE_SETS.keys())}"
    )


def list_gene_sets() -> list[str]:
    """Return names of all built-in gene sets."""
    return list(_GENE_SETS.keys())


def get_cell_markers(cell_type: str | None = None) -> dict[str, list[str]] | list[str]:
    """Return cell type marker genes.

    Parameters
    ----------
    cell_type : str, optional
        If given, return markers for that cell type.
        If None, return the full dictionary.
    """
    if cell_type is None:
        return dict(_CELL_MARKERS)
    key = cell_type.replace(" ", "_")
    for k, v in _CELL_MARKERS.items():
        if k.lower() == key.lower():
            return list(v)
    raise KeyError(
        f"Unknown cell type {cell_type!r}. Available: {list(_CELL_MARKERS.keys())}"
    )
