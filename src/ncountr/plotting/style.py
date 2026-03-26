"""Shared matplotlib style settings."""

from __future__ import annotations

import matplotlib.pyplot as plt

_NCOUNTR_STYLE = {
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}


def set_style(dpi: int | None = None) -> None:
    """Apply ncountr matplotlib style.

    Parameters
    ----------
    dpi : int, optional
        Override default DPI (150 display / 300 save).
    """
    plt.rcParams.update(_NCOUNTR_STYLE)
    if dpi is not None:
        plt.rcParams["figure.dpi"] = dpi
        plt.rcParams["savefig.dpi"] = dpi
