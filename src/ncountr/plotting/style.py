"""Shared matplotlib style settings."""

from __future__ import annotations

import matplotlib.pyplot as plt

_NCOUNTR_STYLE = {
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
}


def set_style(dpi: int | None = None) -> None:
    """Apply ncountr matplotlib style.

    Parameters
    ----------
    dpi : int, optional
        Override default DPI (150 display / 200 save).
    """
    plt.rcParams.update(_NCOUNTR_STYLE)
    if dpi is not None:
        plt.rcParams["figure.dpi"] = dpi
        plt.rcParams["savefig.dpi"] = dpi
