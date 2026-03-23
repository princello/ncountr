"""Parse Nanostring RCC files into a NanostringExperiment."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Union

import pandas as pd

from ncountr.experiment import NanostringExperiment


def parse_rcc(filepath: Union[str, Path]) -> dict:
    """Parse a single Nanostring RCC file.

    Parameters
    ----------
    filepath : str or Path
        Path to a ``.RCC`` file.

    Returns
    -------
    dict
        Keys: ``sample`` (sample attributes), ``lane`` (lane attributes),
        ``counts`` (dict of ``(CodeClass, GeneName) -> count``).
    """
    result: dict = {"sample": {}, "lane": {}, "counts": {}}
    current_section: str | None = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith("<") and not line.startswith("</"):
                current_section = line.strip("<>").lower()
                continue
            if line.startswith("</") or not line:
                current_section = None
                continue
            if current_section == "sample_attributes":
                parts = line.split(",", 1)
                if len(parts) == 2:
                    result["sample"][parts[0]] = parts[1]
            elif current_section == "lane_attributes":
                parts = line.split(",", 1)
                if len(parts) == 2:
                    result["lane"][parts[0]] = parts[1]
            elif current_section == "code_summary":
                parts = line.split(",")
                if len(parts) >= 4 and parts[0] != "CodeClass":
                    code_class, name, _accession, count = (
                        parts[0],
                        parts[1],
                        parts[2],
                        int(parts[3]),
                    )
                    result["counts"][(code_class, name)] = count
    return result


def read_rcc(
    rcc_dirs: Union[str, Path, list[Union[str, Path]]],
    *,
    file_pattern: str = "*.RCC",
    sample_id_pattern: str = r"(\d+)",
    sample_id_field: str = "ID",
    sample_id_from: str = "field",
    sample_meta: dict[str, dict] | None = None,
) -> NanostringExperiment:
    """Read RCC files from one or more directories into a NanostringExperiment.

    Parameters
    ----------
    rcc_dirs : str, Path, or list thereof
        Directory or directories containing ``.RCC`` files.
    file_pattern : str
        Glob pattern to match RCC files within each directory.
    sample_id_pattern : str
        Regex applied to extract a clean sample ID.
        The first capture group is used.
    sample_id_field : str
        Which field in the ``<Sample_Attributes>`` section holds the sample ID.
        Only used when ``sample_id_from="field"``.
    sample_id_from : str
        Where to extract the sample ID from.  ``"field"`` (default) uses the
        ``sample_id_field`` from the RCC file.  ``"filename"`` applies the
        regex to the filename instead — useful when internal IDs are
        inconsistent across files.
    sample_meta : dict[str, dict] | None
        Optional per-sample metadata, keyed by sample ID.

    Returns
    -------
    NanostringExperiment
    """
    if isinstance(rcc_dirs, (str, Path)):
        rcc_dirs = [rcc_dirs]

    # Collect RCC file paths
    rcc_files: list[Path] = []
    for d in rcc_dirs:
        d = Path(d)
        rcc_files.extend(sorted(d.glob(file_pattern)))

    if not rcc_files:
        raise FileNotFoundError(
            f"No files matching '{file_pattern}' found in {rcc_dirs}"
        )

    all_data: dict[str, dict] = {}
    lane_rows: dict[str, dict] = {}

    for fp in rcc_files:
        parsed = parse_rcc(fp)
        if sample_id_from == "filename":
            raw_id = fp.stem  # filename without extension
        else:
            raw_id = parsed["sample"].get(sample_id_field, "")
        match = re.search(sample_id_pattern, raw_id)
        if not match:
            continue
        sid = match.group(1)

        all_data[sid] = parsed["counts"]
        lane_rows[sid] = {
            "FovCount": int(parsed["lane"].get("FovCount", 0)),
            "FovCounted": int(parsed["lane"].get("FovCounted", 0)),
            "BindingDensity": float(parsed["lane"].get("BindingDensity", 0)),
            "CartridgeID": parsed["lane"].get("CartridgeID", ""),
        }

    if not all_data:
        raise ValueError("No samples could be parsed from the RCC files")

    # Collect all code class / gene pairs
    all_keys: set[tuple[str, str]] = set()
    for counts in all_data.values():
        all_keys.update(counts.keys())

    endogenous = sorted({name for cls, name in all_keys if cls == "Endogenous"})
    positive = sorted({name for cls, name in all_keys if cls == "Positive"})
    negative = sorted({name for cls, name in all_keys if cls == "Negative"})
    housekeeping = sorted({name for cls, name in all_keys if cls == "Housekeeping"})

    samples = sorted(all_data, key=lambda x: (0, int(x)) if x.isdigit() else (1, x))

    def _build_df(genes: list[str], code_class: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                sid: {g: all_data[sid].get((code_class, g), 0) for g in genes}
                for sid in samples
            }
        )

    raw_counts = _build_df(endogenous, "Endogenous")
    pos_counts = _build_df(positive, "Positive")
    neg_counts = _build_df(negative, "Negative")
    hk_counts = _build_df(housekeeping, "Housekeeping")

    lane_df = pd.DataFrame(lane_rows).T
    lane_df.index.name = "sample"
    # Ensure proper dtypes for numeric columns
    for col in ("FovCount", "FovCounted"):
        if col in lane_df.columns:
            lane_df[col] = pd.to_numeric(lane_df[col], errors="coerce").astype(int)
    if "BindingDensity" in lane_df.columns:
        lane_df["BindingDensity"] = pd.to_numeric(lane_df["BindingDensity"], errors="coerce").astype(float)

    # Build sample_meta DataFrame
    if sample_meta:
        meta_df = pd.DataFrame.from_dict(sample_meta, orient="index")
        meta_df.index.name = "sample"
    else:
        meta_df = pd.DataFrame(index=pd.Index(samples, name="sample"))

    return NanostringExperiment(
        raw_counts=raw_counts,
        pos_counts=pos_counts,
        neg_counts=neg_counts,
        hk_counts=hk_counts,
        sample_meta=meta_df,
        lane_info=lane_df,
    )
