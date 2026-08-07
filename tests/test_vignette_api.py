"""Static checks on the vignette notebooks.

These exist because the package test suite can be entirely green while an
advertised tutorial calls a function signature that no longer exists — which is
exactly what happened to four of the five vignettes: they called
``plot_pathway_scores(df, score_column=..., group_column=...)`` and
``plot_heatmap(exp, genes=..., samples=..., z_score=...)``, neither of which the
package has ever accepted.

The checks here are static (ast + inspect.signature), so they are fast enough to
run on every push. Actually executing the notebooks needs network access and the
GEO downloads; that runs separately in .github/workflows/vignettes.yml.
"""

from __future__ import annotations

import ast
import inspect
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
VIGNETTE_DIR = REPO_ROOT / "examples" / "vignettes"


def _notebook_paths() -> list[Path]:
    return sorted(VIGNETTE_DIR.glob("*.ipynb"))


def _notebook_id(path: Path) -> str:
    return path.name


NOTEBOOKS = _notebook_paths()


def _public_signatures() -> dict[str, inspect.Signature]:
    """Map public ncountr callable name -> live signature."""
    import ncountr
    from ncountr.core import de, gsea, normalize, pathway, qc
    from ncountr.io import export, geo, rcc
    from ncountr.plotting import (
        correlation_plots,
        de_plots,
        gsea_plots,
        heatmaps,
        pathway_plots,
        qc_plots,
    )

    modules = [
        ncountr, de, gsea, normalize, pathway, qc, export, geo, rcc,
        correlation_plots, de_plots, gsea_plots, heatmaps, pathway_plots, qc_plots,
    ]

    sigs: dict[str, inspect.Signature] = {}
    for module in modules:
        for name, obj in vars(module).items():
            if name.startswith("_") or not callable(obj):
                continue
            # Only functions defined inside ncountr — not re-exported third-party
            # callables such as pd.DataFrame that happen to be module globals.
            if not getattr(obj, "__module__", "").startswith("ncountr"):
                continue
            try:
                sigs.setdefault(name, inspect.signature(obj))
            except (TypeError, ValueError):
                continue
    return sigs


SIGNATURES = _public_signatures()


def _code_cells(nb: dict):
    """Yield (cell_index, source) for code cells, with notebook magics stripped."""
    for idx, cell in enumerate(nb["cells"]):
        if cell.get("cell_type") != "code":
            continue
        lines = "".join(cell.get("source", [])).splitlines()
        # %pip / !git etc. are not valid Python; blank them but keep line numbers
        cleaned = "\n".join(
            "" if line.lstrip().startswith(("%", "!", "?")) else line for line in lines
        )
        yield idx, cleaned


def _called_name(node: ast.Call) -> str | None:
    func = node.func
    if isinstance(func, ast.Attribute):
        return func.attr
    if isinstance(func, ast.Name):
        return func.id
    return None


def test_signature_map_is_populated():
    """Guard against the audit silently passing because it found no functions."""
    assert len(SIGNATURES) > 20, f"only discovered {len(SIGNATURES)} ncountr callables"
    for expected in ("read_rcc", "de", "normalize", "to_anndata",
                     "plot_heatmap", "plot_pathway_scores"):
        assert expected in SIGNATURES, f"{expected} missing from signature map"


def test_vignettes_exist():
    assert NOTEBOOKS, f"no notebooks found under {VIGNETTE_DIR}"


@pytest.mark.parametrize("path", NOTEBOOKS, ids=_notebook_id)
def test_notebook_is_parseable(path: Path):
    """Every code cell must be syntactically valid Python."""
    nb = json.loads(path.read_text())
    for idx, source in _code_cells(nb):
        try:
            ast.parse(source)
        except SyntaxError as exc:  # pragma: no cover - failure path
            pytest.fail(f"{path.name} cell {idx}: syntax error: {exc}")


@pytest.mark.parametrize("path", NOTEBOOKS, ids=_notebook_id)
def test_notebook_calls_match_signatures(path: Path):
    """Keyword arguments passed to ncountr functions must actually exist.

    Only names that unambiguously resolve to an ncountr callable are checked, so
    a same-named method on a third-party object does not produce a false
    positive beyond the shared-name case.
    """
    nb = json.loads(path.read_text())
    problems: list[str] = []

    for idx, source in _code_cells(nb):
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue  # reported by test_notebook_is_parseable

        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = _called_name(node)
            if name is None or name not in SIGNATURES:
                continue

            sig = SIGNATURES[name]
            params = sig.parameters
            accepts_var_kw = any(p.kind is p.VAR_KEYWORD for p in params.values())
            accepts_var_pos = any(p.kind is p.VAR_POSITIONAL for p in params.values())

            for kw in node.keywords:
                if kw.arg is None:  # **kwargs splat — cannot check statically
                    continue
                if kw.arg not in params and not accepts_var_kw:
                    problems.append(
                        f"cell {idx}: {name}() got unexpected keyword "
                        f"{kw.arg!r} — real signature: {name}{sig}"
                    )

            positional = [
                p for p in params.values()
                if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
            ]
            if not accepts_var_pos and len(node.args) > len(positional):
                problems.append(
                    f"cell {idx}: {name}() given {len(node.args)} positional "
                    f"args but accepts {len(positional)} — "
                    f"real signature: {name}{sig}"
                )

    if problems:
        pytest.fail(
            f"{path.name}: {len(problems)} API mismatch(es):\n  "
            + "\n  ".join(problems)
        )


@pytest.mark.parametrize("path", NOTEBOOKS, ids=_notebook_id)
def test_notebook_has_no_committed_output(path: Path):
    """Source notebooks must be stored clean: no outputs, no execution state."""
    nb = json.loads(path.read_text())
    dirty: list[str] = []

    for idx, cell in enumerate(nb["cells"]):
        if cell.get("cell_type") != "code":
            continue
        outputs = cell.get("outputs") or []
        if outputs:
            kinds = sorted({o.get("output_type", "?") for o in outputs})
            dirty.append(f"cell {idx}: {len(outputs)} output(s) [{', '.join(kinds)}]")
        if cell.get("execution_count") is not None:
            dirty.append(f"cell {idx}: execution_count={cell['execution_count']}")

    if dirty:
        pytest.fail(
            f"{path.name} has committed execution state "
            f"(strip it before committing):\n  " + "\n  ".join(dirty)
        )


@pytest.mark.parametrize("path", NOTEBOOKS, ids=_notebook_id)
def test_notebook_has_no_traceback(path: Path):
    """A committed traceback means a broken tutorial was published."""
    nb = json.loads(path.read_text())
    for idx, cell in enumerate(nb["cells"]):
        for output in cell.get("outputs") or []:
            if output.get("output_type") == "error":
                pytest.fail(
                    f"{path.name} cell {idx} contains an error output: "
                    f"{output.get('ename')}: {output.get('evalue')}"
                )
