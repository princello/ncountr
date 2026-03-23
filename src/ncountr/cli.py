"""Command-line interface for ncountr."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ncountr import __version__


@click.group()
@click.version_option(__version__)
def cli():
    """ncountr — Nanostring nCounter analysis pipeline."""


@cli.command()
def init():
    """Print a template YAML config to stdout."""
    from ncountr.config import generate_template
    click.echo(generate_template())


@cli.command()
@click.argument("config_path", type=click.Path(exists=True))
@click.option("--output-dir", "-o", default=None, help="Override output directory.")
def run(config_path: str, output_dir: str | None):
    """Run the full analysis pipeline from a config file."""
    import matplotlib
    matplotlib.use("Agg")

    from ncountr.config import load_config
    from ncountr.io.rcc import read_rcc
    from ncountr.io.export import export_counts, export_qc, export_de
    from ncountr.core.qc import qc
    from ncountr.core.normalize import normalize
    from ncountr.core.de import de
    from ncountr.core.pathway import score_gene_set
    from ncountr.plotting.qc_plots import plot_qc
    from ncountr.plotting.de_plots import plot_volcano
    import matplotlib.pyplot as plt

    cfg = load_config(config_path)
    outdir = Path(output_dir or cfg.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    click.echo(f"ncountr v{__version__} — running pipeline")
    click.echo(f"Config: {config_path}")
    click.echo(f"Output: {outdir}")

    # Step 1: Parse
    if not cfg.rcc_dirs:
        click.echo("Error: no rcc_dirs specified in config.", err=True)
        sys.exit(1)

    click.echo("\n[1/5] Parsing RCC files...")
    experiment = read_rcc(
        cfg.rcc_dirs,
        file_pattern=cfg.file_pattern,
        sample_id_pattern=cfg.sample_id_pattern,
        sample_meta=cfg.sample_meta or None,
    )
    click.echo(f"  {experiment}")

    # Step 2: QC
    click.echo("\n[2/5] Running QC...")
    qc_results = qc(
        experiment,
        fov_ratio_threshold=cfg.fov_ratio_threshold,
        pos_r2_threshold=cfg.pos_r2_threshold,
    )
    click.echo(f"  FOV pass: {qc_results['FovPass'].sum()}/{len(qc_results)}")
    click.echo(f"  PosCtrl pass: {qc_results['PosCtrlPass'].sum()}/{len(qc_results)}")

    fig = plot_qc(experiment, output=outdir / f"qc_summary.{cfg.figure_format}",
                  fov_threshold=cfg.fov_ratio_threshold)
    plt.close(fig)
    click.echo(f"  Saved: qc_summary.{cfg.figure_format}")

    # Step 3: Normalize
    click.echo("\n[3/5] Normalizing...")
    normalize(experiment, method=cfg.normalization_method)
    click.echo(f"  Method: {cfg.normalization_method}")

    # Step 4: DE (if comparison groups are defined)
    if len(cfg.comparison) == 2 and cfg.sample_meta:
        click.echo("\n[4/5] Differential expression...")
        group_a_name, group_b_name = cfg.comparison
        group_col = cfg.group_column

        group_a = [s for s, m in cfg.sample_meta.items()
                   if m.get(group_col) == group_a_name and s not in cfg.exclude_samples]
        group_b = [s for s, m in cfg.sample_meta.items()
                   if m.get(group_col) == group_b_name and s not in cfg.exclude_samples]

        if group_a and group_b:
            de_results = de(
                experiment,
                group_a=group_a,
                group_b=group_b,
                test=cfg.de_test,
                correction=cfg.de_correction,
            )
            n_sig = (de_results["padj"] < 0.05).sum()
            click.echo(f"  {group_a_name} ({len(group_a)}) vs {group_b_name} ({len(group_b)})")
            click.echo(f"  Significant (padj < 0.05): {n_sig}")

            # Resolve highlight genes
            hl_genes = None
            hl_label = cfg.highlight_label
            if cfg.highlight_genes:
                if isinstance(cfg.highlight_genes, str):
                    # Treat as builtin gene set name
                    from ncountr.datasets import get_gene_set
                    try:
                        hl_genes = get_gene_set(cfg.highlight_genes)
                        if hl_label == "Highlighted":
                            hl_label = cfg.highlight_genes.replace("_", "/")
                        click.echo(f"  Highlighting: {hl_label} ({len(hl_genes)} genes)")
                    except KeyError:
                        click.echo(f"  Warning: gene set '{cfg.highlight_genes}' not found, skipping highlight")
                elif isinstance(cfg.highlight_genes, list):
                    hl_genes = cfg.highlight_genes
                    click.echo(f"  Highlighting: {hl_label} ({len(hl_genes)} genes)")

            fig = plot_volcano(
                de_results,
                highlight_genes=hl_genes,
                highlight_label=hl_label,
                highlight_color=cfg.highlight_color,
                output=outdir / f"volcano.{cfg.figure_format}",
                title=f"Nanostring DE: {group_a_name} vs {group_b_name}\n"
                      f"({len(group_a)} {group_a_name.lower()}, {len(group_b)} {group_b_name.lower()})",
            )
            plt.close(fig)
            click.echo(f"  Saved: volcano.{cfg.figure_format}")
        else:
            click.echo(f"  Skipped: no samples in comparison groups")
    else:
        click.echo("\n[4/5] DE: skipped (no comparison defined)")

    # Step 5: Gene set scoring
    click.echo("\n[5/5] Gene set scoring...")
    for gs_name, gs_value in cfg.gene_sets.items():
        try:
            if gs_value == "builtin":
                scores = score_gene_set(experiment, gene_set=gs_name)
            elif isinstance(gs_value, list):
                scores = score_gene_set(experiment, gene_set=gs_value)
            else:
                continue
            click.echo(f"  {gs_name}: scored {len(scores)} samples")
            scores.to_csv(outdir / f"pathway_score_{gs_name}.csv")
        except (KeyError, ValueError) as e:
            click.echo(f"  {gs_name}: {e}")

    # Export
    export_counts(experiment, outdir)
    export_qc(experiment, outdir)
    export_de(experiment, outdir)

    click.echo(f"\nDone! Results in {outdir}")


@cli.command()
@click.option("--rcc-dir", "-d", required=True, type=click.Path(exists=True),
              help="Directory containing RCC files.")
@click.option("--pattern", "-p", default="*.RCC", help="File glob pattern.")
@click.option("--id-pattern", default=r"(\d+)", help="Regex for sample ID extraction.")
@click.option("--output", "-o", default="raw_counts.csv", help="Output CSV path.")
def parse(rcc_dir: str, pattern: str, id_pattern: str, output: str):
    """Parse RCC files into a raw count matrix."""
    from ncountr.io.rcc import read_rcc

    experiment = read_rcc(rcc_dir, file_pattern=pattern, sample_id_pattern=id_pattern)
    import pandas as pd
    combined = pd.concat([experiment.raw_counts, experiment.hk_counts])
    combined.to_csv(output)
    click.echo(f"Parsed {experiment.n_genes} genes x {experiment.n_samples} samples → {output}")


@cli.command(name="qc")
@click.option("--counts", "-c", required=True, type=click.Path(exists=True),
              help="Raw count CSV (genes x samples).")
@click.option("--output", "-o", default="qc_results.csv", help="Output CSV path.")
def qc_cmd(counts: str, output: str):
    """Run QC checks on a count matrix."""
    import pandas as pd
    from ncountr.experiment import NanostringExperiment
    from ncountr.core.qc import qc as run_qc

    df = pd.read_csv(counts, index_col=0)
    exp = NanostringExperiment(
        raw_counts=df,
        pos_counts=pd.DataFrame(),
        neg_counts=pd.DataFrame(),
        hk_counts=pd.DataFrame(),
    )
    click.echo("Note: QC from CSV has limited checks (no pos/neg controls).")
    click.echo(f"  {exp.n_genes} genes x {exp.n_samples} samples")


@cli.command(name="fetch-geo")
@click.argument("accession")
@click.option("--output-dir", "-o", default=".", help="Parent directory for downloaded files.")
def fetch_geo_cmd(accession: str, output_dir: str):
    """Download RCC files from GEO (e.g. ncountr fetch-geo GSE275334)."""
    from ncountr.io.geo import fetch_geo

    rcc_dir = fetch_geo(accession, output_dir=output_dir)
    n_files = len(list(rcc_dir.glob("*.RCC")))
    click.echo(f"Downloaded {n_files} RCC files → {rcc_dir}")


@cli.command()
@click.option("--counts", "-c", required=True, type=click.Path(exists=True),
              help="Normalized count CSV (genes x samples).")
@click.option("--groups", "-g", required=True, multiple=True,
              help="Group definition: name:S1,S2,S3")
@click.option("--test", "-t", default="mannwhitneyu", help="Statistical test.")
@click.option("--output", "-o", default="de_results.csv", help="Output CSV path.")
def de_cmd(counts: str, groups: tuple[str], test: str, output: str):
    """Run differential expression from a count CSV."""
    import pandas as pd
    from ncountr.experiment import NanostringExperiment
    from ncountr.core.de import de as run_de

    df = pd.read_csv(counts, index_col=0)

    parsed_groups: dict[str, list[str]] = {}
    for g in groups:
        name, samples_str = g.split(":", 1)
        parsed_groups[name] = samples_str.split(",")

    group_names = list(parsed_groups.keys())
    if len(group_names) != 2:
        click.echo("Error: exactly 2 groups required.", err=True)
        sys.exit(1)

    exp = NanostringExperiment(
        raw_counts=df,
        pos_counts=pd.DataFrame(),
        neg_counts=pd.DataFrame(),
        hk_counts=pd.DataFrame(),
    )

    de_results = run_de(
        exp,
        group_a=parsed_groups[group_names[0]],
        group_b=parsed_groups[group_names[1]],
        counts=df,
        test=test,
    )
    de_results.to_csv(output, index=False)
    n_sig = (de_results["padj"] < 0.05).sum()
    click.echo(f"DE: {len(de_results)} genes, {n_sig} significant (padj < 0.05) → {output}")
