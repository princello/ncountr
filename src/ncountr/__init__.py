"""ncountr — Python pipeline for Nanostring nCounter data analysis."""

__version__ = "0.2.0"

from ncountr.experiment import NanostringExperiment
from ncountr.io.rcc import read_rcc, parse_rcc
from ncountr.io.export import to_anndata
from ncountr.core.qc import qc
from ncountr.core.normalize import normalize
from ncountr.core.de import de, effect_sizes
from ncountr.core.pathway import score_gene_set
from ncountr.core.gsea import gsea, gsea_preranked, competitive_test, self_contained_test
from ncountr.datasets import get_gene_set, list_gene_sets
from ncountr.genesets import load_gmt, list_hallmark_sets

__all__ = [
    "NanostringExperiment",
    "read_rcc",
    "parse_rcc",
    "qc",
    "normalize",
    "de",
    "effect_sizes",
    "score_gene_set",
    "gsea",
    "gsea_preranked",
    "competitive_test",
    "self_contained_test",
    "get_gene_set",
    "list_gene_sets",
    "list_hallmark_sets",
    "load_gmt",
    "to_anndata",
]
