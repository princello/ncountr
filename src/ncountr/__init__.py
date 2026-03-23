"""ncountr — Python pipeline for Nanostring nCounter data analysis."""

__version__ = "0.1.0"

from ncountr.experiment import NanostringExperiment
from ncountr.io.rcc import read_rcc, parse_rcc
from ncountr.io.export import to_anndata
from ncountr.core.qc import qc
from ncountr.core.normalize import normalize
from ncountr.core.de import de
from ncountr.core.pathway import score_gene_set
from ncountr.datasets import get_gene_set, list_gene_sets

__all__ = [
    "NanostringExperiment",
    "read_rcc",
    "parse_rcc",
    "qc",
    "normalize",
    "de",
    "score_gene_set",
    "get_gene_set",
    "list_gene_sets",
    "to_anndata",
]
