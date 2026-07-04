"""cscbench — functional benchmarking of cancer stem cell (CSC) marker-ranking methods.

Score any gene-ranking method (differential expression, transformer attention,
supervised models, network centrality, ...) against independent *functional*
CSC assays (sorted ALDH+/CD44+CD24-/CD133+ populations, sphere self-renewal),
breaking the circularity of validating signature-derived markers against the
same literature that defined them.

Quick start
-----------
>>> from cscbench import FunctionalGate, benchmark_ranking
>>> gate = FunctionalGate.from_matrix(expr, csc_cols, noncsc_cols, name="ALDH prostate")
>>> res  = benchmark_ranking(de_scores, gate)           # de_scores: {gene: score}
>>> res.auroc, res.precision_at_k, res.p_value

Multiple methods x gates
------------------------
>>> from cscbench import run_benchmark
>>> table = run_benchmark({"DE": de_scores, "attention": gf_scores}, [gate1, gate2])
"""

from .functional import FunctionalGate
from .benchmark import BenchmarkResult, benchmark_ranking, run_benchmark

__version__ = "0.1.0"
__all__ = [
    "FunctionalGate",
    "BenchmarkResult",
    "benchmark_ranking",
    "run_benchmark",
]
