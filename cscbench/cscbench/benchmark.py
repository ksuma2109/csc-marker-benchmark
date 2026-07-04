"""Benchmark metrics — score a gene ranking against a functional gate.

Metrics (all label-independent, computed against a ``FunctionalGate``):

- **precision@k**     fraction of a method's top-k genes that are up-regulated
                      in the functional CSC population, with an empirical p-value
                      vs. size-matched random gene sets.
- **AUROC**           ability of the method's *continuous* score to rank
                      top-functional-quartile genes above the rest, genome-wide.
- **mean_func_log2fc**mean functional log2FC of the top-k genes (effect size).
- **enrichment**      precision@k / mean random precision.

These are exactly the metrics used in the accompanying study; see the paper's
Methods and ``examples/reproduce_benchmark.py``.
"""

from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Mapping, Sequence, Union, Optional, List
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from .functional import FunctionalGate

Ranking = Union[Mapping[str, float], Sequence[str]]


def _as_scores(ranking: Ranking) -> "dict[str, float]":
    """Accept either {gene: score} or an ordered list (best first)."""
    if isinstance(ranking, Mapping):
        return dict(ranking)
    # ordered list -> descending integer scores
    n = len(ranking)
    return {g: float(n - i) for i, g in enumerate(ranking)}


@dataclass
class BenchmarkResult:
    method: str
    gate: str
    cancer: Optional[str]
    criterion: Optional[str]
    n_topk_in_gate: int
    precision_at_k: float
    random_precision: float
    enrichment: float
    p_value: float
    mean_func_log2fc: float
    auroc: float

    def as_dict(self) -> dict:
        return asdict(self)


def benchmark_ranking(
    ranking: Ranking,
    gate: FunctionalGate,
    topk: int = 100,
    n_random: int = 1000,
    up_threshold: float = 0.0,
    random_state: int = 0,
) -> Optional[BenchmarkResult]:
    """Score one ranking against one functional gate.

    Parameters
    ----------
    ranking : dict[str, float] | list[str]
        Gene -> score (higher = stronger CSC marker), or an ordered gene list.
    gate : FunctionalGate
    topk : int
        Number of top genes to evaluate for precision / effect size.
    n_random : int
        Random gene sets for the precision null.
    up_threshold : float
        A gene counts as "up in CSC" if its functional log2FC > this.
    random_state : int
        Seed for the random baseline (reproducible).

    Returns
    -------
    BenchmarkResult, or None if fewer than 10 of the method's ranked genes are
    present in the gate (too little overlap to score).
    """
    rng = np.random.default_rng(random_state)
    scores = _as_scores(ranking)
    func = gate.log2fc
    universe = list(func.index)
    uni_set = set(universe)
    vals = func.values
    n = len(vals)

    # method's ranked genes, ordered by score desc, restricted to the gate
    ordered = [g for g, _ in sorted(scores.items(), key=lambda kv: kv[1], reverse=True)]
    topk_genes = [g for g in ordered if g in uni_set][:topk]
    if len(topk_genes) < 10:
        return None

    prec = float((func.loc[topk_genes] > up_threshold).mean())
    mean_lfc = float(func.loc[topk_genes].mean())

    # genome-wide AUROC: method score predicts top-quartile functional genes
    truth = gate.top_quartile_labels()
    mscore = np.array([scores.get(g, 0.0) for g in universe])
    auroc = float(roc_auc_score(truth.values, mscore)) if truth.nunique() == 2 else float("nan")

    # random-set null for precision@k (sample by position; label-safe)
    k = min(topk, n)
    rand = np.fromiter(
        ((vals[rng.choice(n, k, replace=False)] > up_threshold).mean() for _ in range(n_random)),
        dtype=float, count=n_random,
    )
    rand_mean = float(rand.mean())
    p_value = float((rand >= prec).mean())          # empirical one-sided p
    enrichment = float(prec / rand_mean) if rand_mean > 0 else float("nan")

    return BenchmarkResult(
        method=getattr(ranking, "name", None) or "ranking",
        gate=gate.name, cancer=gate.cancer, criterion=gate.criterion,
        n_topk_in_gate=len(topk_genes),
        precision_at_k=round(prec, 4), random_precision=round(rand_mean, 4),
        enrichment=round(enrichment, 3), p_value=round(p_value, 4),
        mean_func_log2fc=round(mean_lfc, 4), auroc=round(auroc, 4),
    )


def run_benchmark(
    rankings: Mapping[str, Ranking],
    gates: Sequence[FunctionalGate],
    topk: int = 100,
    n_random: int = 1000,
    up_threshold: float = 0.0,
    random_state: int = 0,
) -> pd.DataFrame:
    """Score every ranking against every gate; return a tidy DataFrame.

    Parameters
    ----------
    rankings : dict[str, ranking]
        method name -> ranking (dict or ordered list).
    gates : list[FunctionalGate]

    Returns
    -------
    pandas.DataFrame with one row per (method, gate) and all metric columns.
    """
    rows: List[dict] = []
    for method, ranking in rankings.items():
        for gate in gates:
            res = benchmark_ranking(ranking, gate, topk=topk, n_random=n_random,
                                    up_threshold=up_threshold, random_state=random_state)
            if res is None:
                continue
            d = res.as_dict()
            d["method"] = method
            rows.append(d)
    cols = ["method", "gate", "cancer", "criterion", "precision_at_k", "random_precision",
            "enrichment", "p_value", "mean_func_log2fc", "auroc", "n_topk_in_gate"]
    df = pd.DataFrame(rows)
    return df[[c for c in cols if c in df.columns]] if len(df) else df
