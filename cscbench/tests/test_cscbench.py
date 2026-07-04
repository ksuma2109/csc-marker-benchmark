"""Tests for cscbench — synthetic data, no external files."""
import numpy as np
import pandas as pd
import pytest
from cscbench import FunctionalGate, benchmark_ranking, run_benchmark


def _synthetic_gate(n_genes=200, n_csc_up=50, seed=0):
    """A gate where genes g0..g49 are strongly up in CSC; the rest are noise.

    n_genes=200 with 50 up-genes means the top functional quartile (top 50 by
    log2FC) is exactly the up-set, so a ranking that scores the up-genes first
    is a genuinely 'perfect' ranking with respect to the quartile truth.
    """
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    base = rng.normal(5, 1, (n_genes, 6))          # 3 CSC + 3 non-CSC columns
    base[:n_csc_up, :3] += 4.0                       # up-regulate first 50 in CSC cols
    expr = pd.DataFrame(base, index=genes,
                        columns=["c1", "c2", "c3", "n1", "n2", "n3"])
    return expr, genes, n_csc_up


def test_gate_from_matrix_direction():
    expr, genes, n_up = _synthetic_gate()
    gate = FunctionalGate.from_matrix(expr, ["c1", "c2", "c3"], ["n1", "n2", "n3"],
                                      name="synthetic")
    # the up-regulated genes should have positive log2FC
    assert (gate.log2fc.loc[[f"g{i}" for i in range(n_up)]] > 0).all()
    assert len(gate) == len(genes)


def test_from_log2fc_roundtrip():
    s = pd.Series({"CD44": 2.0, "SOX9": 1.5, "ACTB": -0.1})
    gate = FunctionalGate.from_log2fc(s, name="x")
    assert gate.log2fc["CD44"] == 2.0
    assert "ACTB" in gate.genes


def test_perfect_ranking_beats_random():
    expr, genes, n_up = _synthetic_gate()
    gate = FunctionalGate.from_matrix(expr, ["c1", "c2", "c3"], ["n1", "n2", "n3"], name="s")
    # a ranking that puts the true up-genes first should score high
    good = {f"g{i}": float(n_up - i) for i in range(n_up)}
    res = benchmark_ranking(good, gate, topk=50, n_random=500)
    assert res is not None
    assert res.precision_at_k == 1.0          # all top-50 are up in CSC
    assert res.auroc > 0.9                    # discriminates strongly
    assert res.p_value < 0.05                 # beats random
    assert res.enrichment > 1.0


def test_random_ranking_is_near_chance():
    expr, genes, n_up = _synthetic_gate()
    gate = FunctionalGate.from_matrix(expr, ["c1", "c2", "c3"], ["n1", "n2", "n3"], name="s")
    rng = np.random.default_rng(1)
    rand_scores = {g: float(rng.normal()) for g in genes}
    res = benchmark_ranking(rand_scores, gate, topk=50, n_random=500)
    assert res is not None
    assert 0.3 < res.auroc < 0.7              # near chance
    assert res.p_value > 0.05                 # not enriched


def test_ordered_list_input():
    expr, genes, n_up = _synthetic_gate()
    gate = FunctionalGate.from_matrix(expr, ["c1", "c2", "c3"], ["n1", "n2", "n3"], name="s")
    ordered = [f"g{i}" for i in range(n_up)] + [f"g{i}" for i in range(n_up, 500)]
    res = benchmark_ranking(ordered, gate, topk=50, n_random=200)
    assert res.precision_at_k == 1.0


def test_insufficient_overlap_returns_none():
    gate = FunctionalGate.from_log2fc(pd.Series({"A": 1.0, "B": -1.0}), name="tiny")
    res = benchmark_ranking({"X": 1.0, "Y": 2.0}, gate)   # no overlap
    assert res is None


def test_run_benchmark_shape():
    expr, genes, n_up = _synthetic_gate()
    g1 = FunctionalGate.from_matrix(expr, ["c1", "c2", "c3"], ["n1", "n2", "n3"],
                                    name="g1", cancer="A", criterion="ALDH")
    g2 = FunctionalGate.from_matrix(expr, ["c1", "c2", "c3"], ["n1", "n2", "n3"],
                                    name="g2", cancer="B", criterion="sphere")
    good = {f"g{i}": float(n_up - i) for i in range(n_up)}
    rng = np.random.default_rng(2)
    rand = {g: float(rng.normal()) for g in genes}
    table = run_benchmark({"good": good, "random": rand}, [g1, g2],
                          topk=50, n_random=200)
    assert set(table["method"]) == {"good", "random"}
    assert len(table) == 4                     # 2 methods x 2 gates
    # the good ranking should out-AUROC the random one
    good_auc = table[table.method == "good"]["auroc"].mean()
    rand_auc = table[table.method == "random"]["auroc"].mean()
    assert good_auc > rand_auc


def test_reproducible_with_seed():
    gate = FunctionalGate.from_log2fc(
        pd.Series({f"g{i}": np.sin(i) for i in range(200)}), name="s")
    scores = {f"g{i}": float(i) for i in range(200)}
    a = benchmark_ranking(scores, gate, random_state=42, n_random=300)
    b = benchmark_ranking(scores, gate, random_state=42, n_random=300)
    assert a.p_value == b.p_value and a.precision_at_k == b.precision_at_k
