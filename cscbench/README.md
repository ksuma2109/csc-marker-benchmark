# cscbench

**Functional benchmarking of cancer stem cell (CSC) marker-ranking methods.**

Computational CSC marker discovery is evaluated circularly: stemness labels are
defined from published markers, and candidate genes are then judged against the
same literature. `cscbench` breaks the circularity by scoring *any* gene-ranking
method against **independent functional CSC assays** — sorted ALDH⁺, CD44⁺CD24⁻,
or CD133⁺ populations and sphere-versus-adherent self-renewal cultures — using
metrics that need no marker list.

It is method-agnostic (differential expression, transformer attention, supervised
models, network centrality — anything that produces a gene ranking) and
cancer-agnostic (any tumor type with a functional dataset).

## Install

```bash
pip install -e .            # from this directory
# optional extras:
pip install -e ".[plot,test]"
```

## Concepts

- **`FunctionalGate`** — a functional assay as a per-gene log2 fold change
  (CSC vs. non-CSC). Build it from an expression matrix + group columns, or from
  a precomputed fold change (e.g. a DESeq2 table).
- **`benchmark_ranking`** — score one ranking against one gate.
- **`run_benchmark`** — score many methods against many gates → tidy DataFrame.

## Metrics

| Metric | Meaning |
|---|---|
| `precision_at_k` | fraction of the method's top-k genes up-regulated in the functional CSC population |
| `p_value` | empirical one-sided p vs. size-matched random gene sets |
| `enrichment` | precision@k / mean random precision |
| `auroc` | genome-wide: does the method's continuous score rank top-functional-quartile genes above the rest? |
| `mean_func_log2fc` | mean functional log2FC of the top-k genes (effect size) |

`precision@k` rewards **shortlist quality**; `auroc` rewards **genome-wide
ranking**. In the accompanying study these dissociate: differential expression
wins precision, transformer attention wins AUROC.

## Quick start

```python
import pandas as pd
from cscbench import FunctionalGate, benchmark_ranking, run_benchmark

# 1. a functional gate from an expression matrix (rows=genes, group columns)
gate = FunctionalGate.from_matrix(
    expr, csc_cols=["ALDHpos_1", "ALDHpos_2", "ALDHpos_3"],
    noncsc_cols=["ALDHneg_1", "ALDHneg_2", "ALDHneg_3"],
    name="Prostate ALDH+", cancer="Prostate", criterion="ALDH",
)

# ...or from a precomputed per-gene log2FC (e.g. DESeq2)
# gate = FunctionalGate.from_log2fc(deseq2_series, name="Mammosphere")

# 2. score a ranking (dict gene->score, or an ordered gene list)
de_scores = {"CD44": 4.1, "SOX9": 3.8, "MYC": 3.5, ...}
res = benchmark_ranking(de_scores, gate)
print(res.auroc, res.precision_at_k, res.p_value)

# 3. many methods x many gates
table = run_benchmark(
    {"DE": de_scores, "attention": attention_scores}, [gate, gate2],
)
print(table.groupby("method")["auroc"].mean())
```

## Reproducing the study benchmark

`examples/reproduce_benchmark.py` rebuilds the multi-cancer benchmark from the
paper using the package API (four methods × functional gates across cancers),
given the processed functional matrices.

## Citation

If you use `cscbench`, please cite the accompanying paper (see repository root).

## License

MIT
