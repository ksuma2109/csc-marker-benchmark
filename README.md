# Functional benchmarking of cancer stem cell marker-ranking methods

A computational study and reusable framework for identifying cancer stem cell
(CSC) marker genes from single-cell RNA-seq, benchmarking differential
expression against transformer attention (and supervised baselines) using
**independent functional CSC assays** as a label-independent standard.

- **Manuscript:** [`manuscript/manuscript_draft.md`](manuscript/manuscript_draft.md) · rendered PDF [`manuscript/CSC_manuscript.pdf`](manuscript/CSC_manuscript.pdf)
- **Framework (installable package):** [`cscbench/`](cscbench/)
- **Analysis scripts:** [`notebooks/`](notebooks/)
- **Derived results:** [`results/tables/`](results/tables/), figures in [`manuscript/figures/`](manuscript/figures/)

## What this is

Computational CSC marker discovery is evaluated circularly: stemness labels are
defined from published markers, and candidate genes are then judged against the
same literature. This project breaks the circularity by scoring gene rankings
against functional assays that separate cells by phenotype (sorted ALDH⁺,
CD44⁺CD24⁻, CD133⁺; sphere self-renewal). Applied to a breast cancer atlas and
extended across six cancers, the benchmark shows a systematic dissociation —
differential expression produces higher-precision shortlists, transformer
attention is a stronger genome-wide ranker — while revealing that marker
rankings are substantially tissue-specific.

## The `cscbench` package

```bash
pip install -e cscbench
```

```python
from cscbench import FunctionalGate, run_benchmark
gate = FunctionalGate.from_matrix(expr, csc_cols, noncsc_cols, name="ALDH prostate")
table = run_benchmark({"DE": de_scores, "attention": attn_scores}, [gate])
```

Method- and cancer-agnostic; scores any gene ranking against any functional
gate (precision@k with a random-set null, genome-wide AUROC, effect-size
concordance). See [`cscbench/README.md`](cscbench/README.md).

## Reproducing

Analysis scripts run in phase order (`notebooks/phase_*.py`). Raw data are
public (GEO / cBioPortal) and excluded from the repository; accessions are in
the manuscript's Data Availability statement. `cscbench/examples/reproduce_benchmark.py`
rebuilds the multi-cancer benchmark from the processed functional matrices.

## Repository layout

| Path | Contents |
|---|---|
| `cscbench/` | Installable benchmarking framework + tests |
| `notebooks/` | Analysis scripts (single-cell pipeline, Geneformer, benchmarks, survival, figures) |
| `manuscript/` | Manuscript, figures, rendered PDF |
| `results/tables/` | Derived result tables (CSV) |

Large data, model weights, and the Python environment are git-ignored (see
`.gitignore`).

## License

Code and the `cscbench` package are released under the MIT License.
