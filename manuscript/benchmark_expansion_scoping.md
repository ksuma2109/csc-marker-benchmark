# Benchmark-Expansion Scoping — strengthening the functional CSC benchmark before submission

**Purpose.** The current functional benchmark uses 3 gates (GSE115302, GSE36643, GSE182532), concentrated in breast cancer, with thin replication (one PDX set at n=1/group) and one negative gate (mammosphere). The two most likely reviewer objections at *Briefings in Bioinformatics* are (1) **benchmark scale/scope** and (2) **generalizability beyond one cancer**. This document scopes three expansions that address them, with effort and impact, drawn from 14 GEO datasets verified in scouting (10 cancer types; ALDH, sphere, CD133 criteria).

---

## Expansion 1 — Broaden the functional benchmark (HIGH impact, LOW–MODERATE effort) ⭐ recommended first

**What:** Add ready-to-use RNA-seq functional datasets to the existing benchmark's contrast set. Takes the benchmark from 3 gates / 1 cancer to **8 gates across ~5–6 cancers with proper n=3 replicates**, spanning all three functional criteria.

**Datasets (already verified, processed matrices downloadable, n=3/group):**

| Accession | Cancer | Criterion | Processed file | Effort |
|---|---|---|---|---|
| GSE270565 | Prostate | ALDH high/low | `RNAseq_processed_data.csv.gz` | download + ID map |
| GSE243840 | Melanoma | ALDH1A3 high/low | `High_vs_Low_tpm.csv.gz` | download + ID map |
| GSE228203 | Prostate | sphere vs adherent | `gene_fpkm.txt.gz` | download + ID map |
| GSE232783 | Ovarian | sphere vs adherent | `processed_data_file.xls.gz` | download + ID map |
| GSE166947 | Bladder | ALDH+/- | `Gene_expression_cpm.txt.gz` | download + ID map |

**Why it works:** the F1/F4 benchmark already computes everything from a per-gene functional log2FC. Each new dataset is just another entry in the `contrasts` dict: load matrix → compute log2(CSC+1)/(nonCSC+1) → done. No model retraining. The DE/Geneformer/LR/RF rankings are fixed (breast-derived); we simply test them against more, and more diverse, functional ground truths.

**Effort:** ~half a day. Main work is per-dataset I/O (each has a different column layout) and gene-ID mapping. Picks 1–5 need no normalization.

**Impact on the paper:**
- Directly rebuts "benchmark too small / breast-only."
- Tests whether the shortlist-vs-ranker dissociation **holds across cancers** — if Geneformer's AUROC lead persists on prostate/ovarian/bladder/melanoma gates, that is a much stronger, more general claim (and a headline BiB result).
- Lets the mammosphere-negative result be contextualized against 2 more sphere gates (prostate, ovarian), showing whether sphere assays are systematically weaker anchors.
- Risk: the breast-derived rankings may score *lower* against non-breast gates (cancer-specific CSC biology). That is itself a publishable, honest finding ("marker rankings are partly cancer-specific; the method-level dissociation generalizes even where individual genes do not").

**Optional add-ons (more breadth, more effort):** GSE207333 (osteosarcoma CD133) and GSE283717 (GBM CD133) are RNA-seq but ship per-sample files in `RAW.tar` (add a concatenation step). GSE52262 is a breast dual-gate (CD44+CD24- and ALDH) but CEL-only (needs array normalization).

---

## Expansion 2 — Second discovery cancer (MODERATE–HIGH impact, HIGH effort)

**What:** Re-run the full four-method pipeline (DE + Geneformer + LR + RF) on a *second* cancer's scRNA-seq to derive an independent set of rankings, then benchmark those against that cancer's functional gates. We already have processed GBM (GSE84465) and melanoma (GSE72056) scRNA-seq from the cross-cancer analysis, and matching functional gates exist (GSE283717 GBM CD133; GSE243840 melanoma ALDH).

**Why it's stronger but heavier:** this shows the *method-selection conclusions themselves* replicate in a second tumor, not just that breast markers appear elsewhere. That is the most convincing generalizability evidence.

**Effort:** HIGH — the bottleneck is **Geneformer fine-tuning on the new cancer**, which is compute-heavy (the breast run took hours; realistically needs GPU/Colab). DE, LR, and RF on a new cancer are cheap (minutes). So the cost is almost entirely one more transformer fine-tune per cancer.

**Impact:** moderate–high. Valuable, but Expansion 1 delivers most of the generalizability benefit at a fraction of the cost. Do this only if Expansion 1's cross-cancer scoring leaves reviewers wanting a full second-cancer replication.

---

## Expansion 3 — Package the framework as a released tool (HIGH impact for framing, MODERATE effort)

**What:** Refactor the F1/F3/F4 scripts into a documented, installable package with a clean API — input: a gene ranking + a functional dataset (matrix + group labels); output: precision@k, AUROC, effect-size concordance, random-baseline p-values, and the summary figure. Add a repo, README, example notebook, and minimal tests.

**Why it matters for BiB (and opens Bioinformatics):** the manuscript now *claims* a "functional-benchmarking framework." A released tool makes that claim concrete rather than rhetorical — the single biggest lever for a benchmark paper's credibility. It would also open the **Bioinformatics Applications Note** route (tool-focused), which was previously closed because nothing was packaged.

**Effort:** MODERATE. The logic exists and is tested-by-use; the work is API design, packaging (pip/conda), docs, and a small test suite. No new science.

**Impact:** HIGH for the "framework" positioning; converts a scripts-based analysis into a reusable community resource, which is exactly what BiB benchmark papers are expected to provide.

---

## Recommended sequence

1. **Expansion 1** (broaden the benchmark) — highest impact-to-effort; do first. Rewrites §3.3/§2.6 with an 8-gate, multi-cancer result and turns the biggest weakness into a strength.
2. **Expansion 3** (package the framework) — makes the paper's central claim real and is BiB's expectation for a benchmark; moderate effort, no new science.
3. **Expansion 2** (second discovery cancer) — only if reviewers demand full replication; high compute cost.

Expansions 1 and 3 together would move the paper from "modest single-cancer benchmark" to "multi-cancer benchmark with a released framework" — squarely in BiB's wheelhouse and, with the packaged tool, potentially viable as a Bioinformatics Applications Note.

**Immediate next action if approved:** implement Expansion 1 — the 5 ready datasets plug directly into the existing benchmark code (`phase_F1`/`phase_F4` `contrasts` dict) with no model retraining.
