# Phase F6 — CRISPR-dependency validation (DepMap)
#
# Independent perturbation-based validation of the CSC marker rankings, using
# public genome-wide CRISPR-knockout fitness screens (DepMap 24Q2). For each
# gene, DepMap reports a Chronos "gene effect" per cell line (more negative =
# knockout reduces fitness = the line depends on the gene). We ask, in BREAST
# cancer cell lines: are the predicted CSC markers enriched for CRISPR
# dependencies, and does the Geneformer (attention) ranking capture more
# dependency signal than differential expression?
#
# We frame CRISPR dependency as another functional gate for cscbench:
#   gate value = -(mean gene effect across breast cancer lines)
#   => higher value = more essential; "up" (> 0.5) = a DepMap dependency.
#
# HONEST SCOPE: DepMap measures proliferation/fitness essentiality in 2D
# culture, NOT stemness/self-renewal. Many bona fide CSC markers are
# context-specific and are NOT common-essential, so LOW absolute precision is
# expected; the informative comparison is RELATIVE (methods vs. random, and
# vs. each other) plus per-candidate dependency.
#
# Output: results/tables/F6_crispr_validation.csv
#         results/tables/F6_candidate_dependencies.csv
#         results/figures/F6_crispr_validation.png

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cscbench"))
from cscbench import FunctionalGate, benchmark_ranking, run_benchmark

os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

print("=" * 64)
print("F6 — CRISPR-DEPENDENCY VALIDATION (DepMap 24Q2, breast lines)")
print("=" * 64)

# ── breast cancer cell lines ────────────────────────────────────────────────
model = pd.read_csv("data/raw/depmap/Model.csv")
breast = model[(model["OncotreeLineage"] == "Breast") &
               (model["OncotreePrimaryDisease"] != "Non-Cancerous")]
breast_ids = set(breast["ModelID"])
print(f"Breast cancer cell lines: {len(breast_ids)}")

# ── CRISPR gene-effect matrix (rows = ModelID, cols = 'SYMBOL (entrez)') ─────
print("Loading CRISPRGeneEffect.csv (large)...")
ge = pd.read_csv("data/raw/depmap/CRISPRGeneEffect.csv", index_col=0)
ge = ge.loc[ge.index.isin(breast_ids)]
print(f"  Breast lines with CRISPR data: {ge.shape[0]}  |  genes: {ge.shape[1]}")

# parse gene symbols from column headers, mean effect across breast lines
ge.columns = [c.split(" (")[0] for c in ge.columns]
ge = ge.loc[:, ~ge.columns.duplicated()]
mean_effect = ge.mean(axis=0)                       # per-gene mean Chronos effect
dependency = -mean_effect                            # higher = more essential
dependency = dependency[dependency.index.astype(str).str.len() > 0].dropna()
n_dep = int((mean_effect < -0.5).sum())
print(f"  Genes that are dependencies (mean effect < -0.5): {n_dep}")

# ── build the CRISPR gate and load method rankings ──────────────────────────
gate = FunctionalGate.from_log2fc(
    dependency, name="Breast CRISPR dependency (DepMap)",
    cancer="Breast", criterion="CRISPR")

s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
lr = pd.read_csv("results/tables/F4_logreg_ranking.csv"); lr = lr[lr["lr_coef"] > 0]
rf = pd.read_csv("results/tables/F4_rf_ranking.csv");    rf = rf[rf["rf_signed"] > 0]
RANKINGS = {
    "Stage1_DE":         dict(zip(s1["gene_symbol"], s1["wilcoxon_score"])),
    "Stage2_Geneformer": dict(zip(s2["gene_symbol"], s2["attention_score"])),
    "LogReg":            dict(zip(lr["gene_symbol"], lr["lr_coef"])),
    "RandomForest":      dict(zip(rf["gene_symbol"], rf["rf_importance"])),
}

# ── benchmark: dependency as a gate (up_threshold 0.5 = DepMap dependency) ───
print("\nScoring methods against the CRISPR-dependency gate...")
rows = []
for m, r in RANKINGS.items():
    res = benchmark_ranking(r, gate, topk=100, n_random=1000, up_threshold=0.5, random_state=42)
    if res:
        d = res.as_dict(); d["method"] = m; rows.append(d)
bench = pd.DataFrame(rows)[
    ["method", "precision_at_k", "random_precision", "enrichment", "p_value",
     "mean_func_log2fc", "auroc", "n_topk_in_gate"]]
bench = bench.rename(columns={"precision_at_k": "frac_dependency@100",
                              "random_precision": "random_frac",
                              "mean_func_log2fc": "mean_neg_effect@100",
                              "auroc": "AUROC_essentiality"})
bench.to_csv("results/tables/F6_crispr_validation.csv", index=False)
print(bench.to_string(index=False))

print("\n  AUROC (method score predicts top-quartile essentiality):")
for _, r in bench.sort_values("AUROC_essentiality", ascending=False).iterrows():
    print(f"    {r['method']:18s} AUROC={r['AUROC_essentiality']:.3f}  "
          f"dependency-fraction@100={r['frac_dependency@100']:.2f} "
          f"(random {r['random_frac']:.2f}, {r['enrichment']}x, p={r['p_value']:.3f})")

# ── per-candidate dependency (is each key gene a breast dependency?) ─────────
print("\n" + "=" * 64)
print("PER-CANDIDATE CRISPR DEPENDENCY (breast lines)")
print("=" * 64)
WATCH = ["SOX9", "KLF4", "FZD7", "BMPR1B", "KLK5", "FOXI1", "CD44", "MYC", "VIM",
         "FN1", "EPCAM", "SERPINE2", "OSMR", "ALDH1A3", "EGFR"]
cand_rows = []
for g in WATCH:
    if g in mean_effect.index:
        eff = float(mean_effect[g])
        cand_rows.append({"gene": g, "mean_gene_effect": round(eff, 3),
                          "is_dependency": eff < -0.5,
                          "de_rank": (s1["gene_symbol"].tolist().index(g)+1) if g in s1["gene_symbol"].values else None,
                          "gf_rank": (s2["gene_symbol"].tolist().index(g)+1) if g in s2["gene_symbol"].values else None})
cand = pd.DataFrame(cand_rows).sort_values("mean_gene_effect")
cand.to_csv("results/tables/F6_candidate_dependencies.csv", index=False)
print(f"  {'gene':10s} {'mean_effect':>11s} {'dependency?':>11s} {'DE#':>5} {'GF#':>5}")
for _, r in cand.iterrows():
    dep = "YES" if r["is_dependency"] else "no"
    de = str(int(r["de_rank"])) if pd.notna(r["de_rank"]) else "-"
    gf = str(int(r["gf_rank"])) if pd.notna(r["gf_rank"]) else "-"
    print(f"  {r['gene']:10s} {r['mean_gene_effect']:>11.3f} {dep:>11s} {de:>5} {gf:>5}")

# ── figure ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
colors = {"Stage1_DE":"#d62728","Stage2_Geneformer":"#1f77b4","LogReg":"#2ca02c","RandomForest":"#ff7f0e"}
b = bench.set_index("method")
methods = list(RANKINGS.keys())
axes[0].bar(range(len(methods)), [b.loc[m,"AUROC_essentiality"] for m in methods],
            color=[colors[m] for m in methods])
axes[0].axhline(0.5, color="k", ls="--", lw=0.8)
axes[0].set_xticks(range(len(methods))); axes[0].set_xticklabels(methods, rotation=20, ha="right", fontsize=8)
axes[0].set_ylabel("AUROC (predicts CRISPR essentiality)"); axes[0].set_ylim(0,0.75)
axes[0].set_title("Do CSC marker rankings predict\nCRISPR dependency? (DepMap breast lines)")

cand_sorted = cand.sort_values("mean_gene_effect")
bar_c = ["#d62728" if v else "#999" for v in cand_sorted["is_dependency"]]
axes[1].barh(cand_sorted["gene"], cand_sorted["mean_gene_effect"], color=bar_c)
axes[1].axvline(-0.5, color="red", ls="--", lw=0.8, label="dependency cutoff (-0.5)")
axes[1].axvline(0, color="k", lw=0.5)
axes[1].set_xlabel("mean CRISPR gene effect (breast lines; <-0.5 = dependency)")
axes[1].set_title("Per-candidate CRISPR dependency")
axes[1].legend(fontsize=8)
plt.suptitle("Figure. CRISPR-dependency validation of CSC marker predictions (DepMap 24Q2)", fontsize=12)
plt.tight_layout()
plt.savefig("results/figures/F6_crispr_validation.png", dpi=130, bbox_inches="tight")
plt.close()

print("\n" + "=" * 64)
print("PHASE F6 COMPLETE")
print("=" * 64)
print("  Saved: results/tables/F6_crispr_validation.csv")
print("  Saved: results/tables/F6_candidate_dependencies.csv")
print("  Saved: results/figures/F6_crispr_validation.png")
