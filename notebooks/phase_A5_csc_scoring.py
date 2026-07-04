# Phase A5 — CSC Identification & Stemness Scoring
# Input:  data/processed/brca_A4_cancer_epi.h5ad  (24,489 Cancer Epithelial cells)
# Output: data/processed/brca_A5_csc_scored.h5ad
#         results/tables/A5_csc_markers_DE.csv      ← Stage 1 marker list
#
# Steps:
#   1. Score stemness using multiple signatures
#   2. Identify CSC clusters by stemness enrichment
#   3. Run differential expression: CSC-high vs CSC-low
#   4. Produce Stage 1 marker list for G5 comparison

import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os

sc.settings.verbosity = 1
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# ─────────────────────────────────────────────
# STEP 1 — Load Cancer Epithelial subset
# ─────────────────────────────────────────────
print("Loading Cancer Epithelial cells...")
adata = sc.read_h5ad("data/processed/brca_A4_cancer_epi.h5ad")
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Subtypes: {adata.obs['subtype'].value_counts().to_dict()}")
print(f"Minor types: {adata.obs['celltype_minor'].value_counts().to_dict()}")


# ─────────────────────────────────────────────
# STEP 2 — Stemness scoring
# ─────────────────────────────────────────────
# We score each cell using three complementary stemness signatures:
#
# (A) Core pluripotency TF targets — Sox2/Oct4/Nanog/Klf4 regulon
# (B) EMT-stemness signature — epithelial-to-mesenchymal transition genes
#     that co-occur with CSC phenotype in breast cancer
# (C) BRCA-CSC surface markers — experimentally validated in breast CSCs

SIG_PLURIPOTENCY = [
    "SOX2", "POU5F1", "NANOG", "KLF4", "MYC", "LIN28A", "LIN28B",
    "SALL4", "PRDM14", "DNMT3B", "UTF1", "DPPA4", "DPPA2",
    "TDGF1", "ZFP42", "LEFTY1", "LEFTY2",
]

SIG_EMT_STEMNESS = [
    "VIM", "FN1", "CDH2", "TWIST1", "TWIST2", "SNAI1", "SNAI2",
    "ZEB1", "ZEB2", "FOXC2", "GSC", "PDGFRB", "ITGB3",
    "HMGA2", "L1CAM", "CD44", "ALDH1A1", "ALDH1A3",
]

SIG_BRCA_CSC = [
    "CD44", "PROM1", "ITGA6", "EPCAM", "ALDH1A1", "ALDH1A3",
    "BMI1", "EZH2", "NOTCH1", "NOTCH2", "WNT5A", "LGR5",
    "AXIN2", "FZD7", "YAP1", "TEAD1", "FOXM1",
    "SOX9", "KLF5", "ID1", "ID3",
]

def score_signature(adata, gene_list, score_name):
    available = [g for g in gene_list if g in adata.var_names]
    if len(available) < 5:
        print(f"  Warning: only {len(available)} genes available for {score_name}")
    sc.tl.score_genes(adata, available, score_name=score_name, random_state=42)
    print(f"  {score_name}: {len(available)}/{len(gene_list)} genes used")

print("\nScoring stemness signatures...")
score_signature(adata, SIG_PLURIPOTENCY,  "score_pluripotency")
score_signature(adata, SIG_EMT_STEMNESS,  "score_emt_stemness")
score_signature(adata, SIG_BRCA_CSC,      "score_brca_csc")

# Composite stemness score: mean of the three z-scored signatures
for col in ["score_pluripotency", "score_emt_stemness", "score_brca_csc"]:
    mu  = adata.obs[col].mean()
    std = adata.obs[col].std()
    adata.obs[f"{col}_z"] = (adata.obs[col] - mu) / (std + 1e-9)

adata.obs["stemness_composite"] = (
    adata.obs["score_pluripotency_z"] +
    adata.obs["score_emt_stemness_z"] +
    adata.obs["score_brca_csc_z"]
) / 3.0

print("\nStemness composite score summary:")
print(adata.obs["stemness_composite"].describe().to_string())


# ─────────────────────────────────────────────
# STEP 3 — Assign CSC labels
# ─────────────────────────────────────────────
# Top 25% stemness = CSC-high, bottom 25% = CSC-low
# This labelling feeds directly into Phase G1 as well

q75 = adata.obs["stemness_composite"].quantile(0.75)
q25 = adata.obs["stemness_composite"].quantile(0.25)

adata.obs["csc_label"] = "middle"
adata.obs.loc[adata.obs["stemness_composite"] >= q75, "csc_label"] = "csc_high"
adata.obs.loc[adata.obs["stemness_composite"] <= q25, "csc_label"] = "csc_low"

label_counts = adata.obs["csc_label"].value_counts()
print(f"\nCSC label counts: {label_counts.to_dict()}")

# Stemness distribution plot
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for ax, score, label in zip(axes,
    ["score_pluripotency", "score_emt_stemness", "score_brca_csc"],
    ["Pluripotency", "EMT-Stemness", "BRCA-CSC surface"]):
    ax.hist(adata.obs[score], bins=80, color="steelblue", alpha=0.7)
    ax.set_title(f"{label} signature score")
    ax.set_xlabel("Score")
    ax.set_ylabel("Cells")

plt.tight_layout()
plt.savefig("results/figures/A5_stemness_distributions.png", dpi=150)
plt.close()
print("Saved: results/figures/A5_stemness_distributions.png")

# Composite score on UMAP
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.umap(adata, color="stemness_composite",
           color_map="RdYlBu_r", title="Composite Stemness Score",
           ax=axes[0], show=False)

sc.pl.umap(adata, color="csc_label",
           title="CSC Label (top/bottom 25%)",
           palette={"csc_high": "#d62728", "middle": "#aec7e8", "csc_low": "#1f77b4"},
           ax=axes[1], show=False)

sc.pl.umap(adata, color="celltype_minor",
           title="Minor Cell Type",
           legend_loc="right margin", legend_fontsize=8,
           ax=axes[2], show=False)

plt.tight_layout()
plt.savefig("results/figures/A5_stemness_UMAP.png", dpi=150)
plt.close()
print("Saved: results/figures/A5_stemness_UMAP.png")

# Stemness score per minor cell type (violin)
fig, ax = plt.subplots(figsize=(10, 5))
order = (adata.obs.groupby("celltype_minor")["stemness_composite"]
         .median().sort_values(ascending=False).index.tolist())
adata.obs.boxplot(column="stemness_composite", by="celltype_minor",
                  ax=ax, rot=30, figsize=(10, 5))
ax.set_title("Stemness by Cancer Epithelial minor type")
ax.set_xlabel("")
plt.suptitle("")
plt.tight_layout()
plt.savefig("results/figures/A5_stemness_by_minor_type.png", dpi=150)
plt.close()
print("Saved: results/figures/A5_stemness_by_minor_type.png")


# ─────────────────────────────────────────────
# STEP 4 — Differential Expression: CSC-high vs CSC-low
# ─────────────────────────────────────────────
# This is the core of Stage 1. We find genes that are significantly
# upregulated in CSC-high cells compared to CSC-low cells.
# These become the "Stage 1 marker list" for comparison with Geneformer.

print("\nRunning differential expression: CSC-high vs CSC-low...")
adata_labelled = adata[adata.obs["csc_label"] != "middle"].copy()
print(f"DE input: {adata_labelled.n_obs} cells (CSC-high + CSC-low)")

sc.tl.rank_genes_groups(
    adata_labelled,
    groupby="csc_label",
    groups=["csc_high"],
    reference="csc_low",
    method="wilcoxon",
    n_genes=adata_labelled.n_vars,  # rank all genes
    key_added="de_csc",
)

# Extract results into a DataFrame
de_result = sc.get.rank_genes_groups_df(
    adata_labelled,
    group="csc_high",
    key="de_csc",
)
de_result = de_result.rename(columns={
    "names":    "gene_symbol",
    "scores":   "wilcoxon_score",
    "pvals":    "pval",
    "pvals_adj":"pval_adj",
    "logfoldchanges": "log2fc",
})

# Filter to significant upregulated CSC markers
csc_markers = de_result[
    (de_result["pval_adj"] < 0.05) &
    (de_result["log2fc"] > 0.5)
].copy()

print(f"\nSignificant CSC markers (Stage 1): {len(csc_markers)}")
print("\nTop 30 by Wilcoxon score:")
print(csc_markers.head(30)[["gene_symbol", "log2fc", "pval_adj", "wilcoxon_score"]].to_string())


# ─────────────────────────────────────────────
# STEP 5 — Visualize top markers
# ─────────────────────────────────────────────
top_markers = csc_markers["gene_symbol"].head(20).tolist()

# Volcano plot
fig, ax = plt.subplots(figsize=(9, 6))
sig = de_result["pval_adj"] < 0.05
up  = de_result["log2fc"] > 0.5
down= de_result["log2fc"] < -0.5

ax.scatter(de_result.loc[~sig, "log2fc"],
           -np.log10(de_result.loc[~sig, "pval_adj"] + 1e-300),
           s=4, alpha=0.4, color="grey", label="Not significant")
ax.scatter(de_result.loc[sig & up, "log2fc"],
           -np.log10(de_result.loc[sig & up, "pval_adj"] + 1e-300),
           s=6, alpha=0.7, color="#d62728", label=f"CSC-high ({(sig&up).sum()})")
ax.scatter(de_result.loc[sig & down, "log2fc"],
           -np.log10(de_result.loc[sig & down, "pval_adj"] + 1e-300),
           s=6, alpha=0.7, color="#1f77b4", label=f"CSC-low ({(sig&down).sum()})")

# Label top 15 genes
for _, row in csc_markers.head(15).iterrows():
    ax.annotate(row["gene_symbol"], (row["log2fc"],
                -np.log10(row["pval_adj"] + 1e-300)),
                fontsize=7, ha="left")

ax.axvline(0.5,  color="grey", linestyle="--", linewidth=0.8)
ax.axvline(-0.5, color="grey", linestyle="--", linewidth=0.8)
ax.axhline(-np.log10(0.05), color="grey", linestyle=":", linewidth=0.8)
ax.set_xlabel("log2 Fold Change (CSC-high / CSC-low)")
ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("Stage 1 CSC markers — Differential Expression\n(Wilcoxon rank-sum test)")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig("results/figures/A5_volcano_csc_markers.png", dpi=150)
plt.close()
print("Saved: results/figures/A5_volcano_csc_markers.png")

# Dot plot of top markers
sc.pl.dotplot(
    adata,
    var_names=top_markers,
    groupby="celltype_minor",
    standard_scale="var",
    return_fig=True,
).savefig("results/figures/A5_dotplot_top_csc_markers.png",
          dpi=150, bbox_inches="tight")
plt.close()
print("Saved: results/figures/A5_dotplot_top_csc_markers.png")

# UMAP: top 4 CSC markers
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, gene in zip(axes.flat, top_markers[:4]):
    if gene in adata.var_names:
        sc.pl.umap(adata, color=gene, color_map="Reds",
                   title=gene, ax=ax, show=False)
plt.tight_layout()
plt.savefig("results/figures/A5_top4_markers_UMAP.png", dpi=150)
plt.close()
print("Saved: results/figures/A5_top4_markers_UMAP.png")


# ─────────────────────────────────────────────
# STEP 6 — Save Stage 1 marker list + scored object
# ─────────────────────────────────────────────
de_result.to_csv("results/tables/A5_all_DE_results.csv", index=False)
csc_markers.to_csv("results/tables/A5_csc_markers_DE.csv", index=False)

adata.write_h5ad("data/processed/brca_A5_csc_scored.h5ad")

print(f"\nSaved: results/tables/A5_csc_markers_DE.csv  ({len(csc_markers)} Stage 1 markers)")
print(f"Saved: results/tables/A5_all_DE_results.csv")
print(f"Saved: data/processed/brca_A5_csc_scored.h5ad")
print("\n✓ Phase A5 complete.")
print("  Stage 1 pipeline done — CSC marker list produced.")
print("  Next: Phase A6 (trajectory) OR Phase G1 (Geneformer pseudo-labels)")
