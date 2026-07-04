# Phase A6 — Trajectory & Pseudotime Analysis
# Input:  data/processed/brca_A5_csc_scored.h5ad  (Cancer Epithelial, 24,489 cells)
# Output: data/processed/brca_A6_trajectory.h5ad
#
# Approach:
#   Diffusion pseudotime (DPT) — places each cell on a continuous
#   differentiation axis from stem-like to differentiated.
#   DPT requires only spliced counts so works without velocyto preprocessing.
#
#   RNA velocity (scVelo) requires unspliced/spliced read counts from the
#   BAM files (STARsolo or velocyto pipeline). Since we only have the count
#   matrix (no BAM), we run DPT here and note where to add scVelo when
#   BAM files are available.

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

sc.settings.verbosity = 1
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

# ─────────────────────────────────────────────
# STEP 1 — Load scored data
# ─────────────────────────────────────────────
print("Loading CSC-scored Cancer Epithelial cells...")
adata = sc.read_h5ad("data/processed/brca_A5_csc_scored.h5ad")
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")


# ─────────────────────────────────────────────
# STEP 2 — Recompute PCA + neighbors on subset
# ─────────────────────────────────────────────
# The neighbor graph was built on all 100k cells. Subsetting to Cancer
# Epithelial cells breaks connectivity (neighbors from other cell types
# are removed), causing DPT to return inf. We must rebuild the graph
# within this subset.

print("Recomputing PCA on Cancer Epithelial subset...")
# Scale HVGs and rerun PCA
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=30, svd_solver="arpack")
adata.obsm["X_pca_subset"] = adata_hvg.obsm["X_pca"]

print("Recomputing neighbors on Cancer Epithelial subset...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30,
                use_rep="X_pca_subset", random_state=42)

print("Recomputing UMAP on subset...")
sc.tl.umap(adata, random_state=42)

# ─────────────────────────────────────────────
# STEP 3a — Diffusion Map
# ─────────────────────────────────────────────
print("Computing diffusion map...")
sc.tl.diffmap(adata, n_comps=15)
print("Diffusion map complete.")


# ─────────────────────────────────────────────
# STEP 3b — Diffusion Pseudotime (DPT)
# ─────────────────────────────────────────────
# DPT assigns each cell a pseudotime value representing how far it is
# from a chosen root cell along the differentiation trajectory.
#
# Root cell selection: the cell with the highest composite stemness score
# is the most stem-like — it is the natural root of the trajectory.

root_idx = adata.obs["stemness_composite"].values.argmax()
print(f"Root cell: index {root_idx}  "
      f"(stemness = {adata.obs['stemness_composite'].iloc[root_idx]:.3f}, "
      f"minor type = {adata.obs['celltype_minor'].iloc[root_idx]})")

adata.uns["iroot"] = root_idx
sc.tl.dpt(adata, n_dcs=10)
print("DPT complete.")

print("\nPseudotime summary:")
print(adata.obs["dpt_pseudotime"].describe().to_string())


# ─────────────────────────────────────────────
# STEP 4 — Visualise trajectory
# ─────────────────────────────────────────────
print("Generating trajectory plots...")

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Row 1
sc.pl.umap(adata, color="dpt_pseudotime",
           color_map="viridis",
           title="Diffusion Pseudotime",
           ax=axes[0, 0], show=False)

sc.pl.umap(adata, color="stemness_composite",
           color_map="RdYlBu_r",
           title="Composite Stemness Score",
           ax=axes[0, 1], show=False)

sc.pl.umap(adata, color="csc_label",
           palette={"csc_high": "#d62728", "middle": "#aec7e8", "csc_low": "#1f77b4"},
           title="CSC Label",
           ax=axes[0, 2], show=False)

# Row 2
sc.pl.umap(adata, color="celltype_minor",
           title="Minor Cell Type",
           legend_loc="right margin", legend_fontsize=8,
           ax=axes[1, 0], show=False)

sc.pl.umap(adata, color="subtype",
           title="Breast Cancer Subtype",
           ax=axes[1, 1], show=False)

# Pseudotime vs stemness scatter
ax = axes[1, 2]
sample = adata.obs.sample(min(5000, adata.n_obs), random_state=42)
sc.pl.scatter(adata, x="dpt_pseudotime", y="stemness_composite",
              color="celltype_minor",
              title="Pseudotime vs Stemness",
              ax=ax, show=False)

plt.tight_layout()
plt.savefig("results/figures/A6_trajectory_overview.png", dpi=150)
plt.close()
print("Saved: results/figures/A6_trajectory_overview.png")


# ─────────────────────────────────────────────
# STEP 5 — Pseudotime distribution per cell type
# ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))

minor_types = adata.obs["celltype_minor"].unique()
colors = plt.cm.tab10(np.linspace(0, 1, len(minor_types)))

for ct, color in zip(minor_types, colors):
    subset = adata.obs.loc[adata.obs["celltype_minor"] == ct, "dpt_pseudotime"]
    ax.hist(subset, bins=50, alpha=0.6, label=ct, color=color, density=True)

ax.set_xlabel("Diffusion Pseudotime")
ax.set_ylabel("Density")
ax.set_title("Pseudotime distribution by Cancer Epithelial minor type\n"
             "(left = stem-like, right = differentiated)")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig("results/figures/A6_pseudotime_distribution.png", dpi=150)
plt.close()
print("Saved: results/figures/A6_pseudotime_distribution.png")


# ─────────────────────────────────────────────
# STEP 6 — Pseudotime vs stemness correlation
# ─────────────────────────────────────────────
from scipy.stats import spearmanr

corr, pval = spearmanr(
    adata.obs["dpt_pseudotime"],
    adata.obs["stemness_composite"]
)
print(f"\nPseudotime vs stemness correlation: Spearman r = {corr:.3f}, p = {pval:.2e}")

# Median pseudotime per minor cell type
pt_by_type = (adata.obs.groupby("celltype_minor")["dpt_pseudotime"]
              .median().sort_values())
print("\nMedian pseudotime per minor cell type (lower = more stem-like):")
print(pt_by_type.to_string())
pt_by_type.to_csv("results/tables/A6_pseudotime_by_celltype.csv")
print("Saved: results/tables/A6_pseudotime_by_celltype.csv")


# ─────────────────────────────────────────────
# STEP 7 — Genes correlated with pseudotime
# ─────────────────────────────────────────────
# Find genes whose expression changes monotonically along pseudotime.
# These are the genes driving the CSC → differentiated transition.

print("\nComputing gene-pseudotime correlations...")
from scipy.stats import spearmanr as sr

# Use HVGs only (2000 genes) for speed
hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()
pt_vals   = adata.obs["dpt_pseudotime"].values

X = adata[:, hvg_genes].X
if hasattr(X, "toarray"):
    X = X.toarray()

correlations = []
for i, gene in enumerate(hvg_genes):
    r, p = sr(X[:, i], pt_vals)
    correlations.append({"gene": gene, "spearman_r": r, "pval": p})
    if i % 200 == 0:
        print(f"  {i}/{len(hvg_genes)} genes processed...")

corr_df = pd.DataFrame(correlations).sort_values("spearman_r")

# Genes decreasing along pseudotime = expressed in stem-like cells
stem_genes_pt = corr_df[corr_df["pval"] < 0.001].head(30)
# Genes increasing along pseudotime = expressed in differentiated cells
diff_genes_pt = corr_df[corr_df["pval"] < 0.001].tail(30)

print(f"\nTop 10 genes DECREASING along pseudotime (CSC-associated):")
print(stem_genes_pt.head(10)[["gene", "spearman_r", "pval"]].to_string())

print(f"\nTop 10 genes INCREASING along pseudotime (differentiation-associated):")
print(diff_genes_pt.tail(10).sort_values("spearman_r", ascending=False)
      [["gene", "spearman_r", "pval"]].to_string())

corr_df.to_csv("results/tables/A6_gene_pseudotime_correlation.csv", index=False)
print("\nSaved: results/tables/A6_gene_pseudotime_correlation.csv")

# Plot expression of top 4 pseudotime-correlated CSC genes along trajectory
top_stem_pt = stem_genes_pt["gene"].head(4).tolist()

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
for ax, gene in zip(axes.flat, top_stem_pt):
    if gene in adata.var_names:
        sc.pl.scatter(adata, x="dpt_pseudotime", y=gene,
                      color="celltype_minor",
                      title=f"{gene} vs Pseudotime",
                      ax=ax, show=False)
plt.tight_layout()
plt.savefig("results/figures/A6_csc_genes_vs_pseudotime.png", dpi=150)
plt.close()
print("Saved: results/figures/A6_csc_genes_vs_pseudotime.png")


# ─────────────────────────────────────────────
# STEP 8 — Save
# ─────────────────────────────────────────────
adata.write_h5ad("data/processed/brca_A6_trajectory.h5ad")
print("\nSaved: data/processed/brca_A6_trajectory.h5ad")
print("\n✓ Phase A6 complete — Stage 1 pipeline finished.")
print("\nStage 1 outputs:")
print("  results/tables/A5_csc_markers_DE.csv     — 3028 CSC markers (DE method)")
print("  results/tables/A6_gene_pseudotime_correlation.csv — trajectory gene dynamics")
print("\nReady for Stage 2: Phase G1 (Geneformer pseudo-labels)")
