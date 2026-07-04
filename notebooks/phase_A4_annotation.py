# Phase A4 — Cell Type Annotation
# Input:  data/processed/brca_A3_clustered.h5ad
# Output: data/processed/brca_A4_annotated.h5ad
#         data/processed/brca_A4_cancer_epi.h5ad  (Cancer Epithelial subset)
#
# Wu et al. 2021 already provides expert annotations in metadata.csv
# (celltype_major, celltype_minor, celltype_subset).
# This phase:
#   1. Verifies those annotations align with our Leiden clusters
#   2. Maps each cluster to its dominant cell type
#   3. Checks known marker gene expression per cluster
#   4. Subsets to Cancer Epithelial cells for downstream CSC analysis

import scanpy as sc
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
# STEP 1 — Load clustered data
# ─────────────────────────────────────────────
print("Loading clustered data...")
adata = sc.read_h5ad("data/processed/brca_A3_clustered.h5ad")
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"\nCelltype major counts:")
print(adata.obs["celltype_major"].value_counts().to_string())


# ─────────────────────────────────────────────
# STEP 2 — Cross-tabulate clusters vs annotations
# ─────────────────────────────────────────────
# For each Leiden cluster, calculate the fraction of cells from each
# annotated cell type. A well-formed cluster should be dominated by
# one cell type (>60%), confirming that our clustering aligns with biology.

print("\n--- Cluster purity check ---")
crosstab = pd.crosstab(
    adata.obs["leiden_0.3"],
    adata.obs["celltype_major"],
    normalize="index",      # row-normalize: fraction per cluster
)

# Dominant cell type per cluster and its purity
crosstab["dominant_celltype"] = crosstab.idxmax(axis=1)
crosstab["purity"]            = crosstab.drop(columns=["dominant_celltype"]).max(axis=1)

print("\nCluster → dominant cell type (purity %):")
for idx, row in crosstab[["dominant_celltype", "purity"]].iterrows():
    print(f"  Cluster {idx:>2}: {row['dominant_celltype']:<30} {row['purity']*100:.1f}%")

# Save to table
crosstab.to_csv("results/tables/A4_cluster_celltype_crosstab.csv")
print("\nSaved: results/tables/A4_cluster_celltype_crosstab.csv")


# ─────────────────────────────────────────────
# STEP 3 — Verify with canonical marker genes
# ─────────────────────────────────────────────
# Check that known cell type markers are expressed in the expected clusters.
# This confirms that Wu et al.'s annotations are consistent with gene expression.

MARKERS = {
    "Cancer Epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T cells":            ["CD3D", "CD3E", "CD8A", "CD4"],
    "B cells":            ["CD19", "MS4A1", "CD79A"],
    "Myeloid":            ["CD14", "CD68", "LYZ", "CST3"],
    "Fibroblasts":        ["COL1A1", "COL3A1", "FAP", "PDGFRA"],
    "Endothelial":        ["PECAM1", "VWF", "CDH5"],
    "Plasmablasts":       ["MZB1", "IGHG1", "XBP1"],
    "Normal Epithelial":  ["ALDH1A3", "SLPI"],
    "Mast cells":         ["TPSAB1", "CPA3"],
    "pDCs":               ["LILRA4", "IL3RA"],
}

# Filter to genes present in the dataset
present_markers = {}
for celltype, genes in MARKERS.items():
    avail = [g for g in genes if g in adata.var_names]
    if avail:
        present_markers[celltype] = avail

print(f"\nMarker genes available in dataset:")
for ct, genes in present_markers.items():
    print(f"  {ct}: {genes}")

# Dot plot: average expression of markers per cluster
all_marker_genes = [g for genes in present_markers.values() for g in genes]

fig_dp = sc.pl.dotplot(
    adata,
    var_names=present_markers,
    groupby="leiden_0.3",
    standard_scale="var",
    return_fig=True,
)
fig_dp.savefig("results/figures/A4_dotplot_markers.png", dpi=150, bbox_inches="tight")
plt.close()
print("\nSaved: results/figures/A4_dotplot_markers.png")


# ─────────────────────────────────────────────
# STEP 4 — Heatmap: cluster vs cell type composition
# ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 6))
heatmap_data = pd.crosstab(
    adata.obs["leiden_0.3"],
    adata.obs["celltype_major"],
    normalize="index",
) * 100  # convert to percent

sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".0f",
    cmap="Blues",
    linewidths=0.3,
    ax=ax,
    cbar_kws={"label": "% of cluster"},
)
ax.set_title("Cell type composition per Leiden cluster (res=0.3)\n% cells per cluster")
ax.set_xlabel("Cell Type (major)")
ax.set_ylabel("Leiden Cluster")
plt.tight_layout()
plt.savefig("results/figures/A4_cluster_composition_heatmap.png", dpi=150)
plt.close()
print("Saved: results/figures/A4_cluster_composition_heatmap.png")


# ─────────────────────────────────────────────
# STEP 5 — Annotate clusters on UMAP
# ─────────────────────────────────────────────
# Use Wu et al. annotations directly (they are expert-curated)
# and overlay on our UMAP to confirm spatial coherence.

fig, axes = plt.subplots(1, 3, figsize=(21, 6))

sc.pl.umap(adata, color="celltype_major",
           title="Cell Type Major (Wu et al. annotations)",
           legend_loc="right margin", legend_fontsize=8,
           ax=axes[0], show=False)

sc.pl.umap(adata, color="leiden_0.3",
           title="Leiden Clusters (res=0.3)",
           legend_loc="on data", legend_fontsize=9,
           ax=axes[1], show=False)

sc.pl.umap(adata, color="subtype",
           title="Breast Cancer Subtype",
           ax=axes[2], show=False)

plt.tight_layout()
plt.savefig("results/figures/A4_UMAP_annotations.png", dpi=150)
plt.close()
print("Saved: results/figures/A4_UMAP_annotations.png")


# ─────────────────────────────────────────────
# STEP 6 — Subset to Cancer Epithelial cells
# ─────────────────────────────────────────────
# All downstream CSC analysis (Phases A5, G1-G5) operates on Cancer Epithelial
# cells only. Normal epithelial cells and stromal/immune cells are excluded.
# The CSC population lives within the Cancer Epithelial compartment.

cancer_epi = adata[adata.obs["celltype_major"] == "Cancer Epithelial"].copy()
print(f"\nCancer Epithelial cells: {cancer_epi.n_obs}")
print(f"Subtypes represented:")
print(cancer_epi.obs["subtype"].value_counts().to_string())
print(f"\ncelltype_minor breakdown:")
print(cancer_epi.obs["celltype_minor"].value_counts().to_string())

# UMAP of Cancer Epithelial subset coloured by subtype and minor annotation
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc.pl.umap(cancer_epi, color="subtype",
           title="Cancer Epithelial — breast cancer subtype",
           ax=axes[0], show=False)

sc.pl.umap(cancer_epi, color="celltype_minor",
           title="Cancer Epithelial — minor annotation",
           legend_loc="right margin", legend_fontsize=7,
           ax=axes[1], show=False)

plt.tight_layout()
plt.savefig("results/figures/A4_cancer_epi_UMAP.png", dpi=150)
plt.close()
print("Saved: results/figures/A4_cancer_epi_UMAP.png")

# Save subsets
out_full = "data/processed/brca_A4_annotated.h5ad"
out_epi  = "data/processed/brca_A4_cancer_epi.h5ad"

adata.write_h5ad(out_full)
cancer_epi.write_h5ad(out_epi)

print(f"\nSaved: {out_full}  ({adata.n_obs} cells)")
print(f"Saved: {out_epi}  ({cancer_epi.n_obs} cells)")
print("\n✓ Phase A4 complete. Next: Phase A5 — CSC Identification & Stemness Scoring")
