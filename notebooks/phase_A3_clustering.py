# Phase A3 — Dimensionality Reduction & Clustering
# Input:  data/processed/brca_A2_preprocessed.h5ad
# Output: data/processed/brca_A3_clustered.h5ad
#
# Steps:
#   1. Scale data
#   2. PCA (50 components)
#   3. Harmony batch correction (per patient)
#   4. Nearest-neighbor graph
#   5. UMAP
#   6. Leiden clustering
#   7. Save

import scanpy as sc
import harmonypy as hm
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

sc.settings.verbosity = 1
sc.settings.figdir   = "results/figures/"
os.makedirs("results/figures", exist_ok=True)

# ─────────────────────────────────────────────
# STEP 1 — Load preprocessed data
# ─────────────────────────────────────────────
print("Loading preprocessed data...")
adata = sc.read_h5ad("data/processed/brca_A2_preprocessed.h5ad")
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Patients (batches): {adata.obs['orig.ident'].nunique()}")


# ─────────────────────────────────────────────
# STEP 2 — Scale (HVGs only)
# ─────────────────────────────────────────────
# Scaling centres each gene to mean=0, std=1 across cells.
# We only scale the 2,000 highly variable genes selected in A2.
# Clip at max_value=10 to prevent outlier genes from dominating PCA.

print("Scaling HVG expression...")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)


# ─────────────────────────────────────────────
# STEP 3 — PCA
# ─────────────────────────────────────────────
# PCA compresses 2,000 HVG dimensions → 50 principal components.
# Each PC captures a major axis of variation across cells.
# PC1 often captures cell cycle; PC2-20 typically capture cell type differences.

print("Running PCA (50 components)...")
sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")

# Copy PCA embedding and variance stats back to the main object
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.uns["pca"]    = adata_hvg.uns["pca"]

# Elbow plot — find how many PCs capture meaningful variance
sc.pl.pca_variance_ratio(adata_hvg, n_pcs=50, log=True,
                         save="_A3_elbow.png", show=False)
print("Saved: results/figures/pca_variance_ratio_A3_elbow.png")


# ─────────────────────────────────────────────
# STEP 4 — Harmony batch correction
# ─────────────────────────────────────────────
# This dataset has 26 tumor samples from multiple patients.
# Without correction, cells cluster by patient rather than cell type.
# Harmony adjusts the PCA embedding so patient-level technical variation
# is removed while biological variation (cell types, states) is preserved.

print("Running Harmony batch correction (by patient)...")
# Bypass scanpy's wrapper — it still transposes Z_corr which breaks with harmonypy>=2.0
# where Z_corr is already (n_cells, n_pcs), not (n_pcs, n_cells)
ho = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    vars_use=["orig.ident"],
    max_iter_harmony=20,
    random_state=42,
)
adata.obsm["X_pca_harmony"] = ho.Z_corr   # shape: (n_cells, n_pcs) — no .T needed
print("Harmony complete.")


# ─────────────────────────────────────────────
# STEP 5 — Nearest-Neighbor Graph
# ─────────────────────────────────────────────
# Build a k-nearest-neighbor graph in the Harmony-corrected PCA space.
# Each cell is connected to its 15 most similar cells.
# This graph is the basis for both UMAP layout and Leiden clustering.

print("Building neighbor graph...")
sc.pp.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=30,
    use_rep="X_pca_harmony",
    random_state=42,
)


# ─────────────────────────────────────────────
# STEP 6 — UMAP
# ─────────────────────────────────────────────
# UMAP projects the high-dimensional neighbor graph into 2D for visualization.
# Cells that are transcriptionally similar end up close together.
# Important: distances on UMAP are NOT meaningful — only neighborhoods are.

print("Computing UMAP...")
sc.tl.umap(adata, random_state=42)
print("UMAP complete.")


# ─────────────────────────────────────────────
# STEP 7 — Leiden Clustering
# ─────────────────────────────────────────────
# Leiden finds groups of densely connected cells in the neighbor graph.
# Resolution controls granularity: lower = fewer, bigger clusters.
# We run two resolutions to compare coarse vs fine clustering.

print("Running Leiden clustering...")
sc.tl.leiden(adata, resolution=0.3, key_added="leiden_0.3", random_state=42)
sc.tl.leiden(adata, resolution=0.6, key_added="leiden_0.6", random_state=42)

n_clusters_coarse = adata.obs["leiden_0.3"].nunique()
n_clusters_fine   = adata.obs["leiden_0.6"].nunique()
print(f"Leiden clusters (res=0.3): {n_clusters_coarse}")
print(f"Leiden clusters (res=0.6): {n_clusters_fine}")


# ─────────────────────────────────────────────
# STEP 8 — Visualize
# ─────────────────────────────────────────────
print("Generating UMAP plots...")

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Row 1: biological annotations
sc.pl.umap(adata, color="celltype_major", title="Cell Type (major)",
           legend_loc="on data", legend_fontsize=7,
           ax=axes[0, 0], show=False)

sc.pl.umap(adata, color="subtype", title="Cancer Subtype",
           ax=axes[0, 1], show=False)

sc.pl.umap(adata, color="orig.ident", title="Patient ID (batch)",
           legend_loc="right margin", legend_fontsize=6,
           ax=axes[0, 2], show=False)

# Row 2: cluster results
sc.pl.umap(adata, color="leiden_0.3",
           title=f"Leiden res=0.3 ({n_clusters_coarse} clusters)",
           legend_loc="on data", ax=axes[1, 0], show=False)

sc.pl.umap(adata, color="leiden_0.6",
           title=f"Leiden res=0.6 ({n_clusters_fine} clusters)",
           legend_loc="on data", ax=axes[1, 1], show=False)

# Known CSC marker expression on UMAP
sc.pl.umap(adata, color="CD44", title="CD44 (CSC marker)",
           color_map="Reds", ax=axes[1, 2], show=False)

plt.tight_layout()
plt.savefig("results/figures/A3_UMAP_overview.png", dpi=150)
plt.close()
print("Saved: results/figures/A3_UMAP_overview.png")

# Additional: cluster size table
cluster_sizes = adata.obs["leiden_0.3"].value_counts().sort_index()
print("\nCluster sizes (res=0.3):")
print(cluster_sizes.to_string())


# ─────────────────────────────────────────────
# STEP 9 — Save
# ─────────────────────────────────────────────
out_path = "data/processed/brca_A3_clustered.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved: {out_path}")
print("\n✓ Phase A3 complete. Next: Phase A4 — Cell Type Annotation")
