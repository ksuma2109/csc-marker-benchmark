# Phase A2 — Quality Control & Preprocessing
# Dataset: GSE176078 — Wu et al. 2021 Breast Cancer scRNA-seq (100,064 cells)
#
# Run this script from the stem_cells directory with the venv active:
#   source venv/bin/activate
#   python notebooks/phase_A2_QC.py

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os

sc.settings.verbosity = 1
sc.settings.figdir = "results/figures/"
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

DATA_DIR = "data/raw/Wu_etal_2021_BRCA_scRNASeq/"

# ─────────────────────────────────────────────
# STEP 1 — Load the count matrix
# ─────────────────────────────────────────────
# The .mtx file is stored genes × cells (rows=genes, cols=cells).
# We transpose it (.T) so rows=cells, cols=genes — the standard orientation.

print("Loading count matrix...")
adata = sc.read_mtx(DATA_DIR + "count_matrix_sparse.mtx").T

# Attach cell barcodes as row names
barcodes = pd.read_csv(DATA_DIR + "count_matrix_barcodes.tsv", header=None)[0]
adata.obs_names = barcodes.values

# Attach gene names as column names
genes = pd.read_csv(DATA_DIR + "count_matrix_genes.tsv", header=None)[0]
adata.var_names = genes.values

print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
# Expected: 100064 cells × 29733 genes


# ─────────────────────────────────────────────
# STEP 2 — Attach metadata
# ─────────────────────────────────────────────
# The metadata CSV has one row per cell with cell type annotations,
# QC metrics (nCount_RNA, nFeature_RNA, percent.mito), and cancer subtype.

print("Attaching metadata...")
meta = pd.read_csv(DATA_DIR + "metadata.csv", index_col=0)
adata.obs = meta

print("Metadata columns:", list(adata.obs.columns))
print(adata.obs.head(3))


# ─────────────────────────────────────────────
# STEP 3 — Inspect QC metrics
# ─────────────────────────────────────────────
# The authors already ran QC before publishing, so these metrics reflect
# the cleaned dataset. We inspect them to understand the data quality.

# nCount_RNA    = total RNA molecules detected per cell (library size)
# nFeature_RNA  = number of unique genes detected per cell
# percent.mito  = % of reads from mitochondrial genes (high = dying cell)

print("\n--- QC metric summary ---")
print(adata.obs[["nCount_RNA", "nFeature_RNA", "percent.mito"]].describe().round(2))


# ─────────────────────────────────────────────
# STEP 4 — Visualize QC metrics
# ─────────────────────────────────────────────
# Violin plots show the distribution of QC metrics across all cells.
# This tells you whether the dataset looks clean.

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

axes[0].violinplot(adata.obs["nCount_RNA"], showmedians=True)
axes[0].set_title("Total RNA per cell\n(nCount_RNA)")
axes[0].set_ylabel("UMI count")

axes[1].violinplot(adata.obs["nFeature_RNA"], showmedians=True)
axes[1].set_title("Genes detected per cell\n(nFeature_RNA)")
axes[1].set_ylabel("Gene count")

axes[2].violinplot(adata.obs["percent.mito"], showmedians=True)
axes[2].set_title("Mitochondrial %\n(percent.mito)")
axes[2].set_ylabel("% mito reads")

plt.tight_layout()
plt.savefig("results/figures/A2_QC_violin.png", dpi=150)
plt.show()
print("Saved: results/figures/A2_QC_violin.png")


# ─────────────────────────────────────────────
# STEP 5 — Inspect cell type composition
# ─────────────────────────────────────────────
# See how many cells of each major cell type are in the dataset.
# We need to identify the Cancer Epithelial cells — these are the ones
# we will subset for CSC analysis in later phases.

print("\n--- Cell type composition (celltype_major) ---")
print(adata.obs["celltype_major"].value_counts())

print("\n--- Breast cancer subtype breakdown ---")
print(adata.obs["subtype"].value_counts())

# Plot cell type composition as a bar chart
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

adata.obs["celltype_major"].value_counts().plot(
    kind="barh", ax=axes[0], color="steelblue"
)
axes[0].set_title("Cell type composition")
axes[0].set_xlabel("Number of cells")

adata.obs["subtype"].value_counts().plot(
    kind="bar", ax=axes[1], color="coral"
)
axes[1].set_title("Cancer subtype breakdown")
axes[1].set_ylabel("Number of cells")
axes[1].tick_params(axis="x", rotation=30)

plt.tight_layout()
plt.savefig("results/figures/A2_cell_composition.png", dpi=150)
plt.show()
print("Saved: results/figures/A2_cell_composition.png")


# ─────────────────────────────────────────────
# STEP 6 — Normalize counts
# ─────────────────────────────────────────────
# Different cells capture different amounts of RNA (library size).
# Normalization makes cells comparable by scaling each cell to a
# total count of 10,000, then log-transforming.
#
# After log1p:  value = log(count + 1)
#   count=0  → 0      (silent genes stay 0)
#   count=1  → 0.69   (low expression)
#   count=10 → 2.40   (moderate)
#   count=100 → 4.62  (high expression)

print("\nNormalizing...")
sc.pp.normalize_total(adata, target_sum=1e4)   # scale each cell to 10,000 total
sc.pp.log1p(adata)                              # log-transform
print("Normalization complete.")


# ─────────────────────────────────────────────
# STEP 7 — Identify highly variable genes
# ─────────────────────────────────────────────
# Most of the 29,733 genes are either silent or constitutively expressed
# across all cells (housekeeping genes) — they carry no useful information
# for distinguishing cell types or states.
#
# We select the top 2,000 highly variable genes (HVGs) — genes whose
# expression varies most across cells. These are the informative genes
# used for downstream analysis.

print("\nFinding highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", span=1.0)

n_hvg = adata.var["highly_variable"].sum()
print(f"Highly variable genes selected: {n_hvg}")

sc.pl.highly_variable_genes(adata, save="_A2_HVG.png", show=False)
print("Saved: results/figures/highly_variable_genes_A2_HVG.png")


# ─────────────────────────────────────────────
# STEP 8 — Save processed data
# ─────────────────────────────────────────────
# Save the processed AnnData object so we don't have to redo
# loading + normalization in every subsequent phase.
# .h5ad is the standard single-cell data format (HDF5-based).

out_path = "data/processed/brca_A2_preprocessed.h5ad"
os.makedirs("data/processed", exist_ok=True)
adata.write_h5ad(out_path)
print(f"\nSaved processed data → {out_path}")

print("\n✓ Phase A2 complete. Next: Phase A3 — Dimensionality Reduction & Clustering")
print("  Load the saved file with: adata = sc.read_h5ad('data/processed/brca_A2_preprocessed.h5ad')")
