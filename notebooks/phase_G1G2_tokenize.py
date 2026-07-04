# Phase G1 + G2 — CSC Pseudo-Labels → Geneformer Tokenization
#
# G1: Pull CSC labels from Phase A5 (stemness-scored Cancer Epithelial cells)
#     and merge with raw integer counts (all 29,733 genes).
#
# G2: Run TranscriptomeTokenizer to convert each cell into a ranked gene
#     sequence (Geneformer input format). Output: HuggingFace .dataset
#
# Input:  data/processed/brca_A5_csc_scored.h5ad   ← labels
#         data/raw/Wu_etal_2021_BRCA_scRNASeq/      ← raw counts
# Output: data/geneformer/input/brca_csc_labelled.h5ad
#         data/geneformer/brca_csc_tokenized.dataset

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("data/geneformer/input", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

sc.settings.verbosity = 1

DATA_DIR = "data/raw/Wu_etal_2021_BRCA_scRNASeq/"

# ─────────────────────────────────────────────
# PHASE G1 — Load labels from A5
# ─────────────────────────────────────────────
print("=== Phase G1: Loading CSC pseudo-labels from Phase A5 ===")

scored = sc.read_h5ad("data/processed/brca_A5_csc_scored.h5ad")
label_df = scored.obs[["csc_label", "stemness_composite", "subtype",
                        "celltype_minor", "orig.ident"]].copy()

label_counts = label_df["csc_label"].value_counts()
print(f"Label counts:\n{label_counts.to_string()}")

# Keep only CSC-high and CSC-low (drop middle 50%)
labelled_barcodes = label_df[label_df["csc_label"] != "middle"].index
print(f"\nLabelled cells (csc_high + csc_low): {len(labelled_barcodes)}")
print(f"  csc_high: {(label_df['csc_label'] == 'csc_high').sum()}")
print(f"  csc_low:  {(label_df['csc_label'] == 'csc_low').sum()}")

# ─────────────────────────────────────────────
# PHASE G2a — Load raw counts (all genes, integer)
# ─────────────────────────────────────────────
# Geneformer requires:
#   1. Raw un-normalized integer counts (X)
#   2. ALL genes (no feature selection) — it ranks them internally
#   3. var["ensembl_id"] — gene symbols work; tokenizer maps them to Ensembl
#   4. obs["n_counts"] — total UMI per cell

print("\n=== Phase G2a: Loading raw counts ===")
raw = sc.read_mtx(DATA_DIR + "count_matrix_sparse.mtx").T
raw.obs_names = pd.read_csv(DATA_DIR + "count_matrix_barcodes.tsv", header=None)[0].values
raw.var_names = pd.read_csv(DATA_DIR + "count_matrix_genes.tsv", header=None)[0].values
print(f"Raw matrix: {raw.n_obs} cells × {raw.n_vars} genes")

# Subset to labelled Cancer Epithelial cells only
raw_labelled = raw[labelled_barcodes, :].copy()
print(f"After subsetting to labelled cells: {raw_labelled.n_obs} cells")

# Add required Geneformer fields
raw_labelled.var["ensembl_id"] = raw_labelled.var_names.values   # gene symbols
raw_labelled.obs["n_counts"]   = np.asarray(raw_labelled.X.sum(axis=1)).flatten().astype(int)

# Merge CSC labels into obs
raw_labelled.obs = raw_labelled.obs.join(
    label_df[["csc_label", "stemness_composite", "subtype", "celltype_minor"]],
    how="left"
)

print(f"\nObs columns: {list(raw_labelled.obs.columns)}")
print(f"n_counts range: {raw_labelled.obs['n_counts'].min()} – {raw_labelled.obs['n_counts'].max()}")
print(f"Labels: {raw_labelled.obs['csc_label'].value_counts().to_dict()}")

# Filter out cells with zero counts (safety check)
raw_labelled = raw_labelled[raw_labelled.obs["n_counts"] > 0].copy()
print(f"After removing zero-count cells: {raw_labelled.n_obs}")

# ─────────────────────────────────────────────
# PHASE G2b — Save h5ad for tokenizer
# ─────────────────────────────────────────────
# Geneformer tokenize_data() reads from a DIRECTORY, so we save one file there.
h5ad_path = "data/geneformer/input/brca_csc_labelled.h5ad"
raw_labelled.write_h5ad(h5ad_path)
print(f"\nSaved: {h5ad_path}  ({raw_labelled.n_obs} cells × {raw_labelled.n_vars} genes)")

# Plot n_counts distribution by label
fig, ax = plt.subplots(figsize=(8, 4))
for label, color in [("csc_high", "#d62728"), ("csc_low", "#1f77b4")]:
    subset = raw_labelled.obs.loc[raw_labelled.obs["csc_label"] == label, "n_counts"]
    ax.hist(subset, bins=60, alpha=0.6, color=color, label=label, density=True)
ax.set_xlabel("Total UMI count (n_counts)")
ax.set_ylabel("Density")
ax.set_title("UMI count distribution by CSC label (Geneformer input)")
ax.legend()
plt.tight_layout()
plt.savefig("results/figures/G1_ncounts_distribution.png", dpi=150)
plt.close()
print("Saved: results/figures/G1_ncounts_distribution.png")

# ─────────────────────────────────────────────
# PHASE G2c — Tokenize with Geneformer
# ─────────────────────────────────────────────
# TranscriptomeTokenizer converts each cell's raw counts into a ranked
# sequence of gene tokens (highest expressed genes first, up to 4096).
# custom_attr_name_dict maps obs columns to dataset column names.

print("\n=== Phase G2c: Tokenizing with Geneformer TranscriptomeTokenizer ===")
from geneformer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer(
    custom_attr_name_dict={
        "csc_label":       "csc_label",
        "subtype":         "subtype",
        "celltype_minor":  "celltype_minor",
    },
    nproc=4,
    model_version="V2",
)

tk.tokenize_data(
    data_directory="data/geneformer/input/",
    output_directory="data/geneformer/",
    output_prefix="brca_csc",
    file_format="h5ad",
)
print("Tokenization complete.")

# Inspect the resulting dataset
from datasets import load_from_disk
dataset = load_from_disk("data/geneformer/brca_csc.dataset")
print(f"\nTokenized dataset: {len(dataset)} cells")
print(f"Columns: {dataset.column_names}")
print(f"Sample input_ids length: {len(dataset[0]['input_ids'])}")
print(f"Sample csc_label: {dataset[0]['csc_label']}")

# Label distribution in tokenized dataset
import collections
labels = [d["csc_label"] for d in dataset]
print(f"\nLabel distribution: {collections.Counter(labels)}")

print("\n✓ G1 + G2 complete.")
print("  Output: data/geneformer/brca_csc.dataset")
print("  Next: Phase G3 — Fine-tune Geneformer on CSC-high vs CSC-low")
