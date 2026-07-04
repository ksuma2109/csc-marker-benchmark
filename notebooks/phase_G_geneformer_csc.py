# Stage 2 — Geneformer: Attention-Based CSC Gene Program Discovery
# Dataset: GSE176078 — Wu et al. 2021 Breast Cancer (Cancer Epithelial cells)
#
# Pipeline:
#   G1 — Generate CSC pseudo-labels from stemness scoring
#   G2 — Tokenize data for Geneformer (gene symbol → Ensembl → tokens)
#   G3 — Fine-tune Geneformer (CSC-high vs CSC-low)
#   G4 — Extract attention weights → gene importance ranking
#   G5 — Compare against Stage 1 markers
#
# Run from stem_cells/ with venv active:
#   source venv/bin/activate
#   python notebooks/phase_G_geneformer_csc.py

import os
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import mygene
import loompy
import torch
import pickle
from pathlib import Path
from datasets import Dataset
from geneformer import TranscriptomeTokenizer, Classifier
from transformers import BertForSequenceClassification, BertConfig
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)
os.makedirs("data/geneformer", exist_ok=True)

DEVICE = (
    "mps" if torch.backends.mps.is_available()
    else "cuda" if torch.cuda.is_available()
    else "cpu"
)
print(f"Using device: {DEVICE}")


# ─────────────────────────────────────────────
# PHASE G1 — Generate CSC Pseudo-Labels
# ─────────────────────────────────────────────
# We don't have ground-truth CSC labels. Instead we use the Malta et al.
# mRNAsi stemness signature (a gene list derived from pluripotent stem cells)
# to score each Cancer Epithelial cell. Top 25% = CSC-high, bottom 25% = CSC-low.
# These become the training labels for Geneformer fine-tuning.

print("\n=== Phase G1: Generating CSC pseudo-labels ===")

# Load preprocessed data from Phase A2
adata = sc.read_h5ad("data/processed/brca_A2_preprocessed.h5ad")
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Subset to Cancer Epithelial cells only
cancer_epi = adata[adata.obs["celltype_major"] == "Cancer Epithelial"].copy()
print(f"Cancer Epithelial cells: {cancer_epi.n_obs}")

# Score stemness using decoupler with the PanglaoDB/MSigDB stem cell signature
# We use a curated set of core pluripotency/stemness TF target genes
# (proxy for mRNAsi; full mRNAsi requires the original OCLR model weights)
STEMNESS_GENES = [
    "SOX2", "OCT4", "NANOG", "KLF4", "MYC", "LIN28A", "LIN28B",
    "ALDH1A1", "CD44", "PROM1", "ITGA6", "EPCAM", "SSEA4",
    "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "VIM",
    "BMI1", "EZH2", "HMGA2", "FOXM1", "YAP1", "NOTCH1",
    "WNT5A", "FZD7", "LGR5", "AXIN2", "TCF4",
]

# Filter to genes present in the dataset
available = [g for g in STEMNESS_GENES if g in cancer_epi.var_names]
print(f"Stemness genes available in dataset: {len(available)}/{len(STEMNESS_GENES)}")

# Score each cell: mean normalized expression of available stemness genes
cancer_epi.obs["stemness_score"] = (
    cancer_epi[:, available].X.toarray().mean(axis=1)
    if hasattr(cancer_epi.X, "toarray")
    else cancer_epi[:, available].X.mean(axis=1)
)

# Assign labels: top 25% = csc_high (1), bottom 25% = csc_low (0)
q75 = cancer_epi.obs["stemness_score"].quantile(0.75)
q25 = cancer_epi.obs["stemness_score"].quantile(0.25)

cancer_epi.obs["csc_label"] = "middle"
cancer_epi.obs.loc[cancer_epi.obs["stemness_score"] >= q75, "csc_label"] = "csc_high"
cancer_epi.obs.loc[cancer_epi.obs["stemness_score"] <= q25, "csc_label"] = "csc_low"

# Subset to top/bottom quartile only (discard middle 50%)
labelled = cancer_epi[cancer_epi.obs["csc_label"] != "middle"].copy()
label_counts = labelled.obs["csc_label"].value_counts()
print(f"Labelled cells: {label_counts.to_dict()}")

# Plot stemness score distribution
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(cancer_epi.obs["stemness_score"], bins=80, color="steelblue", alpha=0.7)
ax.axvline(q75, color="red",    linestyle="--", label=f"Q75 (CSC-high cutoff): {q75:.3f}")
ax.axvline(q25, color="orange", linestyle="--", label=f"Q25 (CSC-low cutoff): {q25:.3f}")
ax.set_xlabel("Stemness Score")
ax.set_ylabel("Number of cells")
ax.set_title("Cancer Epithelial cells — stemness score distribution")
ax.legend()
plt.tight_layout()
plt.savefig("results/figures/G1_stemness_distribution.png", dpi=150)
plt.close()
print("Saved: results/figures/G1_stemness_distribution.png")

# Save label mapping for comparison in Phase G5
labelled.obs[["stemness_score", "csc_label", "subtype"]].to_csv(
    "results/tables/G1_csc_labels.csv"
)
print("Saved: results/tables/G1_csc_labels.csv")


# ─────────────────────────────────────────────
# PHASE G2 — Tokenize for Geneformer
# ─────────────────────────────────────────────
# Geneformer requires:
# 1. Gene names as Ensembl IDs (ENSG...) — not gene symbols
# 2. Raw (un-normalized) integer counts — it does its own rank-based encoding
# 3. Data in .loom format with specific metadata fields

print("\n=== Phase G2: Tokenizing for Geneformer ===")

# Step 2a — Convert gene symbols → Ensembl IDs
print("Querying Ensembl IDs via mygene...")
mg = mygene.MyGeneInfo()
gene_symbols = labelled.var_names.tolist()

# Query in batches (mygene has rate limits)
results = mg.querymany(
    gene_symbols,
    scopes="symbol",
    fields="ensembl.gene",
    species="human",
    returnall=False,
)

# Build symbol → Ensembl ID mapping
symbol_to_ensembl = {}
for r in results:
    if "ensembl" in r:
        ens = r["ensembl"]
        if isinstance(ens, list):
            ens = ens[0]
        symbol_to_ensembl[r["query"]] = ens.get("gene", None)
    elif "notfound" not in r:
        symbol_to_ensembl[r["query"]] = None

n_mapped = sum(1 for v in symbol_to_ensembl.values() if v is not None)
print(f"Genes mapped to Ensembl IDs: {n_mapped}/{len(gene_symbols)}")

# Add Ensembl IDs to var
labelled.var["ensembl_id"] = [symbol_to_ensembl.get(g) for g in labelled.var_names]

# Keep only genes with valid Ensembl IDs
has_ensembl = labelled.var["ensembl_id"].notna()
labelled_ensembl = labelled[:, has_ensembl].copy()
print(f"Genes with Ensembl ID: {labelled_ensembl.n_vars}")

# Step 2b — Write .loom file for Geneformer tokenizer
# Geneformer expects:
#   - Row attribute "ensembl_id" on genes
#   - Column attributes "n_counts" (total UMIs per cell) and "cell_type" on cells
print("Writing loom file...")

# Need raw counts — reload from the raw count matrix
raw = sc.read_mtx("data/raw/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx").T
raw.obs_names = pd.read_csv(
    "data/raw/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv", header=None
)[0].values
raw.var_names = pd.read_csv(
    "data/raw/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv", header=None
)[0].values

# Subset to labelled cells and Ensembl-mapped genes
raw_subset = raw[labelled_ensembl.obs_names, labelled_ensembl.var_names].copy()
raw_subset.obs["csc_label"] = labelled_ensembl.obs["csc_label"].values
raw_subset.obs["n_counts"] = np.asarray(raw_subset.X.sum(axis=1)).flatten()
raw_subset.var["ensembl_id"] = labelled_ensembl.var["ensembl_id"].values

# Write loom
loom_path = "data/geneformer/brca_cancer_epi_labelled.loom"
row_attrs = {"ensembl_id": raw_subset.var["ensembl_id"].values}
col_attrs = {
    "n_counts":  raw_subset.obs["n_counts"].values.astype(int),
    "csc_label": raw_subset.obs["csc_label"].values,
    "cell_type": raw_subset.obs["csc_label"].values,
}
matrix = raw_subset.X.T.toarray() if hasattr(raw_subset.X, "toarray") else raw_subset.X.T
loompy.create(loom_path, matrix, row_attrs, col_attrs)
print(f"Saved: {loom_path}")

# Step 2c — Tokenize with Geneformer
print("Tokenizing with Geneformer TranscriptomeTokenizer...")
tokenizer = TranscriptomeTokenizer(
    custom_attr_name_dict={"cell_type": "cell_type"},
    nproc=4,
)
tokenizer.tokenize_data(
    loom_file_paths=[loom_path],
    output_directory="data/geneformer/",
    output_prefix="brca_csc_tokenized",
    file_format="loom",
)
print("Tokenized dataset saved to data/geneformer/")


# ─────────────────────────────────────────────
# PHASE G3 — Fine-Tune Geneformer
# ─────────────────────────────────────────────
# Fine-tune the pre-trained Geneformer model to classify CSC-high vs CSC-low.
# After fine-tuning, the attention heads will have learned to focus on genes
# most predictive of CSC identity.

print("\n=== Phase G3: Fine-tuning Geneformer ===")

# Download pre-trained weights (6-layer Geneformer, ~30M param)
from huggingface_hub import snapshot_download
model_dir = "geneformer_repo"  # already downloaded

# Map string labels → integers
label2id = {"csc_low": 0, "csc_high": 1}
id2label = {0: "csc_low", 1: "csc_high"}

classifier = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "cell_type", "states": "all"},
    filter_data=None,
    training_args={
        "num_train_epochs": 5,
        "learning_rate": 5e-5,
        "per_device_train_batch_size": 12,
        "per_device_eval_batch_size": 12,
        "warmup_steps": 500,
        "weight_decay": 0.01,
        "output_dir": "results/geneformer_model",
        "evaluation_strategy": "epoch",
        "save_strategy": "epoch",
        "load_best_model_at_end": True,
    },
    max_ncells=None,
    freeze_layers=4,        # freeze first 4 layers, fine-tune top 2
    num_crossval_splits=1,
    forward_batch_size=100,
    nproc=4,
)

classifier.prepare_data(
    input_data_file="data/geneformer/brca_csc_tokenized.dataset",
    output_directory="data/geneformer/",
    output_prefix="brca_csc",
    split_id_dict={"attr_key": "cell_type", "train": 0.8, "test": 0.2},
)

all_metrics = classifier.validate(
    model_directory=model_dir,
    prepared_input_data_file="data/geneformer/brca_csc_prep_data.dataset",
    id_class_dict_file="data/geneformer/brca_csc_id_class_dict.pkl",
    output_directory="results/geneformer_model/",
    output_prefix="brca_csc",
    predict_trainer=False,
)
print("Fine-tuning metrics:", all_metrics)


# ─────────────────────────────────────────────
# PHASE G4 — Extract Attention Weights → Gene Ranking
# ─────────────────────────────────────────────
# Forward-pass CSC-high cells through the fine-tuned model with
# output_attentions=True. Aggregate attention across all heads and
# layers to get a per-gene importance score.

print("\n=== Phase G4: Extracting attention weights ===")

from transformers import AutoTokenizer, BertModel
from datasets import load_from_disk
import pickle

# Load fine-tuned model
fine_tuned_model = BertForSequenceClassification.from_pretrained(
    "results/geneformer_model"
)
fine_tuned_model.eval()
fine_tuned_model.to(DEVICE)

# Load token → gene dictionary from Geneformer
with open(os.path.join(model_dir, "geneformer/token_dictionary.pkl"), "rb") as f:
    token_dict = pickle.load(f)
id_to_gene = {v: k for k, v in token_dict.items()}

# Load tokenized dataset and filter to CSC-high cells
dataset = load_from_disk("data/geneformer/brca_csc_tokenized.dataset")
csc_high_dataset = dataset.filter(lambda x: x["cell_type"] == "csc_high")
print(f"CSC-high cells for attention extraction: {len(csc_high_dataset)}")

# Accumulate attention per gene token across cells
gene_attention_sum = {}
gene_attention_count = {}

BATCH_SIZE = 16
for i in range(0, min(len(csc_high_dataset), 500), BATCH_SIZE):
    batch = csc_high_dataset[i : i + BATCH_SIZE]
    input_ids = torch.tensor(batch["input_ids"]).to(DEVICE)
    attention_mask = (input_ids != 0).long()

    with torch.no_grad():
        outputs = fine_tuned_model(
            input_ids=input_ids,
            attention_mask=attention_mask,
            output_attentions=True,
        )

    # outputs.attentions: tuple of [n_layers] each (batch, heads, seq, seq)
    # Mean over layers and heads → (batch, seq, seq)
    stacked = torch.stack(outputs.attentions, dim=0)      # (layers, batch, heads, seq, seq)
    mean_attn = stacked.mean(dim=(0, 2))                  # (batch, seq, seq)

    # Column-sum: how much attention each gene RECEIVES from all others
    gene_importance = mean_attn.sum(dim=1)                # (batch, seq)

    for b_idx in range(gene_importance.shape[0]):
        tokens = input_ids[b_idx].cpu().tolist()
        scores = gene_importance[b_idx].cpu().tolist()
        for token_id, score in zip(tokens, scores):
            if token_id in id_to_gene:
                gene = id_to_gene[token_id]
                gene_attention_sum[gene]   = gene_attention_sum.get(gene, 0) + score
                gene_attention_count[gene] = gene_attention_count.get(gene, 0) + 1

# Average attention score per gene
gene_avg_attention = {
    g: gene_attention_sum[g] / gene_attention_count[g]
    for g in gene_attention_sum
}

# Sort by attention score descending
attention_df = (
    pd.DataFrame.from_dict(gene_avg_attention, orient="index", columns=["attention_score"])
    .sort_values("attention_score", ascending=False)
    .reset_index()
    .rename(columns={"index": "ensembl_id"})
)

# Map Ensembl IDs back to gene symbols
ensembl_to_symbol = {v: k for k, v in symbol_to_ensembl.items() if v}
attention_df["gene_symbol"] = attention_df["ensembl_id"].map(ensembl_to_symbol)

# Save
attention_df.to_csv("results/tables/G4_geneformer_gene_ranking.csv", index=False)
print(f"\nTop 20 genes by Geneformer attention:")
print(attention_df.head(20)[["gene_symbol", "ensembl_id", "attention_score"]].to_string())
print("Saved: results/tables/G4_geneformer_gene_ranking.csv")

# Plot top 30 genes
top30 = attention_df.dropna(subset=["gene_symbol"]).head(30)
fig, ax = plt.subplots(figsize=(10, 8))
ax.barh(top30["gene_symbol"][::-1], top30["attention_score"][::-1], color="steelblue")
ax.set_xlabel("Mean Attention Score")
ax.set_title("Top 30 CSC genes by Geneformer attention\n(Stage 2 — data-driven)")
plt.tight_layout()
plt.savefig("results/figures/G4_top_genes_attention.png", dpi=150)
plt.close()
print("Saved: results/figures/G4_top_genes_attention.png")


# ─────────────────────────────────────────────
# PHASE G5 — Compare Stage 1 vs Stage 2
# ─────────────────────────────────────────────
# Stage 1 markers come from differential expression (Phase A5).
# Stage 2 markers come from Geneformer attention (Phase G4).
# This phase computes the overlap and identifies novel candidates.

print("\n=== Phase G5: Comparing Stage 1 vs Stage 2 markers ===")

# Load Stage 1 marker genes (produced in Phase A5)
stage1_path = "results/tables/A5_csc_markers_DE.csv"
if not os.path.exists(stage1_path):
    print(f"Stage 1 markers not found at {stage1_path}")
    print("Run Phase A5 first, then re-run Phase G5.")
else:
    stage1_df = pd.read_csv(stage1_path)
    stage1_genes = set(stage1_df["gene_symbol"].dropna().head(200))

    # Stage 2 top 200 by attention
    stage2_genes = set(
        attention_df.dropna(subset=["gene_symbol"])["gene_symbol"].head(200)
    )

    shared      = stage1_genes & stage2_genes
    stage2_only = stage2_genes - stage1_genes
    stage1_only = stage1_genes - stage2_genes

    jaccard = len(shared) / len(stage1_genes | stage2_genes)

    print(f"\n--- Comparison (top 200 genes each) ---")
    print(f"Stage 1 (DE)           : {len(stage1_genes)} genes")
    print(f"Stage 2 (Attention)    : {len(stage2_genes)} genes")
    print(f"Shared                 : {len(shared)} genes")
    print(f"Stage 2 only (novel)   : {len(stage2_only)} genes")
    print(f"Stage 1 only           : {len(stage1_only)} genes")
    print(f"Jaccard similarity     : {jaccard:.3f}")

    print(f"\nTop Stage-2-only candidates (novel CSC markers):")
    # Rank stage2_only by attention score
    novel = (
        attention_df[attention_df["gene_symbol"].isin(stage2_only)]
        .dropna(subset=["gene_symbol"])
        .head(20)
    )
    print(novel[["gene_symbol", "attention_score"]].to_string())

    # Save comparison tables
    pd.DataFrame(sorted(shared),      columns=["gene"]).to_csv("results/tables/G5_shared_markers.csv",      index=False)
    pd.DataFrame(sorted(stage2_only), columns=["gene"]).to_csv("results/tables/G5_novel_stage2_markers.csv", index=False)
    pd.DataFrame(sorted(stage1_only), columns=["gene"]).to_csv("results/tables/G5_stage1_only_markers.csv",  index=False)

    # Venn-style bar chart
    fig, ax = plt.subplots(figsize=(7, 4))
    categories = ["Stage 1 only\n(DE)", "Shared\n(both)", "Stage 2 only\n(Attention)"]
    counts     = [len(stage1_only), len(shared), len(stage2_only)]
    colors     = ["#FF6B6B", "#4ECDC4", "#45B7D1"]
    ax.bar(categories, counts, color=colors, edgecolor="white", linewidth=1.5)
    for i, (cat, cnt) in enumerate(zip(categories, counts)):
        ax.text(i, cnt + 1, str(cnt), ha="center", fontweight="bold")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Stage 1 vs Stage 2 CSC marker overlap\n(Jaccard = {jaccard:.3f})")
    plt.tight_layout()
    plt.savefig("results/figures/G5_marker_comparison.png", dpi=150)
    plt.close()
    print("Saved: results/figures/G5_marker_comparison.png")

print("\n✓ Geneformer pipeline complete.")
print("  Key outputs:")
print("  results/tables/G4_geneformer_gene_ranking.csv  — all genes ranked by attention")
print("  results/tables/G5_novel_stage2_markers.csv     — novel CSC candidates")
print("  results/figures/G4_top_genes_attention.png     — attention bar chart")
print("  results/figures/G5_marker_comparison.png       — Stage 1 vs Stage 2 overlap")
