# Phase B1 — Cross-Cancer Validation: GBM (Darmanis et al. 2017, GSE84465)
#
# Goal: Run Stage 1 CSC pipeline on glioblastoma single-cell data,
#       then compare GBM-CSC markers with breast CSC markers to identify
#       a pan-cancer CSC gene program.
#
# Dataset: GSE84465 — 3,589 cells from 4 GBM patients
#          Darmanis et al. 2017 Cell Reports
#
# Pipeline mirrors breast cancer Phase A2–A5:
#   B1a: Download & load
#   B1b: QC + normalization
#   B1c: Clustering + annotation
#   B1d: GSC scoring + labelling
#   B1e: DE (Stage 1 for GBM)
#   B1f: Cross-cancer comparison
#
# Output: results/tables/B1_gbm_csc_markers.csv
#         results/figures/B1_*.png
#         results/tables/B1_pancancer_consensus.csv

import os, gc, requests, gzip, io
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

os.makedirs("data/raw/GSE84465",      exist_ok=True)
os.makedirs("data/processed",         exist_ok=True)
os.makedirs("results/tables",         exist_ok=True)
os.makedirs("results/figures",        exist_ok=True)

# ─────────────────────────────────────────────────────────────────
# B1a — DOWNLOAD GBM DATA FROM GEO
# ─────────────────────────────────────────────────────────────────
print("=" * 60)
print("B1a — DOWNLOADING GSE84465 (GBM)")
print("=" * 60)

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl"
FILES = {
    "GSE84465_GBM_All_data.csv.gz": "data/raw/GSE84465/counts.csv.gz",
    # No separate metadata file on GEO for this dataset
}

for fname, dest in FILES.items():
    if os.path.exists(dest):
        print(f"  Already downloaded: {dest}")
        continue
    url = f"{BASE_URL}/{fname}"
    print(f"  Downloading {fname} ...")
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"  Saved: {dest}")

# ─────────────────────────────────────────────────────────────────
# B1b — LOAD + QC + NORMALIZE
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B1b — LOADING AND PREPROCESSING")
print("=" * 60)

H5AD_PATH = "data/processed/gbm_B1_preprocessed.h5ad"
if os.path.exists(H5AD_PATH):
    print(f"  Loading cached: {H5AD_PATH}")
    adata = sc.read_h5ad(H5AD_PATH)
else:
    print("  Reading count matrix (space-delimited)...")
    with gzip.open("data/raw/GSE84465/counts.csv.gz", "rt") as f:
        counts = pd.read_csv(f, index_col=0, sep=r'\s+', engine="python")  # genes × cells
    # Strip surrounding quotes from gene names and cell IDs (GEO format artifact)
    counts.index   = counts.index.str.strip('"')
    counts.columns = counts.columns.str.strip('"')

    print(f"  Raw shape (genes × cells): {counts.shape}")
    counts = counts.T  # cells × genes
    adata  = sc.AnnData(X=counts.values.astype("float32"),
                        obs=pd.DataFrame(index=counts.index),
                        var=pd.DataFrame(index=counts.columns))

    # Darmanis 2017 cell-type annotations are embedded in cell IDs:
    # format "patientID.CellType_barcode" (e.g. "1001000173.G8")
    # The GEO series only has the count matrix — cell types from paper Table S1
    # We use the first component of the barcode (patient ID) as batch info
    print(f"  Cell IDs (first 5): {list(adata.obs.index[:5])}")
    print(f"  No separate metadata file — will use score-based annotation")

    # QC
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    mito_genes = adata.var_names.str.upper().str.startswith("MT-")
    mito_X = adata[:, mito_genes].X
    mito_sum = np.asarray(mito_X.sum(axis=1)).flatten()
    adata.obs["pct_mito"] = mito_sum / adata.obs["total_counts"] * 100
    before = adata.n_obs
    adata = adata[
        (adata.obs["total_counts"] >= 500) &
        (adata.obs["n_genes_by_counts"] >= 300) &
        (adata.obs["pct_mito"] < 20)
    ].copy()
    print(f"  After QC: {adata.n_obs} cells (removed {before - adata.n_obs})")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    # HVG → PCA → neighbours → UMAP → Leiden
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                 subset=False, layer=None)
    sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, resolution=0.4, key_added="leiden", random_state=42)

    adata.write_h5ad(H5AD_PATH)
    print(f"  Saved: {H5AD_PATH}")

print(f"  Dataset: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")
print(f"  Metadata columns: {list(adata.obs.columns)}")

# ─────────────────────────────────────────────────────────────────
# B1c — CELL TYPE ANNOTATION
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B1c — CELL TYPE ANNOTATION")
print("=" * 60)

# GSE84465 has no separate metadata file on GEO — only the count matrix.
# The dataset contains mixed cell types (tumor, astrocyte, microglia, neuron, etc.)
# We score all cells for GSC vs differentiated identity to find stem-like cells.
# Reference: Darmanis et al. 2017 Cell Reports — ~3,500 cells, 4 GBM patients
print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")
print(f"  All {adata.n_obs} cells — scoring for GSC identity across full dataset")
adata_tumor = adata  # use all cells; GSC scoring separates stem-like subset

# ─────────────────────────────────────────────────────────────────
# B1d — GSC SCORING + LABELLING
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B1d — GSC SCORING")
print("=" * 60)

# Published GBM CSC (GSC) marker signatures
# Source: Lathia et al. 2015 Genes & Dev; Neftel et al. 2019 Cell
GSC_SIGNATURE = ["SOX2","NES","PROM1","CD44","ALDH1A1","BMI1","OLIG2",
                 "MYC","NOTCH1","FZD7","L1CAM","CD133","EZH2","NANOG"]
DIFF_SIGNATURE = ["GFAP","S100B","AQP4",           # astrocyte
                  "MAP2","TUBB3","DCX","NEFL",       # neuron
                  "MBP","MOG","CNP",                 # oligodendrocyte
                  "TMEM119","CX3CR1","P2RY12"]       # microglia

def score_sig(ad, genes, name):
    avail = [g for g in genes if g in ad.var_names]
    print(f"  {name}: {len(avail)}/{len(genes)} genes found")
    if len(avail) == 0: return np.zeros(ad.n_obs)
    expr = pd.DataFrame(
        ad.raw.X[:, [ad.raw.var_names.get_loc(g) for g in avail]].toarray()
        if hasattr(ad.raw.X, "toarray") else
        ad.raw.X[:, [list(ad.raw.var_names).index(g) for g in avail]],
        index=ad.obs.index, columns=avail
    )
    return expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-9)).mean(axis=1).values

adata_tumor.obs["score_gsc"]  = score_sig(adata_tumor, GSC_SIGNATURE,  "GSC signature")
adata_tumor.obs["score_diff"] = score_sig(adata_tumor, DIFF_SIGNATURE, "Differentiation signature")
adata_tumor.obs["stemness"]   = adata_tumor.obs["score_gsc"] - adata_tumor.obs["score_diff"]

q75 = adata_tumor.obs["stemness"].quantile(0.75)
q25 = adata_tumor.obs["stemness"].quantile(0.25)
adata_tumor.obs["gsc_label"] = "gsc_mid"
adata_tumor.obs.loc[adata_tumor.obs["stemness"] >= q75, "gsc_label"] = "gsc_high"
adata_tumor.obs.loc[adata_tumor.obs["stemness"] <= q25, "gsc_label"] = "gsc_low"

n_hi = (adata_tumor.obs["gsc_label"] == "gsc_high").sum()
n_lo = (adata_tumor.obs["gsc_label"] == "gsc_low").sum()
print(f"\n  GSC-high: {n_hi}  |  GSC-low: {n_lo}  |  mid: {adata_tumor.n_obs - n_hi - n_lo}")

# Plot stemness distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(adata_tumor.obs["stemness"], bins=40, color="steelblue", edgecolor="white")
axes[0].axvline(q75, color="red",  linestyle="--", label="Q75 (GSC-high)")
axes[0].axvline(q25, color="blue", linestyle="--", label="Q25 (GSC-low)")
axes[0].set_xlabel("Stemness score"); axes[0].set_ylabel("Cells")
axes[0].set_title("GBM Stemness Score Distribution"); axes[0].legend()
axes[1].scatter(adata_tumor.obsm["X_umap"][:,0], adata_tumor.obsm["X_umap"][:,1],
                c=adata_tumor.obs["stemness"], cmap="RdBu_r", s=5, alpha=0.7)
axes[1].set_title("UMAP — GBM stemness"); axes[1].set_xlabel("UMAP1"); axes[1].set_ylabel("UMAP2")
plt.colorbar(axes[1].collections[0], ax=axes[1], label="stemness")
plt.tight_layout()
plt.savefig("results/figures/B1_gbm_stemness.png", dpi=120, bbox_inches="tight")
plt.close()
print("  Saved: results/figures/B1_gbm_stemness.png")

# ─────────────────────────────────────────────────────────────────
# B1e — DIFFERENTIAL EXPRESSION (Stage 1 for GBM)
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B1e — DIFFERENTIAL EXPRESSION (GSC-high vs GSC-low)")
print("=" * 60)

from scipy.stats import mannwhitneyu

adata_de = adata_tumor[adata_tumor.obs["gsc_label"].isin(["gsc_high","gsc_low"])].copy()
hi_idx   = adata_de.obs["gsc_label"] == "gsc_high"
lo_idx   = adata_de.obs["gsc_label"] == "gsc_low"

X = (adata_de.raw.X.toarray() if hasattr(adata_de.raw.X, "toarray")
     else np.array(adata_de.raw.X))
genes = list(adata_de.raw.var_names)

print(f"  Running Wilcoxon: {hi_idx.sum()} GSC-high vs {lo_idx.sum()} GSC-low")
rows = []
for j, gene in enumerate(genes):
    hi_vals = X[hi_idx.values, j]
    lo_vals = X[lo_idx.values, j]
    if hi_vals.mean() == 0 and lo_vals.mean() == 0: continue
    try:
        u, p = mannwhitneyu(hi_vals, lo_vals, alternative="greater")
    except Exception:
        continue
    n  = len(hi_vals) + len(lo_vals)
    w  = u / (len(hi_vals) * len(lo_vals))
    lfc = np.log2(hi_vals.mean() + 1) - np.log2(lo_vals.mean() + 1)
    rows.append({"gene_symbol": gene, "wilcoxon_score": w, "log2fc": lfc, "pval": p})
    if j % 2000 == 0: print(f"    {j}/{len(genes)} genes...")

de_df = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
_, padj, _, _ = multipletests(de_df["pval"].values, method="fdr_bh")
de_df["pval_adj"] = padj
de_df = de_df.sort_values("wilcoxon_score", ascending=False)

markers = de_df[(de_df["pval_adj"] < 0.05) & (de_df["log2fc"] > 0.5)]
print(f"\n  GBM-CSC markers (padj<0.05, log2FC>0.5): {len(markers)}")
print(f"  Top 15:")
print(markers.head(15)[["gene_symbol","wilcoxon_score","log2fc","pval_adj"]].to_string(index=False))

markers.to_csv("results/tables/B1_gbm_csc_markers.csv", index=False)
print("  Saved: results/tables/B1_gbm_csc_markers.csv")

# ─────────────────────────────────────────────────────────────────
# B1f — CROSS-CANCER COMPARISON
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B1f — CROSS-CANCER COMPARISON (Breast vs GBM)")
print("=" * 60)

brca_s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
brca_s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
CONSENSUS = ['B2M','CD44','EEF1A1','EPCAM','FN1','HLA-E','HMGA1','ID3',
             'KLF4','MYC','RPL22L1','S100A6','SERBP1','SERPINE2','SOX9','VIM']

brca_de_genes  = set(brca_s1["gene_symbol"].head(200))
brca_attn_genes = set(brca_s2["gene_symbol"].head(200))
gbm_genes      = set(markers["gene_symbol"].head(200))

shared_de_gbm    = brca_de_genes  & gbm_genes
shared_attn_gbm  = brca_attn_genes & gbm_genes
shared_all_three = brca_de_genes & brca_attn_genes & gbm_genes
jaccard_de       = len(shared_de_gbm) / len(brca_de_genes | gbm_genes)
jaccard_attn     = len(shared_attn_gbm) / len(brca_attn_genes | gbm_genes)

print(f"\n  Breast DE   ∩ GBM: {len(shared_de_gbm)} genes  (Jaccard={jaccard_de:.3f})")
print(f"  Breast Attn ∩ GBM: {len(shared_attn_gbm)} genes  (Jaccard={jaccard_attn:.3f})")
print(f"  All three methods : {len(shared_all_three)} genes  ← pan-cancer CSC program")

print(f"\n  Shared by ALL THREE (pan-cancer CSC genes):")
pan_cancer = sorted(shared_all_three)
for g in pan_cancer:
    gbm_lfc   = markers.loc[markers["gene_symbol"]==g,"log2fc"].values
    brca_lfc  = brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values
    brca_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
    g_lfc   = f"{gbm_lfc[0]:+.2f}"   if len(gbm_lfc)   else "N/A"
    b_lfc   = f"{brca_lfc[0]:+.2f}"  if len(brca_lfc)  else "N/A"
    b_attn  = f"{brca_attn[0]:.3f}"  if len(brca_attn) else "N/A"
    print(f"    {g:12s}  BRCA_DE log2FC={b_lfc}  BRCA_attn={b_attn}  GBM_log2FC={g_lfc}")

# Also check: does the breast CSC Geneformer attention find known GBM-CSC genes?
known_gbm_csc = {"SOX2","NES","PROM1","CD44","BMI1","OLIG2","L1CAM","EZH2","NANOG","FZD7"}
brca_attn_in_gbm_known = brca_attn_genes & known_gbm_csc
print(f"\n  Known GBM-CSC genes found by Breast Geneformer attention:")
for g in sorted(brca_attn_in_gbm_known):
    rank = brca_s2[brca_s2["gene_symbol"]==g].index[0] + 1 if g in brca_s2["gene_symbol"].values else "?"
    print(f"    {g}  (Geneformer rank #{rank})")

# Save pan-cancer table
pan_df = pd.DataFrame([{
    "gene": g,
    "brca_de_log2fc":    brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values[0] if g in brca_de_genes else None,
    "brca_attn_score":   brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values[0] if g in brca_attn_genes else None,
    "gbm_de_log2fc":     markers.loc[markers["gene_symbol"]==g,"log2fc"].values[0] if g in gbm_genes else None,
} for g in sorted(shared_all_three | (brca_attn_genes & gbm_genes))])
pan_df.to_csv("results/tables/B1_pancancer_consensus.csv", index=False)

# ── Comparison figure ──────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Venn-style bar chart
cats   = ["BRCA DE only","GBM DE only","BRCA Attn ∩ GBM","All three"]
counts = [len(brca_de_genes - gbm_genes - brca_attn_genes),
          len(gbm_genes - brca_de_genes - brca_attn_genes),
          len(shared_attn_gbm - brca_de_genes),
          len(shared_all_three)]
colors = ["#d62728","#ff7f0e","#1f77b4","#2ca02c"]
axes[0].bar(cats, counts, color=colors)
for i,c in enumerate(counts):
    axes[0].text(i, c+0.5, str(c), ha="center", fontweight="bold")
axes[0].set_title("Cross-cancer gene overlap\n(top 200 each)")
axes[0].set_ylabel("Genes"); axes[0].tick_params(axis="x", rotation=20)

# GBM top DE genes
top_gbm = markers.head(20)
axes[1].barh(top_gbm["gene_symbol"][::-1], top_gbm["log2fc"][::-1], color="#ff7f0e")
pan_in_top = set(top_gbm["gene_symbol"]) & shared_all_three
axes[1].set_title(f"Top 20 GBM-CSC markers\n(pan-cancer genes: {pan_in_top})")
axes[1].set_xlabel("log2 Fold Change")

# Scatter: BRCA log2FC vs GBM log2FC for shared genes
shared_both = brca_de_genes & gbm_genes
if shared_both:
    x = [brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values[0] for g in shared_both
         if g in brca_s1["gene_symbol"].values]
    y = [markers.loc[markers["gene_symbol"]==g,"log2fc"].values[0] for g in shared_both
         if g in brca_s1["gene_symbol"].values]
    labs = [g for g in shared_both if g in brca_s1["gene_symbol"].values]
    axes[2].scatter(x, y, s=30, alpha=0.7, color="steelblue")
    for xi, yi, lab in zip(x, y, labs):
        if lab in shared_all_three:
            axes[2].annotate(lab, (xi, yi), fontsize=7, color="red")
    axes[2].axhline(0, color="gray", linewidth=0.5)
    axes[2].axvline(0, color="gray", linewidth=0.5)
    r, p = __import__("scipy.stats",fromlist=["pearsonr"]).pearsonr(x, y)
    axes[2].set_title(f"BRCA vs GBM log2FC\n(shared genes, r={r:.2f}, p={p:.3f})")
    axes[2].set_xlabel("BRCA log2FC"); axes[2].set_ylabel("GBM log2FC")

plt.suptitle("Pan-Cancer CSC Validation: Breast Cancer vs GBM", fontsize=13)
plt.tight_layout()
plt.savefig("results/figures/B1_pancancer_comparison.png", dpi=120, bbox_inches="tight")
plt.close()

print("\nSaved: results/figures/B1_pancancer_comparison.png")
print("Saved: results/tables/B1_pancancer_consensus.csv")
print("\n✓ Phase B1 complete.")
print(f"  Pan-cancer CSC program: {len(shared_all_three)} genes confirmed across breast + GBM")
print(f"  Geneformer breast attention ∩ GBM DE: {len(shared_attn_gbm)} genes")
