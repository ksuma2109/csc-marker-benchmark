# Phase B2 — GBM Refined: Official Cell Type Annotations from GEO + Tumor-only DE
#
# Fix for B1: GBM tumor cells mimic astrocyte/OPC transcriptomes, so
#   marker-based annotation failed. Instead, use official Darmanis 2017
#   annotations from the GEO series matrix (line "cell type: Neoplastic").
#
# Cell IDs are reconstructed as "{plate_id}.{well}" — matching our adata.obs_names.
#
# Pipeline:
#   B2a: Parse GEO series matrix → neoplastic labels → map to cell IDs
#   B2b: Filter adata to neoplastic (tumor) cells only
#   B2c: GSC scoring within tumor cells (Q1 vs Q4)
#   B2d: Wilcoxon DE (GSC-high vs GSC-low, tumor only)
#   B2e: Updated pan-cancer comparison
#   B2f: Summary figures + heatmap
#
# Output: results/tables/B2_gbm_tumor_markers.csv
#         results/tables/B2_pancancer_final.csv
#         results/figures/B2_*.png

import os, gzip, gc, requests
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

os.makedirs("results/tables",  exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# ─────────────────────────────────────────────────────────────────
# B2a — PARSE GEO SERIES MATRIX FOR OFFICIAL CELL TYPE LABELS
# ─────────────────────────────────────────────────────────────────
print("=" * 60)
print("B2a — PARSING GEO SERIES MATRIX FOR CELL LABELS")
print("=" * 60)

MATRIX_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/"
              "matrix/GSE84465_series_matrix.txt.gz")
MATRIX_CACHE = "data/raw/GSE84465/series_matrix.txt.gz"

if not os.path.exists(MATRIX_CACHE):
    print("  Downloading series matrix...")
    r = requests.get(MATRIX_URL, stream=True, timeout=120)
    r.raise_for_status()
    with open(MATRIX_CACHE, "wb") as f:
        for chunk in r.iter_content(8192):
            f.write(chunk)
    print("  Done.")

with gzip.open(MATRIX_CACHE, "rt", errors="ignore") as f:
    content = f.read()

lines = content.split("\n")
plate_line = well_line = celltype_line = neoplastic_line = None
for line in lines:
    if "plate id:" in line:
        plate_line = line
    elif "well:" in line:
        well_line = line
    elif "cell type:" in line:
        celltype_line = line
    elif "neoplastic:" in line:
        neoplastic_line = line

def parse_fields(line):
    parts = line.split("\t")
    return [p.strip().strip('"').split(": ", 1)[-1] for p in parts[1:]]

plates     = parse_fields(plate_line)
wells      = parse_fields(well_line)
celltypes  = parse_fields(celltype_line)
neoplastic = parse_fields(neoplastic_line)

cell_meta = pd.DataFrame({
    "plate":       plates,
    "well":        wells,
    "cell_type":   celltypes,
    "neoplastic":  neoplastic,
})
cell_meta["cell_id"] = cell_meta["plate"] + "." + cell_meta["well"]

print(f"  Total annotated cells: {len(cell_meta)}")
print(f"  Cell type distribution:\n{cell_meta['cell_type'].value_counts().to_string()}")
print(f"\n  Neoplastic breakdown:\n{cell_meta['neoplastic'].value_counts().to_string()}")

neoplastic_ids = set(cell_meta.loc[cell_meta["neoplastic"] == "Neoplastic", "cell_id"])
print(f"\n  Neoplastic cell IDs: {len(neoplastic_ids)}")

# ─────────────────────────────────────────────────────────────────
# LOAD PREPROCESSED GBM DATA
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("LOADING GBM PREPROCESSED DATA")
print("=" * 60)

adata = sc.read_h5ad("data/processed/gbm_B1_preprocessed.h5ad")
print(f"  Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Map cell type labels to adata
adata.obs["cell_id"] = adata.obs_names
cell_meta_idx = cell_meta.set_index("cell_id")
adata.obs["official_celltype"] = adata.obs["cell_id"].map(
    cell_meta_idx["cell_type"]).fillna("Unknown")
adata.obs["is_neoplastic"] = adata.obs["cell_id"].isin(neoplastic_ids)

mapped = adata.obs["is_neoplastic"].sum()
print(f"  Mapped neoplastic cells: {mapped}")
print(f"  Cell type distribution in adata:\n{adata.obs['official_celltype'].value_counts().to_string()}")

# ─────────────────────────────────────────────────────────────────
# B2b — ISOLATE NEOPLASTIC (TUMOR) CELLS
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B2b — ISOLATING NEOPLASTIC CELLS")
print("=" * 60)

adata_tumor = adata[adata.obs["is_neoplastic"]].copy()
print(f"  Neoplastic cells: {adata_tumor.n_obs}")
print(f"  Non-neoplastic excluded: {adata.n_obs - adata_tumor.n_obs}")

# ─────────────────────────────────────────────────────────────────
# B2c — GSC SCORING WITHIN NEOPLASTIC CELLS
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B2c — GSC SCORING (NEOPLASTIC CELLS)")
print("=" * 60)

var_names = list(adata_tumor.raw.var_names)

# GSC signature: neural stem cell / GBM stem cell markers
# Excludes GFAP/AQP4 (normal astrocyte) to avoid confounding
GSC_SIG  = ["SOX2","PROM1","CD44","MYC","NOTCH1","FZD7","L1CAM","EZH2",
             "BMI1","ALDH1A1","ALDH1A3","ITGA6","NES","ABCG2","SOX9"]
DIFF_SIG = ["TUBB3","DCX","NEFL","MAP2",    # neuronal differentiation
            "S100B","MBP","MOG","PLP1"]      # glial differentiation

def z_score_sig(ad, genes, label):
    avail = [g for g in genes if g in var_names]
    print(f"  {label}: {len(avail)}/{len(genes)} genes")
    if not avail: return np.zeros(ad.n_obs)
    idx = [var_names.index(g) for g in avail]
    X   = (ad.raw.X.toarray() if hasattr(ad.raw.X, "toarray")
           else np.array(ad.raw.X))
    expr = X[:, idx].astype(float)
    z    = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-9)
    return z.mean(axis=1)

adata_tumor.obs["score_gsc"]  = z_score_sig(adata_tumor, GSC_SIG,  "GSC signature")
adata_tumor.obs["score_diff"] = z_score_sig(adata_tumor, DIFF_SIG, "Differentiation")
adata_tumor.obs["stemness"]   = adata_tumor.obs["score_gsc"] - adata_tumor.obs["score_diff"]

q75 = adata_tumor.obs["stemness"].quantile(0.75)
q25 = adata_tumor.obs["stemness"].quantile(0.25)
adata_tumor.obs["gsc_label"] = "mid"
adata_tumor.obs.loc[adata_tumor.obs["stemness"] >= q75, "gsc_label"] = "gsc_high"
adata_tumor.obs.loc[adata_tumor.obs["stemness"] <= q25, "gsc_label"] = "gsc_low"

n_hi = (adata_tumor.obs["gsc_label"] == "gsc_high").sum()
n_lo = (adata_tumor.obs["gsc_label"] == "gsc_low").sum()
print(f"\n  Neoplastic GSC-high: {n_hi}  |  GSC-low: {n_lo}")

# ─────────────────────────────────────────────────────────────────
# B2d — DE: GSC-high vs GSC-low (neoplastic cells only)
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B2d — DIFFERENTIAL EXPRESSION (NEOPLASTIC GSC-HIGH vs GSC-LOW)")
print("=" * 60)

adata_de = adata_tumor[adata_tumor.obs["gsc_label"].isin(["gsc_high","gsc_low"])].copy()
hi_idx   = (adata_de.obs["gsc_label"] == "gsc_high").values
lo_idx   = (adata_de.obs["gsc_label"] == "gsc_low").values
X_de = (adata_de.raw.X.toarray() if hasattr(adata_de.raw.X, "toarray")
        else np.array(adata_de.raw.X))

print(f"  {hi_idx.sum()} GSC-high vs {lo_idx.sum()} GSC-low (neoplastic only)")

rows = []
for j, gene in enumerate(var_names):
    hi_vals = X_de[hi_idx, j]
    lo_vals = X_de[lo_idx, j]
    if hi_vals.mean() == 0 and lo_vals.mean() == 0:
        continue
    try:
        u, p = mannwhitneyu(hi_vals, lo_vals, alternative="greater")
    except Exception:
        continue
    w   = u / (hi_idx.sum() * lo_idx.sum())
    lfc = np.log2(hi_vals.mean() + 1) - np.log2(lo_vals.mean() + 1)
    rows.append({"gene_symbol": gene, "wilcoxon_score": w, "log2fc": lfc, "pval": p})
    if j % 3000 == 0:
        print(f"    {j}/{len(var_names)} genes...")

de_df = pd.DataFrame(rows)
_, padj, _, _ = multipletests(de_df["pval"].values, method="fdr_bh")
de_df["pval_adj"] = padj
de_df = de_df.sort_values("wilcoxon_score", ascending=False)

markers_b2 = de_df[(de_df["pval_adj"] < 0.05) & (de_df["log2fc"] > 0.25)]
print(f"\n  GBM-CSC markers (padj<0.05, log2FC>0.25): {len(markers_b2)}")
print(f"  Top 20:")
print(markers_b2.head(20)[["gene_symbol","wilcoxon_score","log2fc","pval_adj"]].to_string(index=False))

markers_b2.to_csv("results/tables/B2_gbm_tumor_markers.csv", index=False)
print("  Saved: results/tables/B2_gbm_tumor_markers.csv")

# ─────────────────────────────────────────────────────────────────
# B2e — PAN-CANCER COMPARISON
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B2e — PAN-CANCER COMPARISON (BRCA vs GBM NEOPLASTIC)")
print("=" * 60)

brca_s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
brca_s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
CONSENSUS = ['B2M','CD44','EEF1A1','EPCAM','FN1','HLA-E','HMGA1','ID3',
             'KLF4','MYC','RPL22L1','S100A6','SERBP1','SERPINE2','SOX9','VIM']

brca_de_top    = set(brca_s1["gene_symbol"].head(200))
brca_attn_top  = set(brca_s2["gene_symbol"].head(200))
gbm_top        = set(markers_b2["gene_symbol"].head(200))

shared_de_gbm   = brca_de_top  & gbm_top
shared_attn_gbm = brca_attn_top & gbm_top
shared_all3     = brca_de_top  & brca_attn_top & gbm_top
consensus_gbm   = set(CONSENSUS) & gbm_top

j_de   = len(shared_de_gbm)   / max(1, len(brca_de_top   | gbm_top))
j_attn = len(shared_attn_gbm) / max(1, len(brca_attn_top | gbm_top))

print(f"  BRCA DE   ∩ GBM neoplastic: {len(shared_de_gbm)} genes  (Jaccard={j_de:.3f})")
print(f"  BRCA Attn ∩ GBM neoplastic: {len(shared_attn_gbm)} genes  (Jaccard={j_attn:.3f})")
print(f"  All three methods (pan-CSC): {len(shared_all3)} genes")
print(f"  Consensus genes in GBM:      {len(consensus_gbm)} / {len(CONSENSUS)}")

print(f"\n  ── PAN-CANCER CSC GENES (all 3 methods) ──")
for g in sorted(shared_all3):
    b_lfc  = brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values
    b_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
    g_lfc  = markers_b2.loc[markers_b2["gene_symbol"]==g,"log2fc"].values
    tag    = "★" if g in CONSENSUS else ""
    print(f"    {tag}{g:12s}  BRCA_DE={b_lfc[0]:+.2f}  BRCA_attn={b_attn[0]:.3f}  GBM={g_lfc[0]:+.2f}")

print(f"\n  ── BRCA GENEFORMER ∩ GBM NEOPLASTIC ──")
for g in sorted(shared_attn_gbm):
    b_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
    g_lfc  = markers_b2.loc[markers_b2["gene_symbol"]==g,"log2fc"].values
    rank   = list(brca_s2["gene_symbol"]).index(g) + 1 if g in brca_s2["gene_symbol"].values else "?"
    print(f"    {g:12s}  BRCA_attn={b_attn[0]:.3f}(rank#{rank})  GBM={g_lfc[0]:+.2f}")

print(f"\n  ── CONSENSUS GENES IN GBM ──")
for g in sorted(consensus_gbm):
    g_lfc = markers_b2.loc[markers_b2["gene_symbol"]==g,"log2fc"].values
    print(f"    {g:12s}  GBM_log2FC={g_lfc[0]:+.2f}")

# Save final table
pan_rows = []
for g in sorted(shared_de_gbm | shared_attn_gbm):
    pan_rows.append({
        "gene":           g,
        "in_brca_de":    g in brca_de_top,
        "in_brca_attn":  g in brca_attn_top,
        "in_gbm_neo":    g in gbm_top,
        "in_consensus":  g in CONSENSUS,
        "brca_de_log2fc":  brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values[0] if g in brca_de_top else None,
        "brca_attn_score": brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values[0] if g in brca_attn_top else None,
        "gbm_de_log2fc":   markers_b2.loc[markers_b2["gene_symbol"]==g,"log2fc"].values[0] if g in gbm_top else None,
    })
pan_df = pd.DataFrame(pan_rows)
pan_df.to_csv("results/tables/B2_pancancer_final.csv", index=False)
print(f"\n  Saved: results/tables/B2_pancancer_final.csv  ({len(pan_df)} genes)")

# ─────────────────────────────────────────────────────────────────
# B2f — FIGURES
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B2f — GENERATING FIGURES")
print("=" * 60)

CT_PALETTE = {
    "Neoplastic":"#d62728","Oligodendrocyte":"#8c564b","Microglia":"#1f77b4",
    "Astrocyte":"#2ca02c","Neuron":"#9467bd","Endothelial":"#7f7f7f","Unknown":"#cccccc",
}

# Figure 1: UMAP with official annotations + stemness
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

colors_ct = [CT_PALETTE.get(ct, "#cccccc") for ct in adata.obs["official_celltype"]]
axes[0].scatter(adata.obsm["X_umap"][:,0], adata.obsm["X_umap"][:,1],
                c=colors_ct, s=3, alpha=0.6)
patches = [mpatches.Patch(color=v, label=k) for k, v in CT_PALETTE.items()
           if k in adata.obs["official_celltype"].values]
axes[0].legend(handles=patches, fontsize=7, markerscale=2)
axes[0].set_title(f"GBM Cell Types (Darmanis 2017)\n"
                  f"Neoplastic: {adata.obs['is_neoplastic'].sum()}")
axes[0].set_xlabel("UMAP1"); axes[0].set_ylabel("UMAP2")

# Tumor only — stemness
axes[1].scatter(adata_tumor.obsm["X_umap"][:,0], adata_tumor.obsm["X_umap"][:,1],
                c=adata_tumor.obs["stemness"], cmap="RdYlGn", s=8, alpha=0.8,
                vmin=-2, vmax=2)
axes[1].set_title(f"Neoplastic cells: stemness score\n({adata_tumor.n_obs} cells)")
axes[1].set_xlabel("UMAP1"); axes[1].set_ylabel("UMAP2")
plt.colorbar(axes[1].collections[0], ax=axes[1], label="stemness (GSC−diff)")

# Top DE genes
top20 = markers_b2.head(20) if len(markers_b2) >= 1 else de_df.head(20)
bar_colors = []
for g in top20["gene_symbol"]:
    if g in shared_all3:        bar_colors.append("#d62728")
    elif g in shared_attn_gbm:  bar_colors.append("#ff7f0e")
    elif g in CONSENSUS:         bar_colors.append("#2ca02c")
    else:                        bar_colors.append("#aaaaaa")
axes[2].barh(top20["gene_symbol"][::-1], top20["log2fc"][::-1],
             color=bar_colors[::-1])
axes[2].set_xlabel("log2 Fold Change (GSC-high vs low)")
axes[2].set_title("Top GBM-CSC Markers (neoplastic cells)")
leg_patches = [mpatches.Patch(color="#d62728",label="Pan-cancer (all 3)"),
               mpatches.Patch(color="#ff7f0e",label="Geneformer∩GBM"),
               mpatches.Patch(color="#2ca02c",label="Consensus"),
               mpatches.Patch(color="#aaaaaa",label="GBM only")]
axes[2].legend(handles=leg_patches, fontsize=7, loc="lower right")

plt.suptitle("Phase B2: GBM Neoplastic Cell CSC Analysis (Darmanis et al. 2017)", fontsize=12)
plt.tight_layout()
plt.savefig("results/figures/B2_gbm_neoplastic.png", dpi=120, bbox_inches="tight")
plt.close()
print("  Saved: results/figures/B2_gbm_neoplastic.png")

# Figure 2: cross-cancer evidence heatmap
evidence_genes = sorted(shared_de_gbm | shared_attn_gbm | (set(CONSENSUS) & gbm_top))
if evidence_genes:
    fig, ax = plt.subplots(figsize=(max(8, len(evidence_genes)*0.45), 5))
    dm = []
    for g in evidence_genes:
        b_lfc  = brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values
        b_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
        g_lfc  = markers_b2.loc[markers_b2["gene_symbol"]==g,"log2fc"].values
        dm.append([
            float(b_lfc[0])  if len(b_lfc)  else 0,
            float(b_attn[0]) if len(b_attn) else 0,
            float(g_lfc[0])  if len(g_lfc)  else 0,
        ])
    dm = np.array(dm)
    dm_norm = (dm - dm.min(axis=0)) / (dm.max(axis=0) - dm.min(axis=0) + 1e-9)
    im = ax.imshow(dm_norm.T, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(evidence_genes)))
    ax.set_xticklabels(evidence_genes, rotation=45, ha="right", fontsize=9)
    ax.set_yticks([0,1,2])
    ax.set_yticklabels(["BRCA DE log2FC","BRCA Geneformer attn","GBM neoplastic log2FC"],
                       fontsize=10)
    for i, g in enumerate(evidence_genes):
        if g in shared_all3:
            ax.add_patch(plt.Rectangle((i-0.5,-0.5),1,3,
                         fill=False, edgecolor="blue", lw=2))
    plt.colorbar(im, ax=ax, label="normalized score")
    ax.set_title("Pan-Cancer CSC Gene Evidence Matrix\n"
                 "(blue border = confirmed by all 3 methods)")
    plt.tight_layout()
    plt.savefig("results/figures/B2_pancancer_heatmap.png", dpi=120, bbox_inches="tight")
    plt.close()
    print("  Saved: results/figures/B2_pancancer_heatmap.png")

print("\n" + "=" * 60)
print("PHASE B2 COMPLETE")
print("=" * 60)
print(f"  GBM neoplastic cells analysed: {adata_tumor.n_obs}")
print(f"  GBM-CSC DE markers (padj<0.05, FC>0.25): {len(markers_b2)}")
print(f"  Pan-cancer program (all 3 methods): {len(shared_all3)} genes — {sorted(shared_all3)}")
print(f"  Geneformer∩GBM: {len(shared_attn_gbm)} genes — {sorted(shared_attn_gbm)}")
print(f"  Consensus validated in GBM: {len(consensus_gbm)} genes — {sorted(consensus_gbm)}")
