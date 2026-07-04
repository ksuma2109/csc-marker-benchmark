# Phase B3 — Cross-Cancer Validation #3: Melanoma (Tirosh et al. 2016, GSE72056)
#
# Dataset: GSE72056 — 4,645 single cells from 19 melanoma patients
#          Tirosh et al. 2016 Science
#          Malignant vs non-malignant annotated by CNV inference (column 2 = 0/1/2)
#          Column layout: cell, malignant(1=no, 2=yes, 0=unresolved), ...
#
# Goal: Identify melanoma cancer stem cell (MelCSC) markers → compare with
#       breast CSC (Stage 1 + Geneformer) and GBM neoplastic markers →
#       define final pan-cancer CSC gene program
#
# Output: results/tables/B3_melanoma_csc_markers.csv
#         results/tables/B3_pancancer_3cancer.csv
#         results/figures/B3_*.png

import os, gzip, io, requests
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

os.makedirs("data/raw/GSE72056",  exist_ok=True)
os.makedirs("data/processed",     exist_ok=True)
os.makedirs("results/tables",     exist_ok=True)
os.makedirs("results/figures",    exist_ok=True)

# ─────────────────────────────────────────────────────────────────
# B3a — DOWNLOAD MELANOMA DATA
# ─────────────────────────────────────────────────────────────────
print("=" * 60)
print("B3a — DOWNLOADING GSE72056 (MELANOMA)")
print("=" * 60)

URL  = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/"
        "suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz")
DEST = "data/raw/GSE72056/melanoma.txt.gz"

if not os.path.exists(DEST):
    print("  Downloading...")
    r = requests.get(URL, stream=True, timeout=180)
    r.raise_for_status()
    with open(DEST, "wb") as f:
        for chunk in r.iter_content(8192):
            f.write(chunk)
    print("  Done.")
else:
    print("  Already downloaded.")

# ─────────────────────────────────────────────────────────────────
# B3b — LOAD + PARSE
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B3b — LOADING AND PARSING MELANOMA DATA")
print("=" * 60)

H5AD = "data/processed/melanoma_B3.h5ad"

if os.path.exists(H5AD):
    print(f"  Loading cached: {H5AD}")
    adata = sc.read_h5ad(H5AD)
else:
    with gzip.open(DEST, "rt") as f:
        raw = pd.read_csv(f, sep="\t", index_col=0)
    # rows = genes (index by name), cols = cells
    # Special metadata rows at top (accessed by index label):
    #   "tumor"                             → patient/tumor ID
    #   "malignant(1=no,2=yes,0=unresolved)"→ 1=non-mal, 2=malignant, 0=unresolved
    #   "non-malignant cell type..."        → immune/stromal subtype
    # Data is already log2(TPM+1) normalized by the paper
    MAL_ROW  = "malignant(1=no,2=yes,0=unresolved)"
    TYPE_ROW = "non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)"
    TUM_ROW  = "tumor"
    meta_idx = [TUM_ROW, MAL_ROW, TYPE_ROW]
    meta_rows = raw.loc[meta_idx].T.copy()
    meta_rows.columns = ["tumor_id", "malignant", "cell_type"]
    counts = raw.drop(index=meta_idx).T   # cells × genes
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)

    print(f"  Shape (cells × genes): {counts.shape}")
    print(f"  Malignant values: {meta_rows['malignant'].value_counts().to_dict()}")
    # malignant: 1=non-malignant, 2=malignant, 0=unresolved

    adata = sc.AnnData(X=counts.values.astype("float32"),
                       obs=pd.DataFrame(index=counts.index),
                       var=pd.DataFrame(index=counts.columns))
    adata.obs["malignant"] = meta_rows["malignant"].astype(str).values
    adata.obs["tumor_id"]  = meta_rows["tumor_id"].astype(str).values
    adata.obs["cell_type"] = meta_rows["cell_type"].astype(str).values

    # Data is already log2(TPM+1) — skip normalize/log1p
    # QC by gene count only
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    before = adata.n_obs
    adata = adata[(adata.obs["n_genes_by_counts"] >= 200)].copy()
    print(f"  After QC: {adata.n_obs} cells (removed {before - adata.n_obs})")
    adata.raw = adata.copy()  # store raw (already log-normalized)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", subset=False)
    sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, random_state=42)

    adata.write_h5ad(H5AD)
    print(f"  Saved: {H5AD}")

print(f"  Dataset: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"  Malignant: {(adata.obs['malignant']=='2').sum()}  Non-malignant: {(adata.obs['malignant']=='1').sum()}")

# ─────────────────────────────────────────────────────────────────
# B3c — ISOLATE MALIGNANT CELLS + SCORE FOR STEMNESS
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B3c — MELANOMA STEM CELL (MelSC) SCORING")
print("=" * 60)

# malignant stored as "2.0" (float→str); normalize to int string
adata.obs["malignant"] = pd.to_numeric(adata.obs["malignant"], errors="coerce").fillna(0).astype(int).astype(str)
adata_mal = adata[adata.obs["malignant"] == "2"].copy()
print(f"  Malignant cells: {adata_mal.n_obs}")
print(f"  Unique malignant values: {adata.obs['malignant'].unique()}")

var_names = list(adata_mal.raw.var_names)

# Melanoma stem cell signature
# Sources: Hoek et al. 2006, 2008; Boiko et al. 2010; Shakhova et al. 2012
MELSC_SIG  = ["ALDH1A1","ALDH3A1","CD44","SOX9","SOX10","NGFR","ABCB5",
              "CD271","JARID1B","KLF4","MYC","NES","NOTCH1","FZD7","VIM",
              "SNAI2","TWIST1","ZEB1","CD133","PROM1"]
DIFF_SIG   = ["MLANA","TYR","DCT","TYRP1","SILV","MC1R",    # melanocyte diff
              "MKI67","TOP2A","PCNA"]                        # proliferating (exclude)

def score_sig(ad, genes, label):
    avail = [g for g in genes if g in var_names]
    print(f"  {label}: {len(avail)}/{len(genes)} genes: {avail}")
    if not avail: return np.zeros(ad.n_obs)
    idx  = [var_names.index(g) for g in avail]
    X    = (ad.raw.X.toarray() if hasattr(ad.raw.X, "toarray") else np.array(ad.raw.X))
    expr = X[:, idx].astype(float)
    z    = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-9)
    return z.mean(axis=1)

adata_mal.obs["score_melsc"] = score_sig(adata_mal, MELSC_SIG, "MelSC signature")
adata_mal.obs["score_diff"]  = score_sig(adata_mal, DIFF_SIG,  "Differentiated")
adata_mal.obs["stemness"]    = adata_mal.obs["score_melsc"] - adata_mal.obs["score_diff"]

q75 = adata_mal.obs["stemness"].quantile(0.75)
q25 = adata_mal.obs["stemness"].quantile(0.25)
adata_mal.obs["sc_label"] = "mid"
adata_mal.obs.loc[adata_mal.obs["stemness"] >= q75, "sc_label"] = "sc_high"
adata_mal.obs.loc[adata_mal.obs["stemness"] <= q25, "sc_label"] = "sc_low"

n_hi = (adata_mal.obs["sc_label"] == "sc_high").sum()
n_lo = (adata_mal.obs["sc_label"] == "sc_low").sum()
print(f"\n  MelSC-high: {n_hi}  |  MelSC-low: {n_lo}")

# ─────────────────────────────────────────────────────────────────
# B3d — WILCOXON DE
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B3d — DE: MelSC-HIGH vs MelSC-LOW")
print("=" * 60)

adata_de = adata_mal[adata_mal.obs["sc_label"].isin(["sc_high","sc_low"])].copy()
hi_idx   = (adata_de.obs["sc_label"] == "sc_high").values
lo_idx   = (adata_de.obs["sc_label"] == "sc_low").values
X_de = (adata_de.raw.X.toarray() if hasattr(adata_de.raw.X, "toarray")
        else np.array(adata_de.raw.X))
print(f"  {hi_idx.sum()} MelSC-high vs {lo_idx.sum()} MelSC-low")

rows = []
for j, gene in enumerate(var_names):
    hi_v, lo_v = X_de[hi_idx, j], X_de[lo_idx, j]
    if hi_v.mean() == 0 and lo_v.mean() == 0: continue
    try:
        u, p = mannwhitneyu(hi_v, lo_v, alternative="greater")
    except Exception:
        continue
    w   = u / (hi_idx.sum() * lo_idx.sum())
    lfc = np.log2(hi_v.mean() + 1) - np.log2(lo_v.mean() + 1)
    rows.append({"gene_symbol": gene, "wilcoxon_score": w, "log2fc": lfc, "pval": p})
    if j % 3000 == 0: print(f"    {j}/{len(var_names)} genes...")

de_df = pd.DataFrame(rows)
_, padj, _, _ = multipletests(de_df["pval"].values, method="fdr_bh")
de_df["pval_adj"] = padj
de_df = de_df.sort_values("wilcoxon_score", ascending=False)

markers_mel = de_df[(de_df["pval_adj"] < 0.05) & (de_df["log2fc"] > 0.25)]
print(f"\n  Melanoma-SC markers (padj<0.05, log2FC>0.25): {len(markers_mel)}")
print(f"  Top 20:")
print(markers_mel.head(20)[["gene_symbol","wilcoxon_score","log2fc","pval_adj"]].to_string(index=False))
markers_mel.to_csv("results/tables/B3_melanoma_csc_markers.csv", index=False)

# ─────────────────────────────────────────────────────────────────
# B3e — THREE-CANCER PAN-CANCER COMPARISON
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B3e — THREE-CANCER PAN-CANCER COMPARISON")
print("=" * 60)

brca_s1  = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
brca_s2  = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
gbm_df   = pd.read_csv("results/tables/B2_gbm_tumor_markers.csv")

CONSENSUS = ['B2M','CD44','EEF1A1','EPCAM','FN1','HLA-E','HMGA1','ID3',
             'KLF4','MYC','RPL22L1','S100A6','SERBP1','SERPINE2','SOX9','VIM']

brca_de   = set(brca_s1["gene_symbol"].head(200))
brca_attn = set(brca_s2["gene_symbol"].head(200))
gbm_neo   = set(gbm_df["gene_symbol"].head(200))
mel_neo   = set(markers_mel["gene_symbol"].head(200))

# Pairwise + three-way overlaps
sh_brca_mel  = brca_de   & mel_neo
sh_attn_mel  = brca_attn & mel_neo
sh_gbm_mel   = gbm_neo   & mel_neo
sh_brca_gbm  = brca_de   & gbm_neo
sh_all3_de   = brca_de   & gbm_neo & mel_neo
sh_all3_attn = brca_attn & gbm_neo & mel_neo
sh_all4      = brca_de   & brca_attn & gbm_neo & mel_neo  # all methods, all cancers

print(f"  BRCA DE   ∩ Mel: {len(sh_brca_mel):3d} genes")
print(f"  BRCA Attn ∩ Mel: {len(sh_attn_mel):3d} genes")
print(f"  GBM       ∩ Mel: {len(sh_gbm_mel):3d} genes")
print(f"  BRCA DE ∩ GBM ∩ Mel (3-cancer DE): {len(sh_all3_de):3d} genes")
print(f"  BRCA Attn ∩ GBM ∩ Mel:             {len(sh_all3_attn):3d} genes")
print(f"  ALL 4 methods (BRCA-DE+Attn+GBM+Mel): {len(sh_all4):3d} genes")

print(f"\n  ═══ 3-CANCER PAN-CSC PROGRAM (BRCA DE + GBM + Mel) ═══")
for g in sorted(sh_all3_de):
    b_lfc  = brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values
    b_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
    g_lfc  = gbm_df.loc[gbm_df["gene_symbol"]==g,"log2fc"].values
    m_lfc  = markers_mel.loc[markers_mel["gene_symbol"]==g,"log2fc"].values
    tag    = "★CON★ " if g in CONSENSUS else "      "
    attn_str = f"{b_attn[0]:.2f}" if len(b_attn) else "N/A"
    print(f"  {tag}{g:12s}  BRCA={b_lfc[0]:+.2f}  Attn={attn_str}  "
          f"GBM={g_lfc[0]:+.2f}  Mel={m_lfc[0]:+.2f}")

print(f"\n  ═══ GENEFORMER BREAST ATTENTION ∩ ALL 3 CANCERS ═══")
for g in sorted(sh_all3_attn):
    b_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
    g_lfc  = gbm_df.loc[gbm_df["gene_symbol"]==g,"log2fc"].values
    m_lfc  = markers_mel.loc[markers_mel["gene_symbol"]==g,"log2fc"].values
    rank   = list(brca_s2["gene_symbol"]).index(g) + 1
    print(f"    {g:12s}  attn={b_attn[0]:.3f}(#{rank})  GBM={g_lfc[0]:+.2f}  Mel={m_lfc[0]:+.2f}")

# Save comprehensive cross-cancer table
all_shared = sh_brca_mel | sh_attn_mel | sh_gbm_mel
pan_rows = []
for g in sorted(all_shared):
    b_lfc  = brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values
    b_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
    g_lfc  = gbm_df.loc[gbm_df["gene_symbol"]==g,"log2fc"].values
    m_lfc  = markers_mel.loc[markers_mel["gene_symbol"]==g,"log2fc"].values
    pan_rows.append({
        "gene":            g,
        "in_brca_de":     g in brca_de,
        "in_brca_attn":   g in brca_attn,
        "in_gbm_neo":     g in gbm_neo,
        "in_mel_csc":     g in mel_neo,
        "in_consensus":   g in CONSENSUS,
        "in_3cancer_de":  g in sh_all3_de,
        "in_all4":        g in sh_all4,
        "brca_de_log2fc": float(b_lfc[0])  if len(b_lfc)  else None,
        "brca_attn_score":float(b_attn[0]) if len(b_attn) else None,
        "gbm_log2fc":     float(g_lfc[0])  if len(g_lfc)  else None,
        "mel_log2fc":     float(m_lfc[0])  if len(m_lfc)  else None,
    })
pan3_df = pd.DataFrame(pan_rows)
if len(pan3_df) > 0:
    pan3_df = pan3_df.sort_values(["in_all4","in_3cancer_de","in_consensus"], ascending=False)
pan3_df.to_csv("results/tables/B3_pancancer_3cancer.csv", index=False)
print(f"\n  Saved: results/tables/B3_pancancer_3cancer.csv ({len(pan3_df)} genes)")

# ─────────────────────────────────────────────────────────────────
# B3f — FIGURES
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("B3f — GENERATING FIGURES")
print("=" * 60)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: Melanoma UMAP coloured by malignant status
mal_color = {"2": "#d62728", "1": "#1f77b4", "0": "#aaaaaa"}
colors_ct = [mal_color.get(v, "#cccccc") for v in adata.obs["malignant"]]
axes[0].scatter(adata.obsm["X_umap"][:,0], adata.obsm["X_umap"][:,1],
                c=colors_ct, s=2, alpha=0.5)
patches = [mpatches.Patch(color="#d62728",label=f"Malignant (n={( adata.obs['malignant']=='2').sum()})"),
           mpatches.Patch(color="#1f77b4",label=f"Non-malignant (n={(adata.obs['malignant']=='1').sum()})"),
           mpatches.Patch(color="#aaaaaa",label="Unresolved")]
axes[0].legend(handles=patches, fontsize=8)
axes[0].set_title("Melanoma — Malignant vs Non-malignant\n(Tirosh et al. 2016, GSE72056)")
axes[0].set_xlabel("UMAP1"); axes[0].set_ylabel("UMAP2")

# Panel 2: stemness within malignant cells
axes[1].scatter(adata_mal.obsm["X_umap"][:,0], adata_mal.obsm["X_umap"][:,1],
                c=adata_mal.obs["stemness"], cmap="RdYlGn", s=5, alpha=0.7,
                vmin=-2, vmax=2)
axes[1].set_title(f"Malignant cells: stemness score\n({adata_mal.n_obs} cells)")
axes[1].set_xlabel("UMAP1"); axes[1].set_ylabel("UMAP2")
plt.colorbar(axes[1].collections[0], ax=axes[1], label="stemness")

# Panel 3: 3-cancer overlap bar chart
cats    = ["BRCA DE\n∩Mel", "BRCA Attn\n∩Mel", "GBM∩Mel",
           "3-cancer\nDE", "Geneformer∩\n3 cancers", "All 4\nmethods"]
vals    = [len(sh_brca_mel), len(sh_attn_mel), len(sh_gbm_mel),
           len(sh_all3_de), len(sh_all3_attn), len(sh_all4)]
clrs    = ["#aaaaaa","#ff7f0e","#aaaaaa","#2ca02c","#1f77b4","#d62728"]
bars    = axes[2].bar(cats, vals, color=clrs, width=0.65)
for bar, val in zip(bars, vals):
    axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                 str(val), ha="center", fontweight="bold", fontsize=11)
axes[2].set_ylabel("Shared genes (top 200 per method)")
axes[2].set_title("Three-Cancer Pan-CSC Gene Overlap")
axes[2].tick_params(axis="x", labelsize=8)

plt.suptitle("Phase B3: Melanoma Cross-Cancer Validation", fontsize=13)
plt.tight_layout()
plt.savefig("results/figures/B3_melanoma_overview.png", dpi=120, bbox_inches="tight")
plt.close()
print("  Saved: results/figures/B3_melanoma_overview.png")

# Figure 2 — 3-cancer evidence heatmap
key_genes = sorted(sh_all3_de | sh_all3_attn)
if key_genes:
    fig, ax = plt.subplots(figsize=(max(10, len(key_genes)*0.6), 5))
    dm, row_labels = [], ["BRCA DE log2FC","BRCA Geneformer attn",
                           "GBM neo log2FC","Mel log2FC"]
    for g in key_genes:
        b_lfc  = brca_s1.loc[brca_s1["gene_symbol"]==g,"log2fc"].values
        b_attn = brca_s2.loc[brca_s2["gene_symbol"]==g,"attention_score"].values
        g_lfc  = gbm_df.loc[gbm_df["gene_symbol"]==g,"log2fc"].values
        m_lfc  = markers_mel.loc[markers_mel["gene_symbol"]==g,"log2fc"].values
        dm.append([float(b_lfc[0])  if len(b_lfc)  else 0,
                   float(b_attn[0]) if len(b_attn) else 0,
                   float(g_lfc[0])  if len(g_lfc)  else 0,
                   float(m_lfc[0])  if len(m_lfc)  else 0])
    dm   = np.array(dm)
    dmn  = (dm - dm.min(axis=0)) / (dm.max(axis=0) - dm.min(axis=0) + 1e-9)
    im   = ax.imshow(dmn.T, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(key_genes)))
    ax.set_xticklabels(key_genes, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(4)); ax.set_yticklabels(row_labels, fontsize=10)
    for i, g in enumerate(key_genes):
        if g in sh_all4:
            ax.add_patch(plt.Rectangle((i-0.5,-0.5),1,4,
                         fill=False, edgecolor="blue", lw=2.5))
        elif g in sh_all3_de:
            ax.add_patch(plt.Rectangle((i-0.5,-0.5),1,4,
                         fill=False, edgecolor="green", lw=1.5, linestyle="--"))
    plt.colorbar(im, ax=ax, label="normalized score")
    from matplotlib.lines import Line2D
    leg = [Line2D([0],[0], color="blue",  lw=2.5, label="All 4 methods"),
           Line2D([0],[0], color="green", lw=1.5, linestyle="--", label="3-cancer DE")]
    ax.legend(handles=leg, fontsize=9, loc="upper right")
    ax.set_title("Three-Cancer CSC Evidence Matrix\n"
                 "(BRCA breast + GBM glioblastoma + Melanoma)")
    plt.tight_layout()
    plt.savefig("results/figures/B3_3cancer_heatmap.png", dpi=120, bbox_inches="tight")
    plt.close()
    print("  Saved: results/figures/B3_3cancer_heatmap.png")

print("\n" + "=" * 60)
print("PHASE B3 COMPLETE — THREE-CANCER VALIDATION")
print("=" * 60)
print(f"  Melanoma SC markers: {len(markers_mel)}")
print(f"  3-cancer pan-CSC (DE): {len(sh_all3_de)} — {sorted(sh_all3_de)}")
print(f"  Geneformer ∩ 3 cancers: {len(sh_all3_attn)} — {sorted(sh_all3_attn)}")
print(f"  All 4 methods: {len(sh_all4)} — {sorted(sh_all4)}")
