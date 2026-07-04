# Clean Figure 1 — regenerated as four single panels directly from data
# (replaces the tiled-composite version). Publication styling: uniform fonts,
# panel labels, single legend per panel.
#
# Panels:
#   A  Cell-type annotated UMAP (all 100k cells)
#   B  Stemness score on cancer-epithelial cells (UMAP)
#   C  DE volcano (CSC-high vs CSC-low)
#   D  Pseudotime on cancer-epithelial cells (UMAP)
#
# Output: manuscript/figures/Figure1.png

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

OUT = "manuscript/figures"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "figure.dpi": 150})

fig, axes = plt.subplots(2, 2, figsize=(13, 12))
fig.subplots_adjust(hspace=0.22, wspace=0.18)

def label(ax, letter):
    ax.text(-0.06, 1.06, letter, transform=ax.transAxes,
            fontsize=20, fontweight="bold", va="top", ha="right")

# ── Panel A — cell-type annotated UMAP ─────────────────────────
print("Panel A: annotated UMAP...")
a4 = sc.read_h5ad("data/processed/brca_A4_annotated.h5ad")
um = a4.obsm["X_umap"]
cts = a4.obs["celltype_major"].astype(str).values
palette = {
    "Cancer Epithelial":"#d62728","T-cells":"#1f77b4","Myeloid":"#2ca02c",
    "CAFs":"#9467bd","Endothelial":"#8c564b","PVL":"#e377c2",
    "B-cells":"#7f7f7f","Plasmablasts":"#bcbd22","Normal Epithelial":"#17becf",
}
ax = axes[0,0]
order = list(palette.keys())
for ct in order:
    m = cts == ct
    ax.scatter(um[m,0], um[m,1], s=1.2, c=palette[ct], alpha=0.55, rasterized=True)
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.set_title(f"Cell-type annotation (n = {a4.n_obs:,})")
ax.set_xticks([]); ax.set_yticks([])
handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor=palette[c],
                  markersize=6, label=c) for c in order]
ax.legend(handles=handles, fontsize=6.5, loc="upper right", framealpha=0.9, ncol=1)
label(ax, "A")
del a4

# ── Panel B — stemness score on cancer-epithelial cells ────────
print("Panel B: stemness UMAP...")
a5 = sc.read_h5ad("data/processed/brca_A5_csc_scored.h5ad")
um = a5.obsm["X_umap"]
sco = a5.obs["stemness_composite"].values
ax = axes[0,1]
scat = ax.scatter(um[:,0], um[:,1], s=2.5, c=sco, cmap="RdYlBu_r",
                  alpha=0.8, rasterized=True, vmin=np.percentile(sco,2),
                  vmax=np.percentile(sco,98))
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.set_title(f"Stemness score, cancer-epithelial (n = {a5.n_obs:,})")
ax.set_xticks([]); ax.set_yticks([])
cb = fig.colorbar(scat, ax=ax, fraction=0.046, pad=0.02)
cb.set_label("stemness score", fontsize=8)
label(ax, "B")
del a5

# ── Panel C — DE volcano ───────────────────────────────────────
print("Panel C: volcano...")
de = pd.read_csv("results/tables/A5_all_DE_results.csv") if os.path.exists(
        "results/tables/A5_all_DE_results.csv") else pd.read_csv("results/tables/A5_csc_markers_DE.csv")
# expected columns: gene_symbol, log2fc, pval_adj
lfc_col = "log2fc" if "log2fc" in de.columns else [c for c in de.columns if "fc" in c.lower()][0]
p_col   = "pval_adj" if "pval_adj" in de.columns else [c for c in de.columns if "adj" in c.lower()][0]
de = de.dropna(subset=[lfc_col, p_col]).copy()
de["neglog10p"] = -np.log10(de[p_col].clip(lower=1e-300))
ax = axes[1,0]
sig_up   = (de[p_col] < 0.05) & (de[lfc_col] > 0.5)
sig_dn   = (de[p_col] < 0.05) & (de[lfc_col] < -0.5)
ns       = ~(sig_up | sig_dn)
ax.scatter(de.loc[ns,lfc_col], de.loc[ns,"neglog10p"], s=4, c="#cccccc", alpha=0.5, rasterized=True)
ax.scatter(de.loc[sig_up,lfc_col], de.loc[sig_up,"neglog10p"], s=6, c="#d62728", alpha=0.6, rasterized=True)
ax.scatter(de.loc[sig_dn,lfc_col], de.loc[sig_dn,"neglog10p"], s=6, c="#1f77b4", alpha=0.6, rasterized=True)
# label top CSC genes
gcol = "gene_symbol" if "gene_symbol" in de.columns else de.columns[0]
for g in ["MYC","CD44","VIM","SOX9","KLF4","SERPINE2"]:
    r = de[de[gcol]==g]
    if len(r):
        ax.annotate(g, (r[lfc_col].values[0], r["neglog10p"].values[0]),
                    fontsize=7, fontweight="bold")
ax.axvline(0, color="k", lw=0.4); ax.axhline(-np.log10(0.05), color="gray", ls="--", lw=0.6)
ax.set_xlabel("log$_2$ fold change (CSC-high vs CSC-low)")
ax.set_ylabel("-log$_{10}$ adjusted p")
ax.set_title("Stage 1 differential expression")
label(ax, "C")

# ── Panel D — pseudotime ───────────────────────────────────────
print("Panel D: pseudotime UMAP...")
a6 = sc.read_h5ad("data/processed/brca_A6_trajectory.h5ad")
um = a6.obsm["X_umap"]
pt = a6.obs["dpt_pseudotime"].values
ax = axes[1,1]
sca = ax.scatter(um[:,0], um[:,1], s=2.5, c=pt, cmap="viridis", alpha=0.8, rasterized=True)
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
ax.set_title("Pseudotime trajectory (stem → differentiated)")
ax.set_xticks([]); ax.set_yticks([])
cb = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
cb.set_label("pseudotime", fontsize=8)
label(ax, "D")
del a6

fig.suptitle("Figure 1. Single-cell identification of breast cancer stem cells (GSE176078)",
             fontsize=13, y=0.965)
fig.savefig(f"{OUT}/Figure1.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved {OUT}/Figure1.png")
