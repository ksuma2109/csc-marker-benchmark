# Phase F2 — Per-Cell Protein-Anchored Benchmark: DE vs Geneformer
#
# The strongest functional anchor in the project. Uses CITE-seq surface
# PROTEIN (antibody-derived tags, ADT) from the SAME breast atlas
# (GSE176078 / SCP1039) to define CSCs the original way the field did —
# CD44+CD24-/low at the PROTEIN level (Al-Hajj et al. 2003) — independent
# of mRNA. Then asks which mRNA gene ranking (Stage 1 DE vs Stage 2
# Geneformer attention) better recovers the protein-defined CSC state,
# per cell.
#
# Only sample CID4515 is usable: it is the one CITE-seq case that has BOTH
# cancer-epithelial cells AND CD44/CD24/EPCAM in its antibody panel.
#   (3838/3946/4040 carry the full panel but contain no malignant cells;
#    4378N is normal tissue.)
#
# Pipeline:
#   F2a  Load ADT (protein) for CID4515, CLR-normalize, gate CD44+CD24-low
#   F2b  Load matched RNA (cancer-epithelial cells), log-normalize
#   F2c  Score each cell with DE list and Geneformer list (scanpy score_genes)
#   F2d  AUROC / AP of each score vs the protein CSC label (per cell)
#   F2e  Figure + summary
#
# Output: results/tables/F2_protein_benchmark.csv
#         results/figures/F2_protein_benchmark.png

import os, gzip, warnings
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

CITE = "data/raw/GSE176078_CITE/CITE/4515/4515_CITE.miniatlas/umi_count"
RNA  = "data/raw/Wu_etal_2021_BRCA_scRNASeq"
SAMPLE = "CID4515"
os.makedirs("results/tables",  exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

print("=" * 64)
print("F2 — PER-CELL PROTEIN-ANCHORED BENCHMARK (CITE-seq, CID4515)")
print("=" * 64)

# ─────────────────────────────────────────────────────────────────
# F2a — LOAD ADT PROTEIN, CLR-NORMALIZE, GATE CD44+CD24-low
# ─────────────────────────────────────────────────────────────────
print("\nF2a — Loading ADT (surface protein)...")
adt_mat = sio.mmread(os.path.join(CITE, "matrix.mtx.gz")).tocsr()   # features × cells
adt_feat = pd.read_csv(os.path.join(CITE, "features.tsv.gz"), header=None)[0].values
adt_bc   = pd.read_csv(os.path.join(CITE, "barcodes.tsv.gz"), header=None)[0].values
adt = pd.DataFrame(adt_mat.toarray(), index=adt_feat, columns=adt_bc)
# feature names look like "CD44-AATCC..." -> take protein symbol before first "-"
adt.index = [f.split("-")[0] for f in adt.index]
adt = adt[~adt.index.duplicated(keep="first")]
print(f"  ADT: {adt.shape[0]} proteins × {adt.shape[1]} cells")

# CLR (centered log-ratio) normalization across proteins, per cell — standard for ADT
def clr(col):
    x = np.log1p(col)
    return x - x.mean()
adt_clr = adt.apply(clr, axis=0)   # per cell

cd44 = adt_clr.loc["CD44"]
cd24 = adt_clr.loc["CD24"]
print(f"  CD44 protein: median CLR={cd44.median():.2f}")
print(f"  CD24 protein: median CLR={cd24.median():.2f}")

# ─────────────────────────────────────────────────────────────────
# F2b — LOAD MATCHED RNA (cancer-epithelial cells)
# ─────────────────────────────────────────────────────────────────
print("\nF2b — Loading matched RNA (cancer epithelial)...")
meta = pd.read_csv(os.path.join(RNA, "metadata.csv"), index_col=0)
meta_s = meta[meta["orig.ident"] == SAMPLE].copy()
meta_s["bc"] = [x.split("_", 1)[1] for x in meta_s.index]   # strip CID4515_ prefix

# RNA sparse matrix (genes × cells)
rna_genes = pd.read_csv(os.path.join(RNA, "count_matrix_genes.tsv"), header=None)[0].values
rna_bc    = pd.read_csv(os.path.join(RNA, "count_matrix_barcodes.tsv"), header=None)[0].values
print(f"  RNA global matrix: {len(rna_genes)} genes × {len(rna_bc)} cells (loading subset)...")

# Build full AnnData once, subset to this sample's cancer-epithelial cells
adata = sc.read_mtx(os.path.join(RNA, "count_matrix_sparse.mtx")).T  # cells × genes
adata.var_names = rna_genes
adata.obs_names = rna_bc
# keep this sample's cancer epithelial cells
epi_ids = meta_s[meta_s["celltype_major"] == "Cancer Epithelial"].index
adata = adata[adata.obs_names.isin(epi_ids)].copy()
adata.obs = adata.obs.join(meta_s)
adata.obs["bc"] = [x.split("_", 1)[1] for x in adata.obs_names]
print(f"  Cancer-epithelial cells (RNA): {adata.n_obs}")

# match to ADT by barcode
adata = adata[adata.obs["bc"].isin(set(adt.columns))].copy()
print(f"  Cells with both RNA + ADT: {adata.n_obs}")

# normalize RNA
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# attach protein values (aligned by barcode)
adata.obs["CD44_prot"] = cd44.reindex(adata.obs["bc"]).values
adata.obs["CD24_prot"] = cd24.reindex(adata.obs["bc"]).values

# ─────────────────────────────────────────────────────────────────
# DEFINE PROTEIN CSC LABEL — CD44+CD24-/low (Al-Hajj gate)
# ─────────────────────────────────────────────────────────────────
print("\nDefining protein CSC gate (CD44+ / CD24-low)...")
cd44_hi = adata.obs["CD44_prot"] >= adata.obs["CD44_prot"].median()
cd24_lo = adata.obs["CD24_prot"] <= adata.obs["CD24_prot"].median()
adata.obs["protein_CSC"] = (cd44_hi & cd24_lo).astype(int)
n_csc = int(adata.obs["protein_CSC"].sum())
print(f"  CD44+CD24-low CSCs: {n_csc} / {adata.n_obs} "
      f"({100*n_csc/adata.n_obs:.1f}%)")

# ─────────────────────────────────────────────────────────────────
# F2c — SCORE CELLS WITH DE AND GENEFORMER GENE LISTS
# ─────────────────────────────────────────────────────────────────
print("\nF2c — Scoring cells with each method's gene list...")
s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])

TOPN = 100
de_genes = [g for g in s1["gene_symbol"].head(TOPN) if g in adata.var_names]
gf_genes = [g for g in s2["gene_symbol"].head(TOPN) if g in adata.var_names]
print(f"  DE genes in data: {len(de_genes)}/{TOPN}")
print(f"  Geneformer genes in data: {len(gf_genes)}/{TOPN}")

sc.tl.score_genes(adata, de_genes, score_name="DE_score", ctrl_size=50)
sc.tl.score_genes(adata, gf_genes, score_name="GF_score", ctrl_size=50)

# ─────────────────────────────────────────────────────────────────
# F2d — AUROC / AP vs PROTEIN CSC LABEL
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 64)
print("F2d — DISCRIMINATION OF PROTEIN-DEFINED CSCs (per cell)")
print("=" * 64)

y = adata.obs["protein_CSC"].values
results = []
for method, col in [("Stage1_DE", "DE_score"), ("Stage2_Geneformer", "GF_score")]:
    score = adata.obs[col].values
    auroc = roc_auc_score(y, score)
    ap    = average_precision_score(y, score)
    # also CD44 protein alone as a positive control
    results.append({"method": method, "AUROC": round(auroc,3),
                    "avg_precision": round(ap,3), "n_cells": adata.n_obs,
                    "n_csc": n_csc})
    print(f"  {method:18s}  AUROC={auroc:.3f}  AP={ap:.3f}")

# positive control: does CD44 mRNA predict CD44+CD24- protein gate?
for ctrl_gene in ["CD44", "CD24"]:
    if ctrl_gene in adata.var_names:
        g = np.asarray(adata[:, ctrl_gene].X.todense()).flatten()
        a = roc_auc_score(y, g if ctrl_gene=="CD44" else -g)
        print(f"  [control] {ctrl_gene} mRNA {'(+)' if ctrl_gene=='CD44' else '(-)'} → AUROC={a:.3f}")

res = pd.DataFrame(results)
res.to_csv("results/tables/F2_protein_benchmark.csv", index=False)
print("\n  Saved: results/tables/F2_protein_benchmark.csv")

winner = res.loc[res["AUROC"].idxmax(), "method"]
print(f"\n  Better per-cell discriminator of protein CSCs: {winner}")

# ─────────────────────────────────────────────────────────────────
# F2e — FIGURE
# ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: CD44 vs CD24 protein scatter, CSC gate highlighted
ax = axes[0]
csc = adata.obs["protein_CSC"] == 1
ax.scatter(adata.obs.loc[~csc,"CD24_prot"], adata.obs.loc[~csc,"CD44_prot"],
           s=8, c="#cccccc", label="non-CSC", alpha=0.6)
ax.scatter(adata.obs.loc[csc,"CD24_prot"], adata.obs.loc[csc,"CD44_prot"],
           s=10, c="#d62728", label=f"CD44+CD24- CSC (n={n_csc})", alpha=0.8)
ax.axhline(adata.obs["CD44_prot"].median(), color="k", ls="--", lw=0.7)
ax.axvline(adata.obs["CD24_prot"].median(), color="k", ls="--", lw=0.7)
ax.set_xlabel("CD24 protein (CLR)"); ax.set_ylabel("CD44 protein (CLR)")
ax.set_title(f"{SAMPLE}: protein CSC gate\n(Al-Hajj CD44+CD24-/low)")
ax.legend(fontsize=8)

# Panel 2: ROC curves
ax = axes[1]
for method, col, color in [("Stage1 DE","DE_score","#d62728"),
                           ("Stage2 Geneformer","GF_score","#1f77b4")]:
    fpr, tpr, _ = roc_curve(y, adata.obs[col].values)
    auc = roc_auc_score(y, adata.obs[col].values)
    ax.plot(fpr, tpr, color=color, lw=2, label=f"{method} (AUROC={auc:.3f})")
ax.plot([0,1],[0,1],"k--",lw=0.8, label="chance")
ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
ax.set_title("Per-cell discrimination of\nprotein-defined CSCs")
ax.legend(fontsize=8, loc="lower right")

# Panel 3: score distributions by protein CSC status
ax = axes[2]
data_de = [adata.obs.loc[~csc,"DE_score"], adata.obs.loc[csc,"DE_score"]]
data_gf = [adata.obs.loc[~csc,"GF_score"], adata.obs.loc[csc,"GF_score"]]
bp1 = ax.boxplot(data_de, positions=[1,2], widths=0.35, patch_artist=True,
                 boxprops=dict(facecolor="#d62728", alpha=0.6))
bp2 = ax.boxplot(data_gf, positions=[4,5], widths=0.35, patch_artist=True,
                 boxprops=dict(facecolor="#1f77b4", alpha=0.6))
ax.set_xticks([1,2,4,5])
ax.set_xticklabels(["non-CSC","CSC","non-CSC","CSC"], fontsize=8)
ax.text(1.5, ax.get_ylim()[1]*0.95, "Stage1 DE", ha="center", color="#d62728", fontweight="bold")
ax.text(4.5, ax.get_ylim()[1]*0.95, "Geneformer", ha="center", color="#1f77b4", fontweight="bold")
ax.set_ylabel("Module score")
ax.set_title("Stemness score by\nprotein CSC status")

plt.suptitle(f"Phase F2: Protein-Anchored CSC Benchmark — {SAMPLE} CITE-seq "
             f"({adata.n_obs} cancer-epithelial cells)", fontsize=12)
plt.tight_layout()
plt.savefig("results/figures/F2_protein_benchmark.png", dpi=120, bbox_inches="tight")
plt.close()

print("  Saved: results/figures/F2_protein_benchmark.png")
print("\n" + "=" * 64)
print("PHASE F2 COMPLETE")
print("=" * 64)
print(f"  {adata.n_obs} cancer-epithelial cells, {n_csc} protein-defined CSCs")
for _, r in res.iterrows():
    print(f"  {r['method']:18s}  AUROC={r['AUROC']}  AP={r['avg_precision']}")
