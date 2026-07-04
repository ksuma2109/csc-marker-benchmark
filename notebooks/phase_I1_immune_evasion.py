# Phase I1 — CSC state and the tumor immune microenvironment (extension B)
#
# Tests whether cancer stem cell (CSC) state is associated with immune evasion,
# using the immune + malignant compartments of the breast atlas (GSE176078;
# 100,064 cells, 26 patients, 35k T-cells). No wet lab.
#
# Two questions:
#   I1a (cell-level, cancer-epithelial): are CSC-high malignant cells
#        intrinsically immune-evasive? Correlate per-cell stemness with
#        (i) MHC-I / antigen-presentation machinery, (ii) checkpoint ligands.
#   I1b (patient-level, 26 patients): do high-CSC tumors have less T-cell
#        infiltration and more T-cell exhaustion?
#
# Output: results/tables/I1_immune_evasion.csv
#         results/tables/I1_patient_level.csv
#         results/figures/I1_immune_evasion.png

import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

print("=" * 64)
print("I1 — CSC STATE AND TUMOR IMMUNE MICROENVIRONMENT")
print("=" * 64)

# gene programs (filtered to present genes at use)
# NB: exclude B2M and HLA-E from the MHC-I score — they are themselves CSC
# consensus genes, so including them would create a shared-gene confound with
# the stemness score. Use only classical antigen-presentation genes.
MHC_APM   = ["HLA-A","HLA-B","HLA-C","TAP1","TAP2","NLRC5"]                        # antigen presentation (classical)
CKPT_LIG  = ["CD274","PDCD1LG2","CD47"]                                            # checkpoint / evasion ligands
EXHAUST   = ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX"]                        # T-cell exhaustion
CYTOTOX   = ["GZMB","GZMA","PRF1","IFNG","NKG7","GZMK"]                            # T-cell effector

def score(ad, genes, name):
    g = [x for x in genes if x in ad.var_names]
    sc.tl.score_genes(ad, g, score_name=name, ctrl_size=50)
    return g

# ─────────────────────────────────────────────────────────────────
# I1a — CELL-LEVEL: are CSC-high cancer cells immune-evasive?
# ─────────────────────────────────────────────────────────────────
print("\nI1a — cancer-epithelial cells (per-cell)...")
epi = sc.read_h5ad("data/processed/brca_A5_csc_scored.h5ad")
score(epi, MHC_APM,  "score_APM")          # antigen-presentation machinery
score(epi, CKPT_LIG, "score_checkpoint")   # checkpoint ligands
stem = epi.obs["stemness_composite"].values

rows = []
for label, col in [("Antigen-presentation (MHC-I)", "score_APM"),
                   ("Checkpoint ligands (PD-L1/L2/CD47)", "score_checkpoint")]:
    r, p = spearmanr(stem, epi.obs[col].values)
    # CSC-high vs CSC-low contrast
    hi = epi.obs.loc[epi.obs["csc_label"]=="csc_high", col]
    lo = epi.obs.loc[epi.obs["csc_label"]=="csc_low",  col]
    u, pmw = mannwhitneyu(hi, lo, alternative="two-sided")
    rows.append({"program": label, "spearman_r_vs_stemness": round(r,3), "spearman_p": p,
                 "mean_csc_high": round(float(hi.mean()),3), "mean_csc_low": round(float(lo.mean()),3),
                 "csc_high_minus_low": round(float(hi.mean()-lo.mean()),3), "mw_p": pmw})
    print(f"  {label}: Spearman r={r:+.3f} (p={p:.1e}); "
          f"CSC-high {hi.mean():+.3f} vs CSC-low {lo.mean():+.3f} (MW p={pmw:.1e})")
cell_df = pd.DataFrame(rows)
cell_df.to_csv("results/tables/I1_immune_evasion.csv", index=False)

# ─────────────────────────────────────────────────────────────────
# I1b — PATIENT-LEVEL: CSC vs T-cell infiltration & exhaustion
# ─────────────────────────────────────────────────────────────────
print("\nI1b — patient-level (26 patients)...")
adata = sc.read_h5ad("data/processed/brca_A4_annotated.h5ad")
# per-patient mean stemness (from cancer-epithelial subset)
epi.obs["patient"] = epi.obs["orig.ident"]
stem_by_pt = epi.obs.groupby("patient")["stemness_composite"].mean()

# T-cell subset + exhaustion / cytotoxic scores
tcells = adata[adata.obs["celltype_major"] == "T-cells"].copy()
cd8 = tcells[tcells.obs["celltype_minor"].astype(str).str.contains("CD8", na=False)].copy()
score(cd8, EXHAUST, "score_exhaustion")
score(cd8, CYTOTOX, "score_cytotoxic")

pt_rows = []
for pt in stem_by_pt.index:
    total = (adata.obs["orig.ident"] == pt).sum()
    n_t   = ((adata.obs["orig.ident"] == pt) & (adata.obs["celltype_major"]=="T-cells")).sum()
    cd8_pt = cd8.obs[cd8.obs["orig.ident"] == pt]
    if total < 50 or len(cd8_pt) < 10:
        continue
    subtype = epi.obs.loc[epi.obs["patient"] == pt, "subtype"].mode()
    pt_rows.append({
        "patient": pt,
        "subtype": subtype.iloc[0] if len(subtype) else "NA",
        "mean_stemness": float(stem_by_pt[pt]),
        "tcell_fraction": n_t / total,
        "n_CD8": len(cd8_pt),
        "cd8_exhaustion": float(cd8_pt["score_exhaustion"].mean()),
        "cd8_cytotoxic":  float(cd8_pt["score_cytotoxic"].mean()),
    })
pt_df = pd.DataFrame(pt_rows)
pt_df.to_csv("results/tables/I1_patient_level.csv", index=False)
print(f"  Patients analysed: {len(pt_df)}  |  subtypes: {pt_df['subtype'].value_counts().to_dict()}")
print(f"  mean stemness by subtype: "
      + ", ".join(f"{k}={v:.2f}" for k,v in pt_df.groupby('subtype')['mean_stemness'].mean().items()))

corr = {}
for col, lab in [("tcell_fraction","T-cell infiltration"),
                 ("cd8_exhaustion","CD8 exhaustion"),
                 ("cd8_cytotoxic","CD8 cytotoxicity")]:
    r, p = spearmanr(pt_df["mean_stemness"], pt_df[col])
    corr[col] = (r, p)
    print(f"  stemness vs {lab:22s}: Spearman r={r:+.3f}  p={p:.3f}")

# ─────────────────────────────────────────────────────────────────
# FIGURE
# ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(20, 5))

# cell-level: stemness vs APM and checkpoint (hexbin for density)
ax = axes[0]
ax.hexbin(stem, epi.obs["score_APM"], gridsize=40, cmap="Blues", mincnt=1)
r0 = cell_df.iloc[0]["spearman_r_vs_stemness"]
ax.set_xlabel("stemness"); ax.set_ylabel("antigen-presentation (MHC-I) score")
ax.set_title(f"CSC vs antigen presentation\nSpearman r={r0:+.3f}")

ax = axes[1]
ax.hexbin(stem, epi.obs["score_checkpoint"], gridsize=40, cmap="Reds", mincnt=1)
r1 = cell_df.iloc[1]["spearman_r_vs_stemness"]
ax.set_xlabel("stemness"); ax.set_ylabel("checkpoint-ligand score")
ax.set_title(f"CSC vs checkpoint ligands\nSpearman r={r1:+.3f}")

# patient-level: stemness vs T-cell fraction
ax = axes[2]
ax.scatter(pt_df["mean_stemness"], pt_df["tcell_fraction"], s=40, color="#2ca02c", alpha=0.8)
r, p = corr["tcell_fraction"]
ax.set_xlabel("patient mean stemness"); ax.set_ylabel("T-cell fraction")
ax.set_title(f"CSC vs T-cell infiltration\nSpearman r={r:+.3f}, p={p:.3f}")

# patient-level: stemness vs CD8 exhaustion
ax = axes[3]
ax.scatter(pt_df["mean_stemness"], pt_df["cd8_exhaustion"], s=40, color="#d62728", alpha=0.8)
r, p = corr["cd8_exhaustion"]
ax.set_xlabel("patient mean stemness"); ax.set_ylabel("CD8 exhaustion score")
ax.set_title(f"CSC vs CD8 exhaustion\nSpearman r={r:+.3f}, p={p:.3f}")

plt.suptitle("Figure. CSC state and the tumor immune microenvironment (GSE176078)", fontsize=13)
plt.tight_layout()
plt.savefig("results/figures/I1_immune_evasion.png", dpi=130, bbox_inches="tight")
plt.close()

print("\n" + "=" * 64)
print("PHASE I1 COMPLETE")
print("=" * 64)
print("  Saved: results/tables/I1_immune_evasion.csv, I1_patient_level.csv")
print("  Saved: results/figures/I1_immune_evasion.png")
