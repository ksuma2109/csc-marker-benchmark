# Phase I8 — Discovery-vs-validation figure for the CSC-immune (second) manuscript
#
# Assembles the key result: TNBC CSCs upregulate CD47 and retain MHC-I, a
# SUBTYPE-SPECIFIC phenotype that REPLICATES in an independent cohort.
#   Panel A  Subtype specificity (discovery, I5): CD47/MHC-I/HLA-E/PD-L1 vs
#            stemness within ER+/TNBC/HER2+ — the effect is TNBC-specific.
#   Panel B  Discovery vs independent validation (I5 TNBC vs I7 Pal 2021):
#            replication of the four programs, with per-tumour validation points.
#
# Output: manuscript2/figures/Figure1.png

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.makedirs("manuscript2/figures", exist_ok=True)

i5 = pd.read_csv("results/tables/I5_subtype_stratified.csv")
i7 = pd.read_csv("results/tables/I7_validation_cohort.csv")
i9 = pd.read_csv("results/tables/I9_validation2_gao.csv")

PROGRAMS = ["MHC-I antigen presentation", "CD47 ('don't eat me')", "HLA-E (NKG2A)", "CD274 (PD-L1)"]
SHORT    = ["MHC-I", "CD47", "HLA-E", "PD-L1"]
i7_cols  = {"MHC-I": "MHC-I_r", "CD47": "CD47_r", "HLA-E": "HLA-E_r", "PD-L1": "PD-L1_r"}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 5.5))

# ── Panel A: subtype specificity (discovery / I5) ──
subs = ["ER+", "TNBC", "HER2+"]
colors = {"ER+": "#1f77b4", "TNBC": "#d62728", "HER2+": "#2ca02c"}
x = np.arange(len(SHORT)); w = 0.25
for i, st in enumerate(subs):
    vals = []
    for prog in PROGRAMS:
        row = i5[(i5["program"] == prog) & (i5["subtype"] == st)]
        vals.append(row["spearman_r"].values[0] if len(row) else np.nan)
    axA.bar(x + (i-1)*w, vals, w, label=st, color=colors[st])
axA.axhline(0, color="k", lw=0.6)
axA.set_xticks(x); axA.set_xticklabels(SHORT)
axA.set_ylabel("Spearman r (program vs stemness)")
axA.set_title("A  Subtype specificity (discovery, Wu 2021)\n"
              "CSC innate-evasion phenotype is TNBC-specific", fontsize=11, loc="left")
axA.legend(title="subtype", fontsize=8)

# ── Panel B: discovery TNBC vs independent validation (Pal 2021) ──
disc = []
for prog in PROGRAMS:
    row = i5[(i5["program"] == prog) & (i5["subtype"] == "TNBC")]
    disc.append(row["spearman_r"].values[0] if len(row) else np.nan)
pal_mean = [i7[i7_cols[s]].mean() for s in SHORT]
gao_mean = [i9[i7_cols[s]].mean() for s in SHORT]
pal_pts  = {s: i7[i7_cols[s]].dropna().values for s in SHORT}
gao_pts  = {s: i9[i7_cols[s]].dropna().values for s in SHORT}

xb = np.arange(len(SHORT)); wb = 0.27
axB.bar(xb - wb, disc, wb, label="Discovery (Wu 2021)", color="#d62728", alpha=0.85)
axB.bar(xb,      pal_mean, wb, label="Validation 1 (Pal 2021)", color="#ff7f0e", alpha=0.85)
axB.bar(xb + wb, gao_mean, wb, label="Validation 2 (Gao 2020)", color="#9467bd", alpha=0.85)
for i, s in enumerate(SHORT):
    axB.scatter(np.full(len(pal_pts[s]), xb[i]),      pal_pts[s], s=14, color="k", zorder=5, alpha=0.6)
    axB.scatter(np.full(len(gao_pts[s]), xb[i] + wb), gao_pts[s], s=14, color="k", zorder=5, alpha=0.6)
axB.axhline(0, color="k", lw=0.6)
axB.set_xticks(xb); axB.set_xticklabels(SHORT)
axB.set_ylabel("Spearman r (program vs stemness, malignant cells)")
axB.set_title("B  Replication across two independent cohorts\n"
              "CD47 & MHC-I replicate; PD-L1 stays ~0", fontsize=11, loc="left")
axB.legend(fontsize=7.5)

fig.suptitle("Figure 1. TNBC cancer stem cells upregulate CD47 and retain antigen presentation "
             "(subtype-specific, independently replicated)", fontsize=12)
plt.tight_layout()
fig.savefig("manuscript2/figures/Figure1.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved manuscript2/figures/Figure1.png")
print("\nDiscovery TNBC r:", dict(zip(SHORT, [round(d,3) for d in disc])))
print("Pal 2021 mean r: ", dict(zip(SHORT, [round(v,3) for v in pal_mean])))
print("Gao 2020 mean r: ", dict(zip(SHORT, [round(v,3) for v in gao_mean])))
