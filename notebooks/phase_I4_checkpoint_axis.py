# Phase I4 — Which immune-checkpoint / evasion axis do CSCs most upregulate?
#
# Follow-up to I1: the CSC-high cancer cells weakly upregulated a *combined*
# checkpoint-ligand score. Here we break it down per-ligand across an expanded
# panel of tumour-intrinsic immune-evasion surface ligands, to nominate the
# specific combination partner (e.g. anti-PD-L1 vs anti-CD47 vs anti-B7-H3)
# for a CSC-directed immunotherapy.
#
# For each ligand: per-cell Spearman vs stemness, and CSC-high vs CSC-low mean
# difference (log-normalized), ranked by effect size.
#
# Output: results/tables/I4_checkpoint_axis.csv
#         results/figures/I4_checkpoint_axis.png

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

# Tumour-intrinsic immune-evasion ligands (surface unless noted).
# label: (gene, pathway / therapeutic axis, note)
LIGANDS = {
    "CD274 (PD-L1)":   ("CD274",    "PD-1 axis",        "anti-PD-L1 (atezolizumab)"),
    "PDCD1LG2 (PD-L2)":("PDCD1LG2", "PD-1 axis",        "anti-PD-1"),
    "CD47":            ("CD47",     "'don't eat me'",   "anti-CD47 (magrolimab)"),
    "CD276 (B7-H3)":   ("CD276",    "B7 family",        "anti-B7-H3 (enoblituzumab)"),
    "VTCN1 (B7-H4)":   ("VTCN1",    "B7 family",        "anti-B7-H4"),
    "LGALS9 (Gal-9)":  ("LGALS9",   "TIM-3 axis",       "anti-Gal-9 / anti-TIM-3"),
    "NT5E (CD73)":     ("NT5E",     "adenosine",        "anti-CD73"),
    "ENTPD1 (CD39)":   ("ENTPD1",   "adenosine",        "anti-CD39"),
    "HLA-E":           ("HLA-E",    "NKG2A (NK/T inhib)","anti-NKG2A (monalizumab); NB CSC-consensus gene"),
    "HLA-G":           ("HLA-G",    "ILT2/4 tolerance", "anti-HLA-G"),
    "CD24":            ("CD24",     "'don't eat me' (Siglec-10)", "CSCs are CD24-low by definition — expect NEGATIVE"),
}

print("=" * 64)
print("I4 — CSC UPREGULATION OF IMMUNE-EVASION LIGANDS (per ligand)")
print("=" * 64)

epi = sc.read_h5ad("data/processed/brca_A5_csc_scored.h5ad")
stem = epi.obs["stemness_composite"].values
hi = epi.obs["csc_label"] == "csc_high"
lo = epi.obs["csc_label"] == "csc_low"

def expr(gene):
    x = epi[:, gene].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()

rows = []
for label, (gene, axis, drug) in LIGANDS.items():
    if gene not in epi.var_names:
        print(f"  {label}: gene not present, skipped")
        continue
    e = expr(gene)
    r, p = spearmanr(stem, e)
    ehi, elo = e[hi.values], e[lo.values]
    delta = float(ehi.mean() - elo.mean())          # log-normalized mean difference
    u, pmw = mannwhitneyu(ehi, elo, alternative="two-sided")
    pct_hi = float((ehi > 0).mean()); pct_lo = float((elo > 0).mean())
    rows.append({"ligand": label, "gene": gene, "axis": axis,
                 "spearman_r_vs_stemness": round(r, 3), "spearman_p": p,
                 "csc_high_minus_low": round(delta, 3), "mw_p": pmw,
                 "pct_expr_high": round(pct_hi, 3), "pct_expr_low": round(pct_lo, 3),
                 "therapeutic": drug})

res = pd.DataFrame(rows).sort_values("csc_high_minus_low", ascending=False)
res.to_csv("results/tables/I4_checkpoint_axis.csv", index=False)

print(f"\n  {'ligand':18s} {'axis':16s} {'CSChi-lo':>9s} {'r':>7s} {'%hi':>5s} {'%lo':>5s}")
for _, r in res.iterrows():
    print(f"  {r['ligand']:18s} {r['axis']:16s} {r['csc_high_minus_low']:>+9.3f} "
          f"{r['spearman_r_vs_stemness']:>+7.3f} {r['pct_expr_high']*100:>4.0f}% {r['pct_expr_low']*100:>4.0f}%")

# nominate top actionable surface ligand (exclude HLA-E confound and CD24)
actionable = res[~res["gene"].isin(["HLA-E", "CD24"])]
top = actionable[actionable["csc_high_minus_low"] > 0].head(1)
print("\n" + "=" * 64)
if len(top):
    t = top.iloc[0]
    print(f"  TOP CSC-upregulated actionable evasion ligand: {t['ligand']}  "
          f"(Δ={t['csc_high_minus_low']:+.3f}, r={t['spearman_r_vs_stemness']:+.3f})")
    print(f"  Nominated combination partner: {t['therapeutic']}")
# CD24 sanity check
cd24 = res[res["gene"] == "CD24"]
if len(cd24):
    print(f"  Sanity check — CD24 (CSC-low by definition): Δ={cd24.iloc[0]['csc_high_minus_low']:+.3f} "
          f"(expected negative)")

# ── figure ──
fig, ax = plt.subplots(figsize=(9, 6))
plot = res.iloc[::-1]
colors = ["#999" if g in ("HLA-E", "CD24") else ("#d62728" if d > 0 else "#1f77b4")
          for g, d in zip(plot["gene"], plot["csc_high_minus_low"])]
ax.barh(plot["ligand"], plot["csc_high_minus_low"], color=colors)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("CSC-high − CSC-low mean expression (log-normalized)")
ax.set_title("Immune-evasion ligand expression: CSC-high vs CSC-low cancer cells\n"
             "(red=upregulated in CSCs; grey=HLA-E/CD24 flagged)")
plt.tight_layout()
plt.savefig("results/figures/I4_checkpoint_axis.png", dpi=130, bbox_inches="tight")
plt.close()
print("\n  Saved: results/tables/I4_checkpoint_axis.csv")
print("  Saved: results/figures/I4_checkpoint_axis.png")
