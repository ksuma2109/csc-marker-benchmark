# Phase I5 — Subtype-stratified robustness of the CSC immune findings
#
# The critical confound for the immune story: TNBC is BOTH the most stem-like
# and the most immunogenic breast subtype. So the cell-level findings (CSC-high
# retain MHC-I; upregulate CD47/HLA-E) could be a cross-subtype artifact rather
# than a CSC effect. This tests whether the associations hold WITHIN each
# subtype — the make-or-break robustness check for the second manuscript.
#
# For each subtype (ER+, TNBC, HER2+): Spearman(stemness, score) and a
# within-subtype CSC-high vs CSC-low contrast (quartiles defined per subtype).
#
# Output: results/tables/I5_subtype_stratified.csv
#         results/figures/I5_subtype_stratified.png

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

# programs / single genes to test (name -> (genes, single?))
PROGRAMS = {
    "MHC-I antigen presentation": (["HLA-A","HLA-B","HLA-C","TAP1","TAP2","NLRC5"], False),
    "CD47 ('don't eat me')":      (["CD47"], True),
    "HLA-E (NKG2A)":              (["HLA-E"], True),
    "CD274 (PD-L1)":              (["CD274"], True),
}

print("=" * 66)
print("I5 — SUBTYPE-STRATIFIED ROBUSTNESS OF CSC IMMUNE FINDINGS")
print("=" * 66)

epi = sc.read_h5ad("data/processed/brca_A5_csc_scored.h5ad")

def get_score(name, genes, single):
    if single:
        g = genes[0]
        x = epi[:, g].X
        return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    sc.tl.score_genes(epi, [x for x in genes if x in epi.var_names], score_name="_s", ctrl_size=50)
    return epi.obs["_s"].values

subtypes = [s for s in ["ER+","TNBC","HER2+"] if s in set(epi.obs["subtype"])]
print(f"Subtypes: {subtypes}  |  cells: "
      + ", ".join(f"{s}={(epi.obs['subtype']==s).sum()}" for s in subtypes))

rows = []
for name, (genes, single) in PROGRAMS.items():
    score = get_score(name, genes, single)
    # overall (all cancer-epithelial)
    r_all, p_all = spearmanr(epi.obs["stemness_composite"], score)
    rows.append({"program": name, "subtype": "ALL", "n": epi.n_obs,
                 "spearman_r": round(r_all,3), "spearman_p": p_all})
    for st in subtypes:
        m = (epi.obs["subtype"] == st).values
        stem_st = epi.obs["stemness_composite"].values[m]
        sc_st = score[m]
        r, p = spearmanr(stem_st, sc_st)
        # within-subtype quartile contrast
        q1, q3 = np.quantile(stem_st, [0.25, 0.75])
        hi = sc_st[stem_st >= q3]; lo = sc_st[stem_st <= q1]
        try:
            _, pmw = mannwhitneyu(hi, lo, alternative="two-sided")
        except Exception:
            pmw = np.nan
        rows.append({"program": name, "subtype": st, "n": int(m.sum()),
                     "spearman_r": round(r,3), "spearman_p": p,
                     "csc_high_minus_low": round(float(hi.mean()-lo.mean()),3), "mw_p": pmw})

res = pd.DataFrame(rows)
res.to_csv("results/tables/I5_subtype_stratified.csv", index=False)

# print
for name in PROGRAMS:
    sub = res[res["program"]==name]
    print(f"\n── {name} ──")
    for _, r in sub.iterrows():
        extra = f"  CSChi-lo={r['csc_high_minus_low']:+.3f}" if r["subtype"]!="ALL" and pd.notna(r.get("csc_high_minus_low")) else ""
        sig = "***" if r["spearman_p"]<0.001 else "**" if r["spearman_p"]<0.01 else "*" if r["spearman_p"]<0.05 else "ns"
        print(f"    {r['subtype']:5s} (n={r['n']:5d})  Spearman r={r['spearman_r']:+.3f} {sig}{extra}")

# robustness verdict
print("\n" + "=" * 66)
print("ROBUSTNESS VERDICT (does the CSC association survive within-subtype?)")
print("=" * 66)
for name in PROGRAMS:
    sub = res[(res["program"]==name) & (res["subtype"]!="ALL")]
    same_dir = (np.sign(sub["spearman_r"]) == np.sign(res[(res["program"]==name)&(res["subtype"]=="ALL")]["spearman_r"].iloc[0])).sum()
    sig_within = (sub["spearman_p"] < 0.05).sum()
    print(f"  {name:30s}: consistent direction in {same_dir}/{len(sub)} subtypes; "
          f"significant in {sig_within}/{len(sub)}")

# ── figure: per-subtype Spearman r for each program ──
fig, ax = plt.subplots(figsize=(10, 5))
progs = list(PROGRAMS.keys())
x = np.arange(len(progs)); w = 0.2
cols = {"ALL":"#333","ER+":"#1f77b4","TNBC":"#d62728","HER2+":"#2ca02c"}
groups = ["ALL"] + subtypes
for i, st in enumerate(groups):
    vals = [res[(res["program"]==p)&(res["subtype"]==st)]["spearman_r"].values for p in progs]
    vals = [v[0] if len(v) else np.nan for v in vals]
    ax.bar(x + (i-(len(groups)-1)/2)*w, vals, w, label=st, color=cols.get(st,"#999"))
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(progs, rotation=15, ha="right", fontsize=8)
ax.set_ylabel("Spearman r (score vs stemness)")
ax.set_title("Subtype-stratified CSC–immune associations\n(consistent within-subtype direction = robust to subtype confound)")
ax.legend(fontsize=8, title="subtype")
plt.tight_layout()
plt.savefig("results/figures/I5_subtype_stratified.png", dpi=130, bbox_inches="tight")
plt.close()
print("\n  Saved: results/tables/I5_subtype_stratified.csv")
print("  Saved: results/figures/I5_subtype_stratified.png")
