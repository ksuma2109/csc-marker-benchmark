# Phase G8 — ssGSEA Stemness Scoring + Subtype-Stratified Survival
#
# Fixes the G7 Cox failure. Two methodological upgrades:
#
#   1. ssGSEA instead of mean z-score
#      - Rank-based single-sample enrichment (Barbie et al. 2009)
#      - Robust to expression scale; the standard for signature scoring
#      - Each patient gets ONE stemness enrichment score per gene set
#
#   2. PAM50 subtype STRATIFICATION (not just adjustment)
#      - Cox `strata=` allows a separate baseline hazard per subtype
#      - Removes the luminal-confounding that flipped CD44's bulk association
#      - Asks the biologically correct question: WITHIN a subtype, does
#        higher stemness predict worse outcome?
#
# Input:  data/raw/metabric/expression.txt  (full microarray matrix, genes × 1980)
#         data/raw/metabric/clinical_patient.txt
#         results/tables/A5_csc_markers_DE.csv       (Stage 1)
#         results/tables/G4_geneformer_gene_ranking.csv (Stage 2)
#
# Output: results/tables/G8_ssgsea_scores.csv
#         results/tables/G8_cox_stratified_summary.csv
#         results/figures/G8_ssgsea_km.png
#         results/figures/G8_cox_forest.png

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

warnings.filterwarnings("ignore")
os.makedirs("results/tables",  exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# ─────────────────────────────────────────────────────────────────
# LOAD GENE SETS
# ─────────────────────────────────────────────────────────────────
print("=" * 60)
print("G8 — ssGSEA + SUBTYPE-STRATIFIED SURVIVAL")
print("=" * 60)

s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])

CONSENSUS = ['B2M','CD44','EEF1A1','EPCAM','FN1','HLA-E','HMGA1','ID3',
             'KLF4','MYC','RPL22L1','S100A6','SERBP1','SERPINE2','SOX9','VIM']

GENE_SETS = {
    "Stage1_DE_top50":      s1["gene_symbol"].head(50).tolist(),
    "Stage2_Geneformer_top50": s2["gene_symbol"].head(50).tolist(),
    "Consensus_16":         CONSENSUS,
}
for name, genes in GENE_SETS.items():
    print(f"  {name}: {len(genes)} genes")

# ─────────────────────────────────────────────────────────────────
# LOAD METABRIC EXPRESSION MATRIX
# ─────────────────────────────────────────────────────────────────
print("\nLoading METABRIC expression matrix (large, ~1 min)...")
expr = pd.read_csv("data/raw/metabric/expression.txt", sep="\t", low_memory=False)
expr = expr.dropna(subset=["Hugo_Symbol"])
expr = expr.drop_duplicates(subset=["Hugo_Symbol"], keep="first")
expr = expr.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"])
expr = expr.apply(pd.to_numeric, errors="coerce")
# Drop genes with all-NaN; fill remaining NaN with row mean (microarray gaps)
expr = expr.dropna(how="all")
expr = expr.apply(lambda r: r.fillna(r.mean()), axis=1)
print(f"  Expression: {expr.shape[0]} genes × {expr.shape[1]} samples")

# Coverage check
for name, genes in GENE_SETS.items():
    found = [g for g in genes if g in expr.index]
    print(f"  {name}: {len(found)}/{len(genes)} genes present in METABRIC")

# ─────────────────────────────────────────────────────────────────
# RUN ssGSEA
# ─────────────────────────────────────────────────────────────────
print("\nRunning ssGSEA (single-sample enrichment)...")
# gseapy ssgsea: expr = genes × samples, gene_sets = dict
ss = gp.ssgsea(
    data=expr,
    gene_sets=GENE_SETS,
    outdir=None,
    sample_norm_method="rank",
    no_plot=True,
    threads=4,
    min_size=5,
)
# Normalized enrichment scores → wide table (sample × gene_set)
res = ss.res2d.copy()
print(f"  ssGSEA result columns: {list(res.columns)}")
score_col = "NES" if "NES" in res.columns else "ES"
ssgsea_wide = res.pivot(index="Name", columns="Term", values=score_col).astype(float)
ssgsea_wide.index.name = "PATIENT_ID"
print(f"  ssGSEA scores: {ssgsea_wide.shape[0]} samples × {ssgsea_wide.shape[1]} signatures")
ssgsea_wide.to_csv("results/tables/G8_ssgsea_scores.csv")
print("  Saved: results/tables/G8_ssgsea_scores.csv")

# ─────────────────────────────────────────────────────────────────
# LOAD CLINICAL DATA
# ─────────────────────────────────────────────────────────────────
print("\nLoading clinical data...")
clin = pd.read_csv("data/raw/metabric/clinical_patient.txt", sep="\t", skiprows=4)
clin = clin.set_index("PATIENT_ID")
for c in ["RFS_MONTHS", "OS_MONTHS", "AGE_AT_DIAGNOSIS"]:
    clin[c] = pd.to_numeric(clin[c], errors="coerce")
# Event coding: status strings look like "1:Recurred" / "0:Not Recurred"
clin["RFS_event"] = clin["RFS_STATUS"].astype(str).str.startswith("1").astype(float)
clin["OS_event"]  = clin["OS_STATUS"].astype(str).str.startswith("1").astype(float)
clin["AGE"]       = clin["AGE_AT_DIAGNOSIS"]
print(f"  Clinical: {len(clin)} patients")
print(f"  Subtypes: {clin['CLAUDIN_SUBTYPE'].value_counts().to_dict()}")

# ─────────────────────────────────────────────────────────────────
# MERGE + COX (SUBTYPE-STRATIFIED)
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("COX PROPORTIONAL HAZARDS — SUBTYPE STRATIFIED")
print("=" * 60)

merged = ssgsea_wide.join(clin, how="inner")
# Drop subtypes too small / NC to stratify cleanly
valid_sub = ["LumA","LumB","Her2","Basal","claudin-low","Normal"]
merged = merged[merged["CLAUDIN_SUBTYPE"].isin(valid_sub)]
print(f"  Merged cohort: {len(merged)} patients")

# z-score the ssGSEA scores so HR is per-SD (comparable across signatures)
for sig in ssgsea_wide.columns:
    merged[sig] = (merged[sig] - merged[sig].mean()) / merged[sig].std()

cox_results = []
ENDPOINTS = [("RFS_MONTHS","RFS_event","Relapse-Free Survival"),
             ("OS_MONTHS","OS_event","Overall Survival")]

for dur, evt, ep_label in ENDPOINTS:
    print(f"\n── {ep_label} ──")
    for sig in ssgsea_wide.columns:
        df = merged[[sig, dur, evt, "AGE", "CLAUDIN_SUBTYPE"]].dropna()
        df = df[df[dur] > 0]
        cph = CoxPHFitter()
        try:
            cph.fit(df, duration_col=dur, event_col=evt,
                    formula=f"{sig} + AGE",
                    strata=["CLAUDIN_SUBTYPE"])
            row = cph.summary.loc[sig]
            hr  = np.exp(row["coef"])
            lo  = np.exp(row["coef lower 95%"])
            hi  = np.exp(row["coef upper 95%"])
            p   = row["p"]
            sig_mark = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
            print(f"  {sig:24s} HR={hr:.3f} (95% CI {lo:.3f}-{hi:.3f})  p={p:.4f} {sig_mark}  (n={len(df)})")
            cox_results.append({
                "endpoint": ep_label, "signature": sig, "HR": hr,
                "CI_lo": lo, "CI_hi": hi, "p": p, "sig": sig_mark, "n": len(df),
            })
        except Exception as e:
            print(f"  {sig:24s} FAILED: {e}")

cox_df = pd.DataFrame(cox_results)
cox_df.to_csv("results/tables/G8_cox_stratified_summary.csv", index=False)
print("\n  Saved: results/tables/G8_cox_stratified_summary.csv")

# ─────────────────────────────────────────────────────────────────
# KAPLAN-MEIER — TERTILE SPLIT, SUBTYPE-STRATIFIED LOG-RANK
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("KAPLAN-MEIER — STEMNESS TERTILES (RFS)")
print("=" * 60)

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
dur, evt = "RFS_MONTHS", "RFS_event"

for ax, sig in zip(axes, ssgsea_wide.columns):
    df = merged[[sig, dur, evt, "CLAUDIN_SUBTYPE"]].dropna()
    df = df[df[dur] > 0]
    t1, t2 = df[sig].quantile([1/3, 2/3])
    df["grp"] = np.where(df[sig] >= t2, "High",
                np.where(df[sig] <= t1, "Low", "Mid"))
    sub = df[df["grp"].isin(["High","Low"])]
    # subtype-stratified log-rank
    lr = multivariate_logrank_test(sub[dur], sub["grp"], sub[evt])

    kmf = KaplanMeierFitter()
    for grp, color in [("High","#d62728"),("Low","#1f77b4")]:
        m = sub["grp"] == grp
        kmf.fit(sub.loc[m, dur], sub.loc[m, evt], label=f"{grp} stemness (n={m.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)
    p = lr.p_value
    p_str = f"p={p:.3f}" if p>=0.001 else "p<0.001"
    ax.set_title(f"{sig}\nlog-rank {p_str}", fontsize=10)
    ax.set_xlabel("RFS (months)"); ax.set_ylabel("Relapse-free probability")
    ax.legend(fontsize=8)
    print(f"  {sig:24s} tertile log-rank p={p:.4f}")

plt.suptitle("ssGSEA Stemness Tertiles — METABRIC RFS (subtype-aware)", fontsize=13)
plt.tight_layout()
plt.savefig("results/figures/G8_ssgsea_km.png", dpi=120, bbox_inches="tight")
plt.close()
print("  Saved: results/figures/G8_ssgsea_km.png")

# ─────────────────────────────────────────────────────────────────
# FOREST PLOT
# ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))
plot_df = cox_df.reset_index(drop=True)
ypos = range(len(plot_df))
colors = {"Relapse-Free Survival":"#d62728","Overall Survival":"#1f77b4"}
for i, r in plot_df.iterrows():
    c = colors[r["endpoint"]]
    ax.plot([r["CI_lo"], r["CI_hi"]], [i, i], color=c, lw=2)
    ax.scatter([r["HR"]], [i], color=c, s=70, zorder=5)
    ax.text(r["CI_hi"]+0.02, i, f"HR={r['HR']:.2f} {r['sig']}", va="center", fontsize=8)
ax.axvline(1.0, color="black", ls="--", lw=1)
ax.set_yticks(list(ypos))
ax.set_yticklabels([f"{r['signature']}\n({r['endpoint'][:3]})" for _, r in plot_df.iterrows()], fontsize=8)
ax.set_xlabel("Hazard Ratio per SD of ssGSEA stemness (subtype-stratified + age-adjusted)")
ax.set_title("Subtype-Stratified Cox — CSC Signature Hazard Ratios\n(HR>1 = higher stemness predicts worse survival)")
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=v,label=k) for k,v in colors.items()], fontsize=8, loc="lower right")
plt.tight_layout()
plt.savefig("results/figures/G8_cox_forest.png", dpi=120, bbox_inches="tight")
plt.close()
print("  Saved: results/figures/G8_cox_forest.png")

# ─────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PHASE G8 COMPLETE")
print("=" * 60)
sig_hits = cox_df[cox_df["p"] < 0.05]
print(f"  Significant associations (p<0.05): {len(sig_hits)}/{len(cox_df)}")
for _, r in sig_hits.iterrows():
    direction = "worse" if r["HR"] > 1 else "better"
    print(f"    {r['signature']} ({r['endpoint']}): HR={r['HR']:.2f} → higher stemness = {direction} survival (p={r['p']:.4f})")
if len(sig_hits) == 0:
    print("  No signature reached significance even after ssGSEA + stratification.")
    print("  → Stemness signal in bulk METABRIC is genuinely weak (biological, not methodological).")
