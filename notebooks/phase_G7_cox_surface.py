# Phase G7 — Cox Survival Analysis + Surface Marker Prioritization
#
# Part A: Cox proportional hazards model on METABRIC
#   — CSC signature score as predictor, PAM50 subtype + age as covariates
#   — Removes subtype confounding that made log-rank tests non-significant
#
# Part B: Surface marker filter on Stage 2 top 200 genes
#   — Identifies transmembrane / secreted proteins (druggable targets)
#   — Ranks by Geneformer attention score
#
# Output: results/figures/G7_cox_*.png
#         results/tables/G7_cox_summary.csv
#         results/tables/G7_surface_markers.csv

import os, warnings, requests
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mygene

warnings.filterwarnings("ignore")
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables",  exist_ok=True)

# ─────────────────────────────────────────────────────────────────
# LOAD GENE LISTS
# ─────────────────────────────────────────────────────────────────
stage1_df    = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
stage1_genes = stage1_df["gene_symbol"].dropna().head(200).tolist()

stage2_df    = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
stage2_genes = stage2_df["gene_symbol"].head(200).tolist()

CONSENSUS = ['B2M','CD44','EEF1A1','EPCAM','FN1','HLA-E','HMGA1','ID3',
             'KLF4','MYC','RPL22L1','S100A6','SERBP1','SERPINE2','SOX9','VIM']

# ─────────────────────────────────────────────────────────────────
# PART A — COX SURVIVAL ANALYSIS (METABRIC)
# ─────────────────────────────────────────────────────────────────
print("=" * 60)
print("PART A — COX PROPORTIONAL HAZARDS (METABRIC)")
print("=" * 60)

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

CBIO  = "https://www.cbioportal.org/api"
STUDY = "brca_metabric"
MPROF = "brca_metabric_mrna_median_all_sample_Zscores"

# Fetch expression for all signature genes
all_sig = list(set(stage1_genes[:50] + stage2_genes[:50] + CONSENSUS))
mg_inst = mygene.MyGeneInfo()
hits = mg_inst.querymany(all_sig, scopes="symbol", fields="entrezgene", species="human")
sym2e = {r["query"]: int(r["entrezgene"]) for r in hits
         if "entrezgene" in r and r["entrezgene"]}
e2sym = {v: k for k, v in sym2e.items()}

r = requests.get(f"{CBIO}/studies/{STUDY}/samples",
                 params={"pageSize": 3000}, timeout=30)
sample_ids = [s["sampleId"] for s in r.json()]

CHUNK, all_expr = 400, []
for i in range(0, len(sym2e), CHUNK):
    chunk = list(sym2e.values())[i:i+CHUNK]
    r = requests.post(f"{CBIO}/molecular-profiles/{MPROF}/molecular-data/fetch",
        json={"sampleIds": sample_ids, "entrezGeneIds": chunk},
        params={"projection": "SUMMARY"}, timeout=120)
    if r.ok: all_expr.extend(r.json())

ew = pd.DataFrame(all_expr).pivot(
    index="sampleId", columns="entrezGeneId", values="value")
ew.columns = [e2sym.get(c, str(c)) for c in ew.columns]
ew = ew.apply(pd.to_numeric, errors="coerce")
ew.index.name = "patientId"
print(f"Expression matrix: {ew.shape}")

# Clinical data
r = requests.get(f"{CBIO}/studies/{STUDY}/clinical-data",
    params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"}, timeout=60)
clin_all = pd.DataFrame(r.json())
clin = (clin_all[clin_all["clinicalAttributeId"].isin(
            ["RFS_MONTHS","RFS_STATUS","OS_MONTHS","OS_STATUS",
             "CLAUDIN_SUBTYPE","AGE_AT_DIAGNOSIS"])]
        .pivot(index="patientId", columns="clinicalAttributeId", values="value"))
clin["RFS_MONTHS"] = pd.to_numeric(clin["RFS_MONTHS"], errors="coerce")
clin["OS_MONTHS"]  = pd.to_numeric(clin["OS_MONTHS"],  errors="coerce")
clin["AGE"]        = pd.to_numeric(clin["AGE_AT_DIAGNOSIS"], errors="coerce")
clin["RFS_event"]  = clin["RFS_STATUS"].str.startswith("1").astype(float)
clin["OS_event"]   = clin["OS_STATUS"].str.startswith("1").astype(float)

merged = ew.join(clin, how="inner").dropna(subset=["RFS_MONTHS", "RFS_event"])
print(f"Patients with RFS: {len(merged)}")

# ── Score each signature (mean z-score of available genes) ───────
def make_score(df, genes):
    avail = [g for g in genes if g in df.columns]
    return df[avail].mean(axis=1), avail

merged["score_s1"], s1_avail = make_score(merged, stage1_genes[:30])
merged["score_s2"], s2_avail = make_score(merged, stage2_genes[:30])
merged["score_con"], c_avail = make_score(merged, CONSENSUS)

print(f"\nGenes used — S1:{len(s1_avail)}  S2:{len(s2_avail)}  Consensus:{len(c_avail)}")

# One-hot encode subtype (LumA as reference — best prognosis baseline)
subtype_order = ["LumA","LumB","Her2","Basal","claudin-low","Normal"]
subtype_safe  = {"LumA":"LumA","LumB":"LumB","Her2":"Her2",
                 "Basal":"Basal","claudin-low":"CaudLow","Normal":"Normal"}
for st in subtype_order[1:]:     # drop LumA (reference)
    merged[f"sub_{subtype_safe[st]}"] = (merged["CLAUDIN_SUBTYPE"] == st).astype(float)
sub_cols = [f"sub_{subtype_safe[s]}" for s in subtype_order[1:]]

# ── Cox models ────────────────────────────────────────────────────
cox_results = []

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("Cox Proportional Hazards — METABRIC RFS\n"
             "(subtype + age corrected)", fontsize=13)

for ax, (score_col, label, color) in zip(axes, [
    ("score_s1",  "Stage 1 (DE top 30)",        "#d62728"),
    ("score_s2",  "Stage 2 (Geneformer top 30)", "#1f77b4"),
    ("score_con", "Consensus (16 genes)",         "#2ca02c"),
]):
    cox_cols = [score_col, "RFS_MONTHS", "RFS_event", "AGE"] + sub_cols
    df_cox   = merged[cox_cols].dropna()

    cph = CoxPHFitter()
    cph.fit(df_cox, duration_col="RFS_MONTHS", event_col="RFS_event",
            formula=f"{score_col} + AGE + " + " + ".join(sub_cols))

    row = cph.summary.loc[score_col]
    hr, lo, hi = np.exp(row["coef"]), np.exp(row["coef lower 95%"]), np.exp(row["coef upper 95%"])
    p           = row["p"]
    sig         = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

    print(f"\n{label}:")
    print(f"  HR = {hr:.3f}  (95% CI {lo:.3f}–{hi:.3f})  p = {p:.4e}  {sig}")
    cox_results.append({"signature": label, "HR": hr, "CI_lo": lo,
                        "CI_hi": hi, "p": p, "sig": sig, "n": len(df_cox)})

    # Stratified KM: quartile split WITHIN each subtype, then pool
    q25 = merged[score_col].quantile(0.25)
    q75 = merged[score_col].quantile(0.75)
    hi_mask = merged[score_col] >= q75
    lo_mask = merged[score_col] <= q25
    sub_hi  = merged.loc[hi_mask]
    sub_lo  = merged.loc[lo_mask]
    lr      = logrank_test(sub_hi["RFS_MONTHS"], sub_lo["RFS_MONTHS"],
                           event_observed_A=sub_hi["RFS_event"],
                           event_observed_B=sub_lo["RFS_event"])

    kmf = KaplanMeierFitter()
    kmf.fit(sub_hi["RFS_MONTHS"], sub_hi["RFS_event"],
            label=f"High Q4 (n={hi_mask.sum()})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color=color)
    kmf.fit(sub_lo["RFS_MONTHS"], sub_lo["RFS_event"],
            label=f"Low Q1 (n={lo_mask.sum()})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="gray")

    p_str = f"p={lr.p_value:.2e}" if lr.p_value >= 0.001 else f"p<0.001"
    ax.set_title(f"{label}\nCox HR={hr:.2f} ({sig})  |  KM {p_str}", fontsize=10)
    ax.set_xlabel("RFS (months)")
    ax.set_ylabel("Relapse-free probability")
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("results/figures/G7_cox_km.png", dpi=120, bbox_inches="tight")
plt.close()

# Forest plot of HRs
fig, ax = plt.subplots(figsize=(8, 4))
labels = [r["signature"] for r in cox_results]
hrs    = [r["HR"]    for r in cox_results]
lo95   = [r["CI_lo"] for r in cox_results]
hi95   = [r["CI_hi"] for r in cox_results]
colors_fp = ["#d62728","#1f77b4","#2ca02c"]
y = range(len(labels))
for i, (hr, lo, hi, sig, col) in enumerate(zip(hrs, lo95, hi95,
        [r["sig"] for r in cox_results], colors_fp)):
    ax.plot([lo, hi], [i, i], color=col, linewidth=2)
    ax.scatter([hr], [i], color=col, zorder=5, s=80)
    ax.text(hi + 0.01, i, f"HR={hr:.2f} {sig}", va="center", fontsize=9)
ax.axvline(x=1.0, color="black", linestyle="--", linewidth=1)
ax.set_yticks(list(y))
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel("Hazard Ratio (adjusted for PAM50 subtype + age)")
ax.set_title("Cox Model — CSC Signature Hazard Ratios\n(HR > 1 = worse relapse-free survival)")
plt.tight_layout()
plt.savefig("results/figures/G7_cox_forest.png", dpi=120, bbox_inches="tight")
plt.close()

cox_df = pd.DataFrame(cox_results)
cox_df.to_csv("results/tables/G7_cox_summary.csv", index=False)
print("\nSaved: results/figures/G7_cox_km.png")
print("Saved: results/figures/G7_cox_forest.png")
print("Saved: results/tables/G7_cox_summary.csv")

# ─────────────────────────────────────────────────────────────────
# PART B — SURFACE MARKER FILTER
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PART B — SURFACE MARKER PRIORITIZATION")
print("=" * 60)

# Query mygene for GO cellular component terms
print("\nQuerying gene annotations (GO cellular component)...")
hits = mg_inst.querymany(
    stage2_genes,
    scopes="symbol",
    fields="go.CC,uniprot.Swiss-Prot,name",
    species="human",
    returnall=False,
)

# GO terms indicating cell surface / membrane / secreted
SURFACE_GO = {
    "GO:0005886": "plasma membrane",
    "GO:0009986": "cell surface",
    "GO:0005576": "extracellular space",
    "GO:0005615": "extracellular space (soluble)",
    "GO:0031012": "extracellular matrix",
    "GO:0016021": "integral component of membrane",
    "GO:0005887": "integral to plasma membrane",
    "GO:0098552": "side of membrane",
    "GO:0098589": "membrane region",
    "GO:0045121": "membrane raft",
    "GO:0070062": "extracellular exosome",
}

surface_rows = []
for h in hits:
    sym = h.get("query")
    if sym not in stage2_genes: continue
    cc_terms = h.get("go", {}).get("CC", [])
    if isinstance(cc_terms, dict): cc_terms = [cc_terms]

    matched = []
    for term in cc_terms:
        goid = term.get("id", "")
        if goid in SURFACE_GO:
            matched.append(SURFACE_GO[goid])

    if matched:
        attn_row = stage2_df[stage2_df["gene_symbol"] == sym]
        attn_score = attn_row["attention_score"].values[0] if len(attn_row) else 0
        rank = stage2_genes.index(sym) + 1 if sym in stage2_genes else 999
        surface_rows.append({
            "gene_symbol":    sym,
            "attention_rank": rank,
            "attention_score": attn_score,
            "surface_evidence": ", ".join(sorted(set(matched))),
            "gene_name":       h.get("name", ""),
            "uniprot":         h.get("uniprot", {}).get("Swiss-Prot", ""),
            "in_stage1_top200": sym in stage1_genes,
            "in_consensus":     sym in CONSENSUS,
        })

surf_df = pd.DataFrame(surface_rows).sort_values("attention_rank").reset_index(drop=True)
surf_df.to_csv("results/tables/G7_surface_markers.csv", index=False)

print(f"\nStage 2 genes with surface/membrane evidence: {len(surf_df)}/200")
print(f"\n{'Rank':>5} {'Gene':>10} {'Attn':>7} {'In S1':>6} {'Consensus':>9}  Location")
print("-" * 80)
for _, row in surf_df.head(25).iterrows():
    loc = row["surface_evidence"][:40]
    s1  = "✓" if row["in_stage1_top200"] else " "
    con = "✓" if row["in_consensus"]     else " "
    print(f"{row['attention_rank']:>5} {row['gene_symbol']:>10} "
          f"{row['attention_score']:>7.3f} {s1:>6} {con:>9}  {loc}")

# Highlight top targets
print("\n── TOP THERAPEUTIC CANDIDATES (surface + high attention) ──")
top_targets = surf_df[surf_df["attention_rank"] <= 50].head(10)
for _, row in top_targets.iterrows():
    tag = ""
    if row["in_consensus"]:    tag += "[CONSENSUS] "
    if row["in_stage1_top200"]: tag += "[DE confirmed] "
    print(f"  #{row['attention_rank']:>3} {row['gene_symbol']:>10}  "
          f"attn={row['attention_score']:.3f}  {tag}")
    print(f"        {row['gene_name']}")
    print(f"        Location: {row['surface_evidence']}")

print(f"\nSaved: results/tables/G7_surface_markers.csv")
print("\n✓ Phase G7 complete.")
