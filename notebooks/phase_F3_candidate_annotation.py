# Phase F3 — Candidate Prioritization with Pathway / Regulon Context
#
# Integrates every line of evidence into one ranked, annotated CSC-candidate
# table for the wet lab and the paper:
#   - method ranks (Stage 1 DE, Stage 2 Geneformer attention)
#   - functional validation (F1 mean log2FC across ALDH/CD44-CD24/sphere gates)
#   - surface / druggability evidence (G7 GO cellular-component)
#   - CSC signaling-pathway membership (curated Wnt/Notch/Hh/BMP-TGFβ/JAK-STAT/Hippo)
#   - cross-cancer support (B2 GBM, B3 melanoma)
#   - a composite priority score
#
# Also runs pathway enrichment (Enrichr via gseapy, network-optional) on the
# Stage 1 vs Stage 2 top-gene lists to show the two methods load different
# pathway space.
#
# Output: results/tables/F3_candidate_priority.csv
#         results/tables/F3_pathway_enrichment.csv   (if Enrichr reachable)
#         results/figures/F3_candidate_landscape.png

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
os.makedirs("results/tables",  exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

print("=" * 64)
print("F3 — CANDIDATE PRIORITIZATION + PATHWAY CONTEXT")
print("=" * 64)

# ─────────────────────────────────────────────────────────────────
# LOAD ALL EVIDENCE LAYERS
# ─────────────────────────────────────────────────────────────────
s1   = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
s2   = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
fval = pd.read_csv("results/tables/F1_gene_level_validation.csv")
surf = pd.read_csv("results/tables/G7_surface_markers.csv")
gbm  = pd.read_csv("results/tables/B2_gbm_tumor_markers.csv")
mel  = pd.read_csv("results/tables/B3_melanoma_csc_markers.csv")

DE_RANK = {g: i+1 for i, g in enumerate(s1["gene_symbol"])}
GF_RANK = {g: i+1 for i, g in enumerate(s2["gene_symbol"])}
GF_SCORE = dict(zip(s2["gene_symbol"], s2["attention_score"]))
SURF_SET = dict(zip(surf["gene_symbol"], surf["surface_evidence"]))
GBM_SET = set(gbm["gene_symbol"].head(200))
MEL_SET = set(mel["gene_symbol"].head(200))

CONSENSUS = ['B2M','CD44','EEF1A1','EPCAM','FN1','HLA-E','HMGA1','ID3',
             'KLF4','MYC','RPL22L1','S100A6','SERBP1','SERPINE2','SOX9','VIM']

# ─────────────────────────────────────────────────────────────────
# CURATED CSC SIGNALING-PATHWAY MEMBERSHIP
# (canonical developmental/stemness pathways; hand-curated from KEGG/Reactome)
# ─────────────────────────────────────────────────────────────────
PATHWAYS = {
    "Wnt":        {"FZD7","FZD1","FZD2","WNT5A","WNT5B","LRP6","CTNNB1","TCF7",
                   "LEF1","AXIN2","DKK1","SFRP1","RSPO3","KLF4","MYC","SOX9","ID3"},
    "Notch":      {"NOTCH1","NOTCH2","NOTCH3","JAG1","JAG2","DLL1","DLL3","HES1",
                   "HEY1","HEYL","DTX1","NRARP"},
    "Hedgehog":   {"GLI1","GLI2","PTCH1","SMO","HHIP","SHH","BMI1"},
    "BMP/TGFβ":   {"BMPR1B","BMPR2","BMP2","BMP4","BMP7","ID1","ID2","ID3","SMAD1",
                   "SMAD4","TGFBR2","SERPINE2","SERPINE1","INHBA","GREM1"},
    "JAK/STAT":   {"OSMR","IL6R","IL6ST","STAT3","STAT1","JAK1","JAK2","SOCS3",
                   "LIFR","CNTFR","IL6","OSM"},
    "Hippo/YAP":  {"YAP1","WWTR1","TEAD1","TEAD4","AMOTL2","CTGF","CYR61","ANKRD1"},
    "EMT/mesench":{"VIM","FN1","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","CDH2",
                   "S100A4","SERPINE2","TIMP1"},
    "Pluripotency":{"SOX2","POU5F1","NANOG","KLF4","MYC","LIN28A","LIN28B",
                    "SALL4","SOX9","HMGA1","HMGA2"},
}

def pathways_of(gene):
    return [p for p, members in PATHWAYS.items() if gene in members]

# ─────────────────────────────────────────────────────────────────
# BUILD MASTER CANDIDATE TABLE
# ─────────────────────────────────────────────────────────────────
# Candidate universe = top-50 of each method ∪ consensus ∪ functionally-validated ∪ surface
cand = set(s1["gene_symbol"].head(50)) | set(s2["gene_symbol"].head(50))
cand |= set(CONSENSUS) | set(fval["gene"]) | set(surf["gene_symbol"].head(40))

fval_idx = fval.set_index("gene")
rows = []
for g in cand:
    de_r = DE_RANK.get(g)
    gf_r = GF_RANK.get(g)
    func = fval_idx.loc[g, "mean_func_log2FC"] if g in fval_idx.index else np.nan
    gates_up = int(fval_idx.loc[g, "n_gates_up"]) if g in fval_idx.index else 0
    gates_n  = int(fval_idx.loc[g, "n_gates_measured"]) if g in fval_idx.index else 0
    paths = pathways_of(g)
    rows.append({
        "gene": g,
        "DE_rank": de_r, "GF_rank": gf_r,
        "GF_attention": round(GF_SCORE.get(g, np.nan), 3) if g in GF_SCORE else np.nan,
        "func_log2FC": round(func, 3) if pd.notna(func) else np.nan,
        "func_gates_up": f"{gates_up}/{gates_n}" if gates_n else "",
        "surface": g in SURF_SET,
        "surface_evidence": SURF_SET.get(g, ""),
        "pathways": "; ".join(paths),
        "in_consensus": g in CONSENSUS,
        "in_GBM": g in GBM_SET,
        "in_melanoma": g in MEL_SET,
        "GF_unique": (gf_r is not None and gf_r <= 50) and (de_r is None or de_r > 200),
    })
cand_df = pd.DataFrame(rows)

# ── composite priority score (transparent, weighted) ──
def priority(r):
    s = 0.0
    if pd.notna(r["GF_rank"]): s += max(0, 50 - r["GF_rank"]) / 50 * 3      # attention prominence
    if pd.notna(r["DE_rank"]) and r["DE_rank"] <= 200: s += 1               # DE support
    if pd.notna(r["func_log2FC"]): s += np.clip(r["func_log2FC"], 0, 2.5)   # functional validation (key)
    if r["surface"]: s += 1.5                                               # druggability
    if r["pathways"]: s += 1                                                # mechanistic context
    s += 0.5 * int(r["in_GBM"]) + 0.5 * int(r["in_melanoma"])               # cross-cancer
    return round(s, 2)

cand_df["priority_score"] = cand_df.apply(priority, axis=1)
cand_df = cand_df.sort_values("priority_score", ascending=False).reset_index(drop=True)
cand_df.to_csv("results/tables/F3_candidate_priority.csv", index=False)

print(f"\nRanked {len(cand_df)} candidates. Top 20 by composite priority:")
print(f"  {'gene':10s} {'DE#':>5} {'GF#':>5} {'funcLFC':>8} {'surf':>5} {'path':>16} {'score':>6}")
for _, r in cand_df.head(20).iterrows():
    de = str(int(r["DE_rank"])) if pd.notna(r["DE_rank"]) else "-"
    gf = str(int(r["GF_rank"])) if pd.notna(r["GF_rank"]) else "-"
    fl = f"{r['func_log2FC']:+.2f}" if pd.notna(r["func_log2FC"]) else "-"
    sf = "Y" if r["surface"] else ""
    pa = (r["pathways"][:15]) if r["pathways"] else ""
    print(f"  {r['gene']:10s} {de:>5} {gf:>5} {fl:>8} {sf:>5} {pa:>16} {r['priority_score']:>6.2f}")

# ── Druggable surface + functional + Geneformer-unique = the wet-lab leads ──
print("\n── WET-LAB LEADS (surface + pathway + functional support) ──")
leads = cand_df[(cand_df["surface"]) & (cand_df["pathways"] != "") &
                (cand_df["func_log2FC"].fillna(-9) > 0)].head(8)
for _, r in leads.iterrows():
    novel = "Geneformer-unique" if r["GF_unique"] else "DE+GF"
    print(f"  {r['gene']:10s}  GF#{int(r['GF_rank']) if pd.notna(r['GF_rank']) else '-':>3}  "
          f"func={r['func_log2FC']:+.2f} ({r['func_gates_up']})  {r['pathways']:18s}  [{novel}]")
    print(f"             surface: {r['surface_evidence'][:60]}")

# ─────────────────────────────────────────────────────────────────
# PATHWAY MEMBERSHIP SUMMARY — DE vs Geneformer top 50
# ─────────────────────────────────────────────────────────────────
print("\n── PATHWAY LOADING: DE top-50 vs Geneformer top-50 ──")
de50 = set(s1["gene_symbol"].head(50))
gf50 = set(s2["gene_symbol"].head(50))
path_rows = []
for p, members in PATHWAYS.items():
    de_hits = sorted(de50 & members)
    gf_hits = sorted(gf50 & members)
    path_rows.append({"pathway": p, "DE_n": len(de_hits), "GF_n": len(gf_hits),
                      "DE_genes": ",".join(de_hits), "GF_genes": ",".join(gf_hits)})
    print(f"  {p:14s}  DE={len(de_hits)} [{','.join(de_hits)[:30]}]   GF={len(gf_hits)} [{','.join(gf_hits)[:30]}]")
path_df = pd.DataFrame(path_rows)

# ─────────────────────────────────────────────────────────────────
# ENRICHR (network-optional)
# ─────────────────────────────────────────────────────────────────
try:
    import gseapy as gp
    print("\nRunning Enrichr (MSigDB Hallmark) on each method's top-100...")
    enr_results = []
    for name, genes in [("Stage1_DE", s1["gene_symbol"].head(100).tolist()),
                        ("Stage2_Geneformer", s2["gene_symbol"].head(100).tolist())]:
        e = gp.enrichr(gene_list=genes, gene_sets="MSigDB_Hallmark_2020",
                       outdir=None, no_plot=True)
        top = e.results.sort_values("Adjusted P-value").head(8)
        top["method"] = name
        enr_results.append(top[["method","Term","Adjusted P-value","Genes"]])
        print(f"\n  {name} — top Hallmark pathways:")
        for _, rr in top.head(5).iterrows():
            print(f"    {rr['Term'][:45]:45s}  padj={rr['Adjusted P-value']:.2e}")
    enr_df = pd.concat(enr_results)
    enr_df.to_csv("results/tables/F3_pathway_enrichment.csv", index=False)
    print("\n  Saved: results/tables/F3_pathway_enrichment.csv")
    enrichr_ok = True
except Exception as e:
    print(f"\n  Enrichr unavailable ({type(e).__name__}); using curated pathway table only.")
    enrichr_ok = False

# ─────────────────────────────────────────────────────────────────
# FIGURE — candidate landscape
# ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(19, 6))

# Panel 1: top candidates by priority, colored by surface/novelty
top15 = cand_df.head(15)[::-1]
colors = ["#2ca02c" if s else "#1f77b4" for s in top15["surface"]]
axes[0].barh(top15["gene"], top15["priority_score"], color=colors)
axes[0].set_xlabel("Composite priority score")
axes[0].set_title("Top 15 CSC candidates\n(green=surface/druggable, blue=other)")
for i,(_,r) in enumerate(top15.iterrows()):
    if r["pathways"]:
        axes[0].text(r["priority_score"]+0.05, i, r["pathways"].split(";")[0][:10],
                     va="center", fontsize=6, color="#555")

# Panel 2: functional validation vs Geneformer rank (the key scatter)
sub = cand_df[cand_df["GF_rank"].notna() & cand_df["func_log2FC"].notna()]
sc_colors = ["#d62728" if u else "#1f77b4" for u in sub["GF_unique"]]
axes[1].scatter(sub["GF_rank"], sub["func_log2FC"], c=sc_colors, s=40, alpha=0.7)
for _, r in sub.iterrows():
    if r["func_log2FC"] > 0.5 or r["GF_rank"] <= 15:
        axes[1].annotate(r["gene"], (r["GF_rank"], r["func_log2FC"]), fontsize=7)
axes[1].axhline(0, color="k", lw=0.5)
axes[1].set_xlabel("Geneformer attention rank (lower=more important)")
axes[1].set_ylabel("Functional validation (mean log2FC)")
axes[1].set_title("Attention rank vs functional support\n(red=Geneformer-unique)")

# Panel 3: pathway loading DE vs GF
x = np.arange(len(path_df)); w = 0.38
axes[2].barh(x - w/2, path_df["DE_n"], w, label="Stage1 DE", color="#d62728")
axes[2].barh(x + w/2, path_df["GF_n"], w, label="Stage2 Geneformer", color="#1f77b4")
axes[2].set_yticks(x); axes[2].set_yticklabels(path_df["pathway"], fontsize=9)
axes[2].set_xlabel("Genes in top-50")
axes[2].set_title("CSC pathway loading\nby method")
axes[2].legend(fontsize=8)

plt.suptitle("Phase F3: Integrated CSC Candidate Prioritization & Pathway Context", fontsize=13)
plt.tight_layout()
plt.savefig("results/figures/F3_candidate_landscape.png", dpi=120, bbox_inches="tight")
plt.close()

print("\n" + "=" * 64)
print("PHASE F3 COMPLETE")
print("=" * 64)
print("  Saved: results/tables/F3_candidate_priority.csv")
print("  Saved: results/figures/F3_candidate_landscape.png")
