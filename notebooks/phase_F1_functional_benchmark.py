# Phase F1 — Functional Ground-Truth Benchmark: DE vs Geneformer
#
# Breaks the circularity. Instead of validating CSC gene lists against the
# same published markers that defined the labels, we score them against
# INDEPENDENT FUNCTIONAL CSC assays (sorted populations / self-renewal):
#
#   Gate 1  ALDH+ (ALDEFLUOR, enzymatic CSC activity)   GSE115302 (PDX)
#   Gate 2  CD44+CD24- (Al-Hajj surface gate)           GSE115302 (PDX) + GSE36643
#   Gate 3  Mammosphere vs adherent (self-renewal)      GSE182532
#
# For each gate we compute a per-gene FUNCTIONAL log2FC (CSC vs non-CSC),
# then ask which method's gene ranking is better supported by function:
#
#   - precision@k   : fraction of method's top-k genes that are UP in CSC
#   - mean func LFC : effect size of method's top-k genes in the assay
#   - AUROC         : does the method's score rank functionally-up genes high?
#   - Spearman rho  : does the method's ranking track functional log2FC?
#   - vs random     : 1000 random gene sets as a null baseline
#
# Methods compared:
#   Stage 1 = DE Wilcoxon ranking      (results/tables/A5_csc_markers_DE.csv)
#   Stage 2 = Geneformer attention     (results/tables/G4_geneformer_gene_ranking.csv)
#
# Output: results/tables/F1_functional_benchmark.csv
#         results/tables/F1_gene_level_validation.csv
#         results/figures/F1_benchmark.png

import os, gzip, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")

os.makedirs("results/tables",  exist_ok=True)
os.makedirs("results/figures", exist_ok=True)
rng = np.random.default_rng(42)

# ─────────────────────────────────────────────────────────────────
# LOAD METHOD RANKINGS
# ─────────────────────────────────────────────────────────────────
print("=" * 64)
print("F1 — FUNCTIONAL GROUND-TRUTH BENCHMARK")
print("=" * 64)

s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])

# Method score per gene (continuous, higher = stronger CSC marker by that method)
DE_SCORE   = dict(zip(s1["gene_symbol"], s1["wilcoxon_score"]))
GF_SCORE   = dict(zip(s2["gene_symbol"], s2["attention_score"]))
DE_RANKED  = s1["gene_symbol"].tolist()        # ordered best→worst
GF_RANKED  = s2["gene_symbol"].tolist()
print(f"  Stage 1 DE: {len(DE_RANKED)} ranked genes")
print(f"  Stage 2 Geneformer: {len(GF_RANKED)} ranked genes")

# ─────────────────────────────────────────────────────────────────
# BUILD FUNCTIONAL CONTRASTS  (gene -> functional log2FC, CSC vs non-CSC)
# ─────────────────────────────────────────────────────────────────
print("\nBuilding functional contrasts...")
contrasts = {}   # name -> pd.Series(index=gene, value=log2FC)

# --- GSE115302 (PDX RNA-seq): ALDH & CD44 vs bulk, 2 models ---
pdx = pd.read_excel("data/raw/GSE115302/matrix.xlsx")
pdx["gene"] = pdx["gene"].astype(str)   # undo Excel date-corruption of gene names
pdx = pdx.dropna(subset=["gene"]).drop_duplicates(subset=["gene"], keep="first").set_index("gene")
for col in ["MC1.bulk","MC1.aldh","MC1.cd44","VAR068.bulk","VAR068.aldh","VAR068.cd44"]:
    pdx[col] = pd.to_numeric(pdx[col], errors="coerce")

def log2fc(num, den):
    return np.log2((pdx[num] + 1) / (pdx[den] + 1))

# Average the two PDX models per gate
contrasts["ALDH+ (PDX, GSE115302)"]      = (log2fc("MC1.aldh","MC1.bulk") + log2fc("VAR068.aldh","VAR068.bulk")) / 2
contrasts["CD44+CD24- (PDX, GSE115302)"] = (log2fc("MC1.cd44","MC1.bulk") + log2fc("VAR068.cd44","VAR068.bulk")) / 2

# --- GSE182532 (mammosphere vs adherent, DESeq2) ---
with gzip.open("data/raw/GSE182532/sphere_vs_adherent.txt.gz","rt") as f:
    sph = pd.read_csv(f)
# single comma-joined column -> split
if sph.shape[1] == 1:
    sph = sph.iloc[:,0].str.split(",", expand=True)
    sph.columns = ["Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
sph["log2FoldChange"] = pd.to_numeric(sph["log2FoldChange"], errors="coerce")
sph = sph.dropna(subset=["Gene","log2FoldChange"]).drop_duplicates("Gene")
contrasts["Mammosphere (GSE182532)"] = sph.set_index("Gene")["log2FoldChange"]

# --- GSE36643 (CD44hi/CD24lo ready DE list; FC is linear) ---
with gzip.open("data/raw/GSE36643/de_44hi24lo.csv.gz","rt") as f:
    g366 = pd.read_csv(f)
g366 = g366.dropna(subset=["gene symbol"]).drop_duplicates("gene symbol")
g366["gene symbol"] = g366["gene symbol"].astype(str)
# FC given as linear fold-change CD44hi/CD24lo vs CD44lo/CD24hi -> log2
contrasts["CD44hi/CD24lo (GSE36643)"] = pd.Series(
    np.log2(pd.to_numeric(g366["FC"], errors="coerce").clip(lower=1e-6).values),
    index=g366["gene symbol"].values)

# Force all contrast indices to clean strings (align safely across assays)
for name in list(contrasts.keys()):
    s = contrasts[name]
    s.index = s.index.astype(str)
    contrasts[name] = s[~s.index.duplicated(keep="first")]
    print(f"  {name}: {contrasts[name].notna().sum()} genes")

# ─────────────────────────────────────────────────────────────────
# BENCHMARK METRICS
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 64)
print("SCORING METHODS AGAINST FUNCTIONAL GATES")
print("=" * 64)

TOPK     = 100      # markers per method to evaluate
UP_THRESH = 0.0     # functional log2FC > 0 => up in CSC
N_RANDOM = 1000

def evaluate(method_name, ranked, score_map, func_lfc, gate):
    func = func_lfc.dropna()
    universe = list(func.index)            # genes measurable in this assay
    uni_set  = set(universe)

    # top-k present in this assay
    topk = [g for g in ranked if g in uni_set][:TOPK]
    if len(topk) < 10:
        return None
    lfc_topk = func.loc[topk]
    prec     = float((lfc_topk > UP_THRESH).mean())         # precision@k (directional)
    mean_lfc = float(lfc_topk.mean())

    # AUROC: method continuous score predicts "functionally up" over the universe
    truth = (func > func.quantile(0.75)).astype(int)        # top-quartile up = functional CSC gene
    mscore = np.array([score_map.get(g, 0.0) for g in universe])
    auroc = roc_auc_score(truth.values, mscore) if truth.nunique() == 2 else np.nan

    # Spearman: method score vs functional log2FC (over method genes present)
    shared = [g for g in ranked if g in uni_set]
    rho, _ = spearmanr([score_map.get(g,0.0) for g in shared], func.loc[shared].values)

    # Random-set null for precision@k
    rand_prec = []
    uni_arr = np.array(universe)
    for _ in range(N_RANDOM):
        samp = rng.choice(uni_arr, size=min(TOPK, len(uni_arr)), replace=False)
        rand_prec.append((func.loc[samp] > UP_THRESH).mean())
    rand_mean = float(np.mean(rand_prec))
    # empirical p: fraction of random sets with precision >= observed
    p_emp = float((np.array(rand_prec) >= prec).mean())

    return {
        "gate": gate, "method": method_name,
        f"precision@{TOPK}": round(prec, 3),
        "random_precision": round(rand_mean, 3),
        "enrichment": round(prec / rand_mean, 2) if rand_mean > 0 else np.nan,
        "p_vs_random": p_emp,
        "mean_func_log2FC": round(mean_lfc, 3),
        "AUROC": round(auroc, 3) if not np.isnan(auroc) else np.nan,
        "spearman_rho": round(rho, 3),
        "n_topk_in_assay": len(topk),
    }

rows = []
for gate, func in contrasts.items():
    for mname, ranked, smap in [("Stage1_DE", DE_RANKED, DE_SCORE),
                                 ("Stage2_Geneformer", GF_RANKED, GF_SCORE)]:
        r = evaluate(mname, ranked, smap, func, gate)
        if r: rows.append(r)

bench = pd.DataFrame(rows)
bench.to_csv("results/tables/F1_functional_benchmark.csv", index=False)

# Print per-gate comparison
for gate in contrasts:
    sub = bench[bench["gate"] == gate]
    if sub.empty: continue
    print(f"\n── {gate} ──")
    for _, r in sub.iterrows():
        print(f"  {r['method']:18s}  precision@{TOPK}={r[f'precision@{TOPK}']:.2f} "
              f"(rand {r['random_precision']:.2f}, {r['enrichment']}x, p={r['p_vs_random']:.3f})  "
              f"meanLFC={r['mean_func_log2FC']:+.2f}  AUROC={r['AUROC']}  rho={r['spearman_rho']}")

# ─────────────────────────────────────────────────────────────────
# HEAD-TO-HEAD SUMMARY  (who wins each gate, each metric)
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 64)
print("HEAD-TO-HEAD: Stage1 DE  vs  Stage2 Geneformer")
print("=" * 64)
metrics = [f"precision@{TOPK}", "mean_func_log2FC", "AUROC", "spearman_rho"]
wins = {"Stage1_DE": 0, "Stage2_Geneformer": 0}
for gate in contrasts:
    sub = bench[bench["gate"] == gate].set_index("method")
    if set(["Stage1_DE","Stage2_Geneformer"]) - set(sub.index): continue
    print(f"\n  {gate}")
    for m in metrics:
        de_v, gf_v = sub.loc["Stage1_DE", m], sub.loc["Stage2_Geneformer", m]
        winner = "DE" if de_v > gf_v else "Geneformer" if gf_v > de_v else "tie"
        if winner == "DE": wins["Stage1_DE"] += 1
        elif winner == "Geneformer": wins["Stage2_Geneformer"] += 1
        print(f"    {m:18s} DE={de_v}   GF={gf_v}   → {winner}")
print(f"\n  Metric wins:  DE={wins['Stage1_DE']}   Geneformer={wins['Stage2_Geneformer']}")

# ─────────────────────────────────────────────────────────────────
# GENE-LEVEL VALIDATION TABLE  (which specific candidates hold up)
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 64)
print("GENE-LEVEL VALIDATION — top candidates across functional gates")
print("=" * 64)

# Mean functional log2FC across all gates for each gene
all_func = pd.DataFrame(contrasts)
mean_func = all_func.mean(axis=1)
n_gates_up = (all_func > 0).sum(axis=1)

CONSENSUS = ['B2M','CD44','EEF1A1','EPCAM','FN1','HLA-E','HMGA1','ID3',
             'KLF4','MYC','RPL22L1','S100A6','SERBP1','SERPINE2','SOX9','VIM']
WATCH = ['SOX9','KLF4','FZD7','CD44','ALDH1A3','MYC','VIM','FN1','EPCAM',
         'OSMR','SERPINE2','ALDH1A1','PROM1','NES','ITGA6']

val_rows = []
for g in sorted(set(WATCH) | set(CONSENSUS) | set(DE_RANKED[:50]) | set(GF_RANKED[:50])):
    if g not in all_func.index: continue
    val_rows.append({
        "gene": g,
        "DE_rank": DE_RANKED.index(g)+1 if g in DE_RANKED else None,
        "GF_rank": GF_RANKED.index(g)+1 if g in GF_RANKED else None,
        "in_consensus": g in CONSENSUS,
        "mean_func_log2FC": round(float(mean_func.get(g, np.nan)), 3),
        "n_gates_up": int(n_gates_up.get(g, 0)),
        "n_gates_measured": int(all_func.loc[g].notna().sum()),
    })
val = pd.DataFrame(val_rows).sort_values("mean_func_log2FC", ascending=False)
val.to_csv("results/tables/F1_gene_level_validation.csv", index=False)

print("\n  Top 20 functionally-validated CSC genes (by mean functional log2FC):")
print(f"  {'gene':12s} {'DE#':>5} {'GF#':>5} {'meanLFC':>8} {'gatesUp':>8}  consensus")
for _, r in val.head(20).iterrows():
    de = str(int(r["DE_rank"])) if pd.notna(r["DE_rank"]) else "-"
    gf = str(int(r["GF_rank"])) if pd.notna(r["GF_rank"]) else "-"
    con = "★" if r["in_consensus"] else ""
    print(f"  {r['gene']:12s} {de:>5} {gf:>5} {r['mean_func_log2FC']:>+8.2f} "
          f"{r['n_gates_up']}/{r['n_gates_measured']:>1}      {con}")

# Geneformer-unique hits that validate functionally (DE missed them)
gf_unique = val[(val["GF_rank"].notna()) & (val["GF_rank"] <= 50) &
                ((val["DE_rank"].isna()) | (val["DE_rank"] > 200)) &
                (val["mean_func_log2FC"] > 0)]
print(f"\n  Geneformer-unique genes (top-50 GF, not DE top-200) that are UP functionally:")
for _, r in gf_unique.iterrows():
    print(f"    {r['gene']:12s} GF#{int(r['GF_rank'])}  meanLFC={r['mean_func_log2FC']:+.2f}  "
          f"up in {r['n_gates_up']}/{r['n_gates_measured']} gates")

# ─────────────────────────────────────────────────────────────────
# FIGURE
# ─────────────────────────────────────────────────────────────────
gates = list(contrasts.keys())
fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

# Panel 1: precision@k DE vs GF per gate, with random baseline
x = np.arange(len(gates)); w = 0.35
de_prec = [bench[(bench.gate==g)&(bench.method=="Stage1_DE")][f"precision@{TOPK}"].values for g in gates]
gf_prec = [bench[(bench.gate==g)&(bench.method=="Stage2_Geneformer")][f"precision@{TOPK}"].values for g in gates]
de_prec = [v[0] if len(v) else 0 for v in de_prec]
gf_prec = [v[0] if len(v) else 0 for v in gf_prec]
rand_b  = [bench[bench.gate==g]["random_precision"].mean() for g in gates]
axes[0].bar(x-w/2, de_prec, w, label="Stage1 DE", color="#d62728")
axes[0].bar(x+w/2, gf_prec, w, label="Stage2 Geneformer", color="#1f77b4")
axes[0].plot(x, rand_b, "k--o", label="Random baseline", markersize=5)
axes[0].set_xticks(x); axes[0].set_xticklabels([g.split(" (")[0] for g in gates], rotation=25, ha="right", fontsize=8)
axes[0].set_ylabel(f"precision@{TOPK} (fraction up in CSC)")
axes[0].set_title("Functional precision of top markers\nvs. random baseline")
axes[0].legend(fontsize=8); axes[0].set_ylim(0,1)

# Panel 2: AUROC per gate
de_auc = [bench[(bench.gate==g)&(bench.method=="Stage1_DE")]["AUROC"].values for g in gates]
gf_auc = [bench[(bench.gate==g)&(bench.method=="Stage2_Geneformer")]["AUROC"].values for g in gates]
de_auc = [v[0] if len(v) else np.nan for v in de_auc]
gf_auc = [v[0] if len(v) else np.nan for v in gf_auc]
axes[1].bar(x-w/2, de_auc, w, label="Stage1 DE", color="#d62728")
axes[1].bar(x+w/2, gf_auc, w, label="Stage2 Geneformer", color="#1f77b4")
axes[1].axhline(0.5, color="k", ls="--", label="chance")
axes[1].set_xticks(x); axes[1].set_xticklabels([g.split(" (")[0] for g in gates], rotation=25, ha="right", fontsize=8)
axes[1].set_ylabel("AUROC (rank functional-CSC genes)")
axes[1].set_title("Discrimination of functional CSC genes")
axes[1].legend(fontsize=8); axes[1].set_ylim(0,1)

# Panel 3: gene-level — mean functional LFC of top validated genes
top_val = val.head(15)[::-1]
colors = ["#2ca02c" if c else ("#1f77b4" if pd.notna(r) else "#d62728")
          for c, r in zip(top_val["in_consensus"], top_val["GF_rank"])]
axes[2].barh(top_val["gene"], top_val["mean_func_log2FC"], color=colors)
axes[2].set_xlabel("Mean functional log2FC (across gates)")
axes[2].set_title("Top functionally-validated CSC genes")
axes[2].axvline(0, color="k", lw=0.5)

plt.suptitle("Functional Ground-Truth Benchmark: DE vs Geneformer (sorted/sphere CSC assays)", fontsize=13)
plt.tight_layout()
plt.savefig("results/figures/F1_benchmark.png", dpi=120, bbox_inches="tight")
plt.close()

print("\n" + "=" * 64)
print("PHASE F1 COMPLETE")
print("=" * 64)
print("  Saved: results/tables/F1_functional_benchmark.csv")
print("  Saved: results/tables/F1_gene_level_validation.csv")
print("  Saved: results/figures/F1_benchmark.png")
print(f"  Metric wins across gates — DE: {wins['Stage1_DE']}  Geneformer: {wins['Stage2_Geneformer']}")
