# Phase F4 — Baseline methods + 4-way functional benchmark
#
# Addresses the "why only DE vs Geneformer?" question by adding two
# supervised baselines that isolate WHAT gives the transformer its edge:
#
#   Logistic regression (LR): supervised, LINEAR, no pretraining
#   Random forest (RF):       supervised, NONLINEAR, no pretraining
#   Stage1 DE:                unsupervised-contrast, univariate
#   Stage2 Geneformer:        pretrained transformer attention
#
# If Geneformer still wins the functional AUROC over LR and RF, its
# advantage is attributable to pretraining/context, not merely to being
# supervised or nonlinear.
#
# LR/RF are trained to classify csc_high vs csc_low cancer-epithelial cells
# (the same labels used for DE and Geneformer), then genes are ranked by
# |coefficient| (LR) and impurity importance (RF). All four rankings are
# scored against the three functional CSC gates from Phase F1.
#
# Output: results/tables/F4_logreg_ranking.csv
#         results/tables/F4_rf_ranking.csv
#         results/tables/F4_benchmark_4methods.csv
#         results/figures/F4_benchmark_4methods.png

import os, gzip, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
rng = np.random.default_rng(42)

os.makedirs("results/tables",  exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# ─────────────────────────────────────────────────────────────────
# BUILD BASELINE RANKINGS (LR, RF) FROM CSC-HIGH/LOW CELLS
# ─────────────────────────────────────────────────────────────────
print("=" * 64)
print("F4 — BASELINE METHODS + 4-WAY FUNCTIONAL BENCHMARK")
print("=" * 64)

print("\nLoading labeled cancer-epithelial cells...")
ad = sc.read_h5ad("data/processed/brca_A5_csc_scored.h5ad")
ad = ad[ad.obs["csc_label"].isin(["csc_high", "csc_low"])].copy()
y = (ad.obs["csc_label"] == "csc_high").astype(int).values
genes = np.array(ad.var_names)
X = ad.X  # log-normalized, sparse
print(f"  {ad.n_obs} cells ({int(y.sum())} high / {int((1-y).sum())} low), {len(genes)} genes")

# --- Logistic regression (linear, L2) ---
print("\nTraining logistic regression (linear baseline)...")
scaler = StandardScaler(with_mean=False)
Xs = scaler.fit_transform(X)
lr = LogisticRegression(penalty="l2", C=1.0, max_iter=1000, solver="saga", n_jobs=-1)
lr.fit(Xs, y)
lr_coef = lr.coef_.ravel()
lr_rank = pd.DataFrame({"gene_symbol": genes, "lr_coef": lr_coef,
                        "lr_abs": np.abs(lr_coef)}).sort_values("lr_coef", ascending=False)
lr_rank.to_csv("results/tables/F4_logreg_ranking.csv", index=False)
print(f"  Top LR genes (positive=CSC-high): {list(lr_rank['gene_symbol'].head(8))}")

# --- Random forest (nonlinear) ---
print("\nTraining random forest (nonlinear baseline)...")
Xd = X.toarray() if sp.issparse(X) else np.asarray(X)
rf = RandomForestClassifier(n_estimators=200, n_jobs=-1, max_depth=None,
                            max_features="sqrt", random_state=42)
rf.fit(Xd, y)
rf_imp = rf.feature_importances_
# sign the importance by direction of mean difference (high - low) so it ranks CSC-up genes
mean_hi = Xd[y == 1].mean(axis=0); mean_lo = Xd[y == 0].mean(axis=0)
rf_dir = np.sign(mean_hi - mean_lo)
rf_rank = pd.DataFrame({"gene_symbol": genes, "rf_importance": rf_imp,
                        "rf_signed": rf_imp * rf_dir}).sort_values("rf_importance", ascending=False)
rf_rank.to_csv("results/tables/F4_rf_ranking.csv", index=False)
print(f"  Top RF genes: {list(rf_rank['gene_symbol'].head(8))}")
del Xd  # free memory

# ─────────────────────────────────────────────────────────────────
# ASSEMBLE ALL FOUR METHOD RANKINGS
# ─────────────────────────────────────────────────────────────────
s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])

# For LR/RF, keep only CSC-UP genes (positive direction) as candidate markers,
# ranked by magnitude — matching the "up in CSC" framing of DE/Geneformer.
lr_up = lr_rank[lr_rank["lr_coef"] > 0].sort_values("lr_coef", ascending=False)
rf_up = rf_rank[rf_rank["rf_signed"] > 0].sort_values("rf_importance", ascending=False)

METHODS = {
    "Stage1_DE":        (s1["gene_symbol"].tolist(),      dict(zip(s1["gene_symbol"], s1["wilcoxon_score"]))),
    "Stage2_Geneformer":(s2["gene_symbol"].tolist(),      dict(zip(s2["gene_symbol"], s2["attention_score"]))),
    "LogReg":           (lr_up["gene_symbol"].tolist(),   dict(zip(lr_up["gene_symbol"], lr_up["lr_coef"]))),
    "RandomForest":     (rf_up["gene_symbol"].tolist(),   dict(zip(rf_up["gene_symbol"], rf_up["rf_importance"]))),
}
for m,(r,_) in METHODS.items():
    print(f"  {m}: {len(r)} ranked genes")

# ─────────────────────────────────────────────────────────────────
# FUNCTIONAL CONTRASTS (identical to Phase F1)
# ─────────────────────────────────────────────────────────────────
print("\nBuilding functional contrasts...")
contrasts = {}
pdx = pd.read_excel("data/raw/GSE115302/matrix.xlsx")
pdx["gene"] = pdx["gene"].astype(str)
pdx = pdx.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
for c in ["MC1.bulk","MC1.aldh","MC1.cd44","VAR068.bulk","VAR068.aldh","VAR068.cd44"]:
    pdx[c] = pd.to_numeric(pdx[c], errors="coerce")
l2 = lambda n,d: np.log2((pdx[n]+1)/(pdx[d]+1))
contrasts["ALDH+ (GSE115302)"]      = (l2("MC1.aldh","MC1.bulk")+l2("VAR068.aldh","VAR068.bulk"))/2
contrasts["CD44+CD24- (GSE115302)"] = (l2("MC1.cd44","MC1.bulk")+l2("VAR068.cd44","VAR068.bulk"))/2
with gzip.open("data/raw/GSE182532/sphere_vs_adherent.txt.gz","rt") as f:
    sph = pd.read_csv(f)
if sph.shape[1]==1:
    sph = sph.iloc[:,0].str.split(",", expand=True)
    sph.columns = ["Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
sph["log2FoldChange"] = pd.to_numeric(sph["log2FoldChange"], errors="coerce")
sph = sph.dropna(subset=["Gene","log2FoldChange"]).drop_duplicates("Gene")
contrasts["Mammosphere (GSE182532)"] = sph.set_index("Gene")["log2FoldChange"]
with gzip.open("data/raw/GSE36643/de_44hi24lo.csv.gz","rt") as f:
    g366 = pd.read_csv(f)
g366 = g366.dropna(subset=["gene symbol"]).drop_duplicates("gene symbol")
contrasts["CD44hi/CD24lo (GSE36643)"] = pd.Series(
    np.log2(pd.to_numeric(g366["FC"],errors="coerce").clip(lower=1e-6).values),
    index=g366["gene symbol"].astype(str).values)
for name in list(contrasts):
    s = contrasts[name]; s.index = s.index.astype(str)
    contrasts[name] = s[~s.index.duplicated(keep="first")]

# ─────────────────────────────────────────────────────────────────
# EVALUATE ALL FOUR METHODS
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 64)
print("4-WAY BENCHMARK vs FUNCTIONAL GATES")
print("=" * 64)
TOPK, N_RAND = 100, 1000

def evaluate(ranked, smap, func, gate, method):
    func = func.dropna(); uni = list(func.index); uset = set(uni)
    topk = [g for g in ranked if g in uset][:TOPK]
    if len(topk) < 10: return None
    prec = float((func.loc[topk] > 0).mean())
    truth = (func > func.quantile(0.75)).astype(int)
    mscore = np.array([smap.get(g,0.0) for g in uni])
    auroc = roc_auc_score(truth.values, mscore) if truth.nunique()==2 else np.nan
    rand = [ (func.loc[rng.choice(np.array(uni), min(TOPK,len(uni)), replace=False)]>0).mean()
             for _ in range(N_RAND)]
    return {"gate":gate,"method":method,f"precision@{TOPK}":round(prec,3),
            "random":round(float(np.mean(rand)),3),
            "enrichment":round(prec/np.mean(rand),2) if np.mean(rand)>0 else np.nan,
            "AUROC":round(auroc,3),
            "mean_func_log2FC":round(float(func.loc[topk].mean()),3)}

rows=[]
for gate,func in contrasts.items():
    for m,(ranked,smap) in METHODS.items():
        r = evaluate(ranked,smap,func,gate,m)
        if r: rows.append(r)
bench = pd.DataFrame(rows)
bench.to_csv("results/tables/F4_benchmark_4methods.csv", index=False)

for gate in contrasts:
    sub = bench[bench.gate==gate]
    if sub.empty: continue
    print(f"\n── {gate} ──")
    for _,r in sub.iterrows():
        print(f"  {r['method']:18s} precision@{TOPK}={r[f'precision@{TOPK}']:.2f} "
              f"({r['enrichment']}x rand)  AUROC={r['AUROC']}  meanLFC={r['mean_func_log2FC']:+.2f}")

# ── AUROC summary: the key comparison ──
print("\n" + "=" * 64)
print("GENOME-WIDE AUROC BY METHOD (mean across sorted gates)")
print("=" * 64)
sorted_gates = ["ALDH+ (GSE115302)","CD44+CD24- (GSE115302)","CD44hi/CD24lo (GSE36643)"]
auroc_summary = (bench[bench.gate.isin(sorted_gates)]
                 .groupby("method")["AUROC"].mean().sort_values(ascending=False))
for m,v in auroc_summary.items():
    print(f"  {m:18s} mean AUROC = {v:.3f}")
prec_summary = (bench[bench.gate.isin(sorted_gates)]
                .groupby("method")[f"precision@{TOPK}"].mean().sort_values(ascending=False))
print("\nPrecision@100 by method (mean across sorted gates):")
for m,v in prec_summary.items():
    print(f"  {m:18s} mean precision = {v:.3f}")

# ─────────────────────────────────────────────────────────────────
# FIGURE — 4-method comparison
# ─────────────────────────────────────────────────────────────────
gates = list(contrasts.keys())
methods = list(METHODS.keys())
colors = {"Stage1_DE":"#d62728","Stage2_Geneformer":"#1f77b4",
          "LogReg":"#2ca02c","RandomForest":"#ff7f0e"}
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

x = np.arange(len(gates)); w = 0.2
for i,m in enumerate(methods):
    vals = [bench[(bench.gate==g)&(bench.method==m)]["AUROC"].values for g in gates]
    vals = [v[0] if len(v) else np.nan for v in vals]
    axes[0].bar(x+(i-1.5)*w, vals, w, label=m, color=colors[m])
axes[0].axhline(0.5, color="k", ls="--", lw=0.8)
axes[0].set_xticks(x); axes[0].set_xticklabels([g.split(" (")[0] for g in gates], rotation=20, ha="right", fontsize=8)
axes[0].set_ylabel("AUROC (rank functional CSC genes)"); axes[0].set_ylim(0,0.85)
axes[0].set_title("Genome-wide discrimination of functional CSC genes\n(4 methods)")
axes[0].legend(fontsize=8, ncol=2)

for i,m in enumerate(methods):
    vals = [bench[(bench.gate==g)&(bench.method==m)][f"precision@{TOPK}"].values for g in gates]
    vals = [v[0] if len(v) else np.nan for v in vals]
    axes[1].bar(x+(i-1.5)*w, vals, w, label=m, color=colors[m])
rand_b = [bench[bench.gate==g]["random"].mean() for g in gates]
axes[1].plot(x, rand_b, "k--o", markersize=4, label="random")
axes[1].set_xticks(x); axes[1].set_xticklabels([g.split(" (")[0] for g in gates], rotation=20, ha="right", fontsize=8)
axes[1].set_ylabel(f"precision@{TOPK}"); axes[1].set_ylim(0,1)
axes[1].set_title("Top-marker functional precision\n(4 methods)")
axes[1].legend(fontsize=8, ncol=2)

plt.suptitle("Figure. Four-method functional benchmark: transformer attention vs. DE, logistic regression, random forest",
             fontsize=12)
plt.tight_layout()
plt.savefig("results/figures/F4_benchmark_4methods.png", dpi=130, bbox_inches="tight")
plt.close()

print("\n" + "=" * 64)
print("PHASE F4 COMPLETE")
print("=" * 64)
print("  Saved: results/tables/F4_benchmark_4methods.csv")
print("  Saved: results/figures/F4_benchmark_4methods.png")
print(f"  Best mean AUROC: {auroc_summary.index[0]} ({auroc_summary.iloc[0]:.3f})")
print(f"  Best mean precision: {prec_summary.index[0]} ({prec_summary.iloc[0]:.3f})")
