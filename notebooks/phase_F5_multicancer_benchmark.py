# Phase F5 — Multi-cancer functional benchmark (Expansion 1)
#
# Extends the 4-method functional benchmark from 4 breast gates to 9 gates
# across 6 cancers and 3 functional criteria, by adding 5 ready RNA-seq
# functional datasets (no model retraining — the DE/Geneformer/LogReg/RF
# rankings are the fixed breast-derived rankings from F4).
#
# New gates:
#   GSE270565  Prostate    ALDH high/low       (RNA-seq, n=3/3)
#   GSE243840  Melanoma    ALDH1A3 high/low    (RNA-seq, n=3/3)
#   GSE228203  Prostate    sphere vs adherent  (RNA-seq, n=3/3)
#   GSE232783  Ovarian     stem vs adherent    (RNA-seq, n=3/3)
#   GSE166947  Bladder     ALDH+/-             (RNA-seq, n=3/3)
#
# Question: does the shortlist-vs-ranker dissociation (DE best precision,
# attention best genome-wide AUROC) hold across cancers, or is it breast-specific?
#
# Output: results/tables/F5_multicancer_benchmark.csv
#         results/figures/F5_multicancer_benchmark.png

import os, gzip, shutil, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
rng = np.random.default_rng(42)
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

print("=" * 66)
print("F5 — MULTI-CANCER FUNCTIONAL BENCHMARK (9 gates, 6 cancers)")
print("=" * 66)

# ─────────────────────────────────────────────────────────────────
# METHOD RANKINGS (fixed, breast-derived — from Stage 1/2 and F4)
# ─────────────────────────────────────────────────────────────────
s1 = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
s2 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
lr = pd.read_csv("results/tables/F4_logreg_ranking.csv")
rf = pd.read_csv("results/tables/F4_rf_ranking.csv")
lr_up = lr[lr["lr_coef"] > 0].sort_values("lr_coef", ascending=False)
rf_up = rf[rf["rf_signed"] > 0].sort_values("rf_importance", ascending=False)

METHODS = {
    "Stage1_DE":         (s1["gene_symbol"].tolist(),    dict(zip(s1["gene_symbol"], s1["wilcoxon_score"]))),
    "Stage2_Geneformer": (s2["gene_symbol"].tolist(),    dict(zip(s2["gene_symbol"], s2["attention_score"]))),
    "LogReg":            (lr_up["gene_symbol"].tolist(), dict(zip(lr_up["gene_symbol"], lr_up["lr_coef"]))),
    "RandomForest":      (rf_up["gene_symbol"].tolist(), dict(zip(rf_up["gene_symbol"], rf_up["rf_importance"]))),
}

# ─────────────────────────────────────────────────────────────────
# CONTRAST BUILDERS
# ─────────────────────────────────────────────────────────────────
def make_lfc(df, csc_cols, non_cols, symbols):
    """Per-gene functional log2FC = log2(mean(CSC)+1) - log2(mean(non)+1)."""
    csc = df[csc_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    non = df[non_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    lfc = np.log2(csc + 1) - np.log2(non + 1)
    s = pd.Series(lfc.values, index=pd.Index(symbols, dtype=str))
    s = s[~s.index.duplicated(keep="first")].dropna()
    s = s[~s.index.isin(["nan", "None", "NA", ""])]   # drop non-gene index labels
    return s

contrasts = {}   # name -> (cancer, criterion, Series)

# ── original 4 breast gates (rebuild exactly as F1/F4) ──
pdx = pd.read_excel("data/raw/GSE115302/matrix.xlsx")
pdx["gene"] = pdx["gene"].astype(str)
pdx = pdx.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
for c in ["MC1.bulk","MC1.aldh","MC1.cd44","VAR068.bulk","VAR068.aldh","VAR068.cd44"]:
    pdx[c] = pd.to_numeric(pdx[c], errors="coerce")
l2 = lambda n,d: np.log2((pdx[n]+1)/(pdx[d]+1))
s = ((l2("MC1.aldh","MC1.bulk")+l2("VAR068.aldh","VAR068.bulk"))/2); s.index = s.index.astype(str)
contrasts["Breast ALDH+ (GSE115302)"] = ("Breast","ALDH", s[~s.index.duplicated()])
s = ((l2("MC1.cd44","MC1.bulk")+l2("VAR068.cd44","VAR068.bulk"))/2); s.index = s.index.astype(str)
contrasts["Breast CD44+CD24- (GSE115302)"] = ("Breast","CD44/CD24", s[~s.index.duplicated()])
with gzip.open("data/raw/GSE182532/sphere_vs_adherent.txt.gz","rt") as f:
    sph = pd.read_csv(f)
if sph.shape[1]==1:
    sph = sph.iloc[:,0].str.split(",", expand=True)
    sph.columns = ["Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
sph["log2FoldChange"] = pd.to_numeric(sph["log2FoldChange"], errors="coerce")
sph = sph.dropna(subset=["Gene","log2FoldChange"]).drop_duplicates("Gene")
ss = pd.Series(sph["log2FoldChange"].values, index=sph["Gene"].astype(str));
contrasts["Breast mammosphere (GSE182532)"] = ("Breast","sphere", ss[~ss.index.duplicated()])
with gzip.open("data/raw/GSE36643/de_44hi24lo.csv.gz","rt") as f:
    g366 = pd.read_csv(f)
g366 = g366.dropna(subset=["gene symbol"]).drop_duplicates("gene symbol")
ss = pd.Series(np.log2(pd.to_numeric(g366["FC"],errors="coerce").clip(lower=1e-6).values),
               index=g366["gene symbol"].astype(str))
contrasts["Breast CD44hi/CD24lo (GSE36643)"] = ("Breast","CD44/CD24", ss[~ss.index.duplicated()])

# ── 5 new multi-cancer gates ──
# GSE270565 prostate ALDH: Pos=CSC, Neg=non; symbol=Gene.name
with gzip.open("data/raw/GSE270565/data.csv.gz","rt") as f:
    d = pd.read_csv(f)
contrasts["Prostate ALDH+ (GSE270565)"] = ("Prostate","ALDH",
    make_lfc(d, ["LNALPosNo1","LNALPosNo2","LNALPosNo3"],
                ["LNALNegNo1","LNALNegNo2","LNALNegNo3"], d["Gene.name"]))

# GSE243840 melanoma ALDH1A3: HIGH=CSC, LOW=non; symbol=geneName
with gzip.open("data/raw/GSE243840/tpm.csv.gz","rt") as f:
    d = pd.read_csv(f)
contrasts["Melanoma ALDH1A3+ (GSE243840)"] = ("Melanoma","ALDH",
    make_lfc(d, ["HIGH1","HIGH2","HIGH3"], ["LOW1","LOW2","LOW3"], d["geneName"]))

# GSE228203 prostate sphere: sphere=CSC, Adh=non; symbol=gene_name
with gzip.open("data/raw/GSE228203/fpkm.txt.gz","rt") as f:
    d = pd.read_csv(f, sep="\t")
contrasts["Prostate sphere (GSE228203)"] = ("Prostate","sphere",
    make_lfc(d, ["sphere_1","sphere_2","sphere_3"], ["Adh_1","Adh_2","Adh_3"], d["gene_name"]))

# GSE232783 ovarian: cols 4-6=Stem(CSC), 1-3=Adherent(non); symbol=gene_ID
with gzip.open("data/raw/GSE232783/data.xls.gz","rb") as fi, open("/tmp/gse232783.xls","wb") as fo:
    shutil.copyfileobj(fi, fo)
d = pd.read_excel("/tmp/gse232783.xls", engine="xlrd")
contrasts["Ovarian sphere (GSE232783)"] = ("Ovarian","sphere",
    make_lfc(d, ["4_FPKM","5_FPKM","6_FPKM"], ["1_FPKM","2_FPKM","3_FPKM"], d["gene_ID"]))

# GSE166947 bladder ALDH: H*=CSC, L*=non, exclude WildType; symbol parsed from Gene
with gzip.open("data/raw/GSE166947/cpm.txt.gz","rt") as f:
    d = pd.read_csv(f, sep="\t")
def seg(c):
    parts = c.split(".")
    return parts[2] if len(parts) > 2 else ""
h_cols = [c for c in d.columns if seg(c).startswith("H")]
l_cols = [c for c in d.columns if seg(c).startswith("L")]
sym = d["Gene"].astype(str).apply(lambda g: g.split("__")[1].rsplit("_",1)[0] if "__" in g else g)
contrasts["Bladder ALDH+ (GSE166947)"] = ("Bladder","ALDH", make_lfc(d, h_cols, l_cols, sym))

print(f"\nBuilt {len(contrasts)} functional gates:")
for name,(cancer,crit,s) in contrasts.items():
    print(f"  {name:38s} [{cancer:8s} {crit:9s}] {len(s)} genes")

# ─────────────────────────────────────────────────────────────────
# EVALUATE
# ─────────────────────────────────────────────────────────────────
TOPK, N_RAND = 100, 1000
def evaluate(ranked, smap, func):
    func = func.dropna()
    uni = list(func.index); uset = set(uni)
    vals = func.values                      # positional array (avoids label lookups)
    n = len(vals)
    topk = [g for g in ranked if g in uset][:TOPK]
    if len(topk) < 10: return None
    prec = float((func.loc[topk] > 0).mean())
    truth = (func > func.quantile(0.75)).astype(int)
    mscore = np.array([smap.get(g,0.0) for g in uni])
    auroc = roc_auc_score(truth.values, mscore) if truth.nunique()==2 else np.nan
    k = min(TOPK, n)
    rand = [(vals[rng.choice(n, k, replace=False)] > 0).mean() for _ in range(N_RAND)]
    return prec, float(np.mean(rand)), auroc

rows = []
for gate,(cancer,crit,func) in contrasts.items():
    for m,(ranked,smap) in METHODS.items():
        r = evaluate(ranked, smap, func)
        if r:
            rows.append({"gate":gate,"cancer":cancer,"criterion":crit,"method":m,
                         "precision@100":round(r[0],3),"random":round(r[1],3),
                         "AUROC":round(r[2],3)})
bench = pd.DataFrame(rows)
bench.to_csv("results/tables/F5_multicancer_benchmark.csv", index=False)

# ── per-gate print ──
for gate in contrasts:
    sub = bench[bench.gate==gate]
    print(f"\n── {gate} ──")
    for _,r in sub.iterrows():
        print(f"  {r['method']:18s} prec@100={r['precision@100']:.2f} (rand {r['random']:.2f})  AUROC={r['AUROC']}")

# ─────────────────────────────────────────────────────────────────
# SUMMARY: does the dissociation generalize?
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 66)
print("GENERALIZATION — mean across all 9 gates, and breast vs non-breast")
print("=" * 66)
def summarize(df, label):
    au = df.groupby("method")["AUROC"].mean().sort_values(ascending=False)
    pr = df.groupby("method")["precision@100"].mean().sort_values(ascending=False)
    print(f"\n[{label}]  best AUROC: {au.index[0]} ({au.iloc[0]:.3f}) | best precision: {pr.index[0]} ({pr.iloc[0]:.3f})")
    print("  AUROC:    " + "  ".join(f"{m}={v:.3f}" for m,v in au.items()))
    print("  prec@100: " + "  ".join(f"{m}={v:.2f}" for m,v in pr.items()))
    return au, pr
au_all, pr_all = summarize(bench, "ALL 9 gates")
summarize(bench[bench.cancer=="Breast"], "Breast only (4 gates)")
summarize(bench[bench.cancer!="Breast"], "Non-breast (5 gates, held-out cancers)")

# ─────────────────────────────────────────────────────────────────
# FIGURE
# ─────────────────────────────────────────────────────────────────
gates = list(contrasts.keys())
methods = list(METHODS.keys())
colors = {"Stage1_DE":"#d62728","Stage2_Geneformer":"#1f77b4","LogReg":"#2ca02c","RandomForest":"#ff7f0e"}
fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
x = np.arange(len(gates)); w = 0.2
for i,m in enumerate(methods):
    vals = [bench[(bench.gate==g)&(bench.method==m)]["AUROC"].values for g in gates]
    vals = [v[0] if len(v) else np.nan for v in vals]
    axes[0].bar(x+(i-1.5)*w, vals, w, label=m, color=colors[m])
axes[0].axhline(0.5, color="k", ls="--", lw=0.8)
axes[0].axvline(3.5, color="gray", ls=":", lw=1.2)
axes[0].text(1.5, 0.80, "Breast (discovery)", ha="center", fontsize=8, style="italic")
axes[0].text(6, 0.80, "Held-out cancers", ha="center", fontsize=8, style="italic")
axes[0].set_xticks(x); axes[0].set_xticklabels([g.split(" (")[0] for g in gates], rotation=35, ha="right", fontsize=7)
axes[0].set_ylabel("AUROC (genome-wide)"); axes[0].set_ylim(0,0.85)
axes[0].set_title("Genome-wide discrimination across 9 functional gates, 6 cancers")
axes[0].legend(fontsize=8, ncol=2)

# mean AUROC breast vs non-breast per method
bre = bench[bench.cancer=="Breast"].groupby("method")["AUROC"].mean()
non = bench[bench.cancer!="Breast"].groupby("method")["AUROC"].mean()
xm = np.arange(len(methods))
axes[1].bar(xm-0.2, [bre[m] for m in methods], 0.4, label="Breast (discovery)", color="#8888cc")
axes[1].bar(xm+0.2, [non[m] for m in methods], 0.4, label="Held-out cancers", color="#cc8888")
axes[1].axhline(0.5, color="k", ls="--", lw=0.8)
axes[1].set_xticks(xm); axes[1].set_xticklabels(methods, rotation=20, ha="right", fontsize=8)
axes[1].set_ylabel("mean AUROC"); axes[1].set_ylim(0,0.75)
axes[1].set_title("Does attention's ranking advantage generalize\nto held-out cancers?")
axes[1].legend(fontsize=8)
for i,m in enumerate(methods):
    axes[1].text(i-0.2, bre[m]+0.01, f"{bre[m]:.2f}", ha="center", fontsize=7)
    axes[1].text(i+0.2, non[m]+0.01, f"{non[m]:.2f}", ha="center", fontsize=7)

plt.suptitle("Figure. Multi-cancer functional benchmark (Expansion 1): 9 gates across 6 cancer types", fontsize=12)
plt.tight_layout()
plt.savefig("results/figures/F5_multicancer_benchmark.png", dpi=130, bbox_inches="tight")
plt.close()

print("\n" + "=" * 66)
print("PHASE F5 COMPLETE")
print("=" * 66)
print("  Saved: results/tables/F5_multicancer_benchmark.csv")
print("  Saved: results/figures/F5_multicancer_benchmark.png")
