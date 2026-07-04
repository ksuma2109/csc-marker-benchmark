# Phase I9 — Second independent-cohort validation (Gao et al. 2020 / GSE148673)
#
# A second, orthogonal replication of the TNBC CSC immune findings, in 5 TNBC
# tumours from Gao et al. 2020 (Navin lab; copyKAT). Malignant cells are gated
# by the provided copyKAT prediction ("T" = aneuploid tumour) — a cleaner gate
# than the marker-based one used for the first validation (I7).
#
# Replication targets (per-cell Spearman vs stemness, malignant cells):
#   MHC-I retention (+), CD47 (+), HLA-E (+), PD-L1 (~0).
#
# Output: results/tables/I9_validation2_gao.csv

import os, gzip, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

STEM  = ["CD44","SOX9","KLF4","MYC","VIM","EPCAM","FN1","SERPINE2","S100A6","ID3","HMGA1","ALDH1A3","PROM1"]
MHC_I = ["HLA-A","HLA-B","HLA-C","TAP1","TAP2","NLRC5"]

print("=" * 64)
print("I9 — SECOND VALIDATION (Gao 2020 / GSE148673, 5 TNBC)")
print("=" * 64)

def gene_vec(ad, g):
    x = ad[:, g].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()

rows = []
for i in range(1, 6):
    with gzip.open(f"data/raw/GSE148673/TNBC{i}.txt.gz", "rt") as f:
        df = pd.read_csv(f, sep="\t", index_col=0)
    copykat = df.loc["copykat.pred"].astype(str)
    expr = df.drop(index="copykat.pred").apply(pd.to_numeric, errors="coerce").fillna(0)
    ad = sc.AnnData(X=expr.T.values.astype("float32"),
                    obs=pd.DataFrame({"copykat": copykat.values}, index=expr.columns),
                    var=pd.DataFrame(index=expr.index.astype(str)))
    ad.var_names_make_unique()
    mal = ad[ad.obs["copykat"] == "T"].copy()          # aneuploid = malignant
    if mal.n_obs < 100:
        print(f"  TNBC{i}: too few malignant cells ({mal.n_obs}), skipped"); continue
    sc.pp.normalize_total(mal, target_sum=1e4); sc.pp.log1p(mal)
    sc.tl.score_genes(mal, STEM, score_name="stemness", ctrl_size=50)
    sc.tl.score_genes(mal, [g for g in MHC_I if g in mal.var_names], score_name="mhc", ctrl_size=50)
    stem = mal.obs["stemness"].values
    out = {"sample": f"TNBC{i}", "n_malignant": mal.n_obs}
    for lab, arr in [("MHC-I", mal.obs["mhc"].values),
                     ("CD47",  gene_vec(mal, "CD47")  if "CD47"  in mal.var_names else None),
                     ("HLA-E", gene_vec(mal, "HLA-E") if "HLA-E" in mal.var_names else None),
                     ("PD-L1", gene_vec(mal, "CD274") if "CD274" in mal.var_names else None)]:
        if arr is None: out[f"{lab}_r"] = np.nan; continue
        r, _ = spearmanr(stem, arr); out[f"{lab}_r"] = round(r, 3)
    rows.append(out)
    print(f"  TNBC{i} (n_mal={mal.n_obs:4d}):  MHC-I r={out['MHC-I_r']:+.3f}  "
          f"CD47 r={out.get('CD47_r',float('nan')):+.3f}  HLA-E r={out.get('HLA-E_r',float('nan')):+.3f}  "
          f"PD-L1 r={out.get('PD-L1_r',float('nan')):+.3f}")

res = pd.DataFrame(rows)
res.to_csv("results/tables/I9_validation2_gao.csv", index=False)

print("\n" + "=" * 64)
print("REPLICATION VERDICT (Gao 2020) — discovery: MHC-I+, CD47+, HLA-E+, PD-L1~0")
print("=" * 64)
for prog, sign in [("MHC-I", 1), ("CD47", 1), ("HLA-E", 1)]:
    v = res[f"{prog}_r"].dropna()
    print(f"  {prog:6s}: mean r={v.mean():+.3f}  | same direction in {(np.sign(v)==sign).sum()}/{len(v)} tumours")
print(f"  PD-L1 : mean r={res['PD-L1_r'].dropna().mean():+.3f} (expected ~0)")
print("\n  Saved: results/tables/I9_validation2_gao.csv")
