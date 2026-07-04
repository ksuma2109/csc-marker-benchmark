# Phase I7 — Independent-cohort validation of the TNBC CSC immune findings
#
# The decisive test: do the within-TNBC findings (I5) replicate in an
# INDEPENDENT cohort — Pal et al. 2021 (GSE161529, Visvader lab), distinct
# from the discovery set (Wu 2021 / GSE176078)? We use the 4 treatment-naive
# sporadic TNBC tumours ("Total cells", both compartments).
#
# Replication targets (within malignant/epithelial cells, per-cell Spearman
# of immune program vs stemness):
#   - MHC-I antigen presentation: expected POSITIVE (CSCs retain MHC-I)
#   - CD47:                        expected POSITIVE (innate evasion)
#   - HLA-E:                       expected POSITIVE
#   - CD274 (PD-L1):               expected ~0 (not the axis)
#
# Malignant epithelial cells are gated by EPCAM/KRT expression (in a TNBC
# tumour the epithelial compartment is overwhelmingly malignant); this avoids
# CNV inference for a validation.
#
# Output: results/tables/I7_validation_cohort.csv
#         results/figures/I7_validation_cohort.png

import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io, scipy.sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

DIR = "data/raw/GSE161529"
SAMPLES = {"GSM4909281":"TN-0126", "GSM4909282":"TN-0135",
           "GSM4909283":"TN-0106", "GSM4909284":"TN-0114"}

# same programs as I5 (discovery)
STEM   = ["CD44","SOX9","KLF4","MYC","VIM","EPCAM","FN1","SERPINE2","S100A6","ID3","HMGA1","ALDH1A3","PROM1"]
MHC_I  = ["HLA-A","HLA-B","HLA-C","TAP1","TAP2","NLRC5"]
EPI    = ["EPCAM","KRT8","KRT18","KRT19"]

feat = pd.read_csv(f"{DIR}/features.tsv.gz", header=None, sep="\t")
GENES = feat[1].astype(str).values   # gene symbols (col 2)

def load(gsm):
    m = scipy.io.mmread(f"{DIR}/{gsm}_matrix.mtx.gz")           # genes x cells (gzipped)
    X = scipy.sparse.csr_matrix(m).T.tocsr()
    bc = pd.read_csv(f"{DIR}/{gsm}_barcodes.tsv.gz", header=None)[0].astype(str).values
    ad = sc.AnnData(X=X, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=GENES))
    ad.var_names_make_unique()
    return ad

def gene_vec(ad, g):
    x = ad[:, g].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()

print("=" * 64)
print("I7 — INDEPENDENT VALIDATION (Pal 2021 / GSE161529, 4 TNBC)")
print("=" * 64)

rows = []
for gsm, name in SAMPLES.items():
    ad = load(gsm)
    # QC
    sc.pp.filter_cells(ad, min_genes=200)
    ad = ad[:, ad.var_names.notnull()].copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    # gate epithelial (malignant) cells
    sc.tl.score_genes(ad, [g for g in EPI if g in ad.var_names], score_name="epi", ctrl_size=50)
    epi = ad[ad.obs["epi"] > ad.obs["epi"].quantile(0.66)].copy()   # top-third epithelial
    if epi.n_obs < 100:
        print(f"  {name}: too few epithelial cells ({epi.n_obs}), skipped"); continue
    # scores within malignant epithelial cells
    sc.tl.score_genes(epi, STEM, score_name="stemness", ctrl_size=50)
    sc.tl.score_genes(epi, [g for g in MHC_I if g in epi.var_names], score_name="mhc", ctrl_size=50)
    stem = epi.obs["stemness"].values
    out = {"sample": name, "n_epi": epi.n_obs}
    for lab, sc_arr in [("MHC-I", epi.obs["mhc"].values),
                        ("CD47", gene_vec(epi, "CD47") if "CD47" in epi.var_names else None),
                        ("HLA-E", gene_vec(epi, "HLA-E") if "HLA-E" in epi.var_names else None),
                        ("PD-L1", gene_vec(epi, "CD274") if "CD274" in epi.var_names else None)]:
        if sc_arr is None:
            out[f"{lab}_r"] = np.nan; continue
        r, p = spearmanr(stem, sc_arr)
        out[f"{lab}_r"] = round(r, 3); out[f"{lab}_p"] = p
    rows.append(out)
    print(f"  {name} (n_epi={epi.n_obs:4d}):  "
          f"MHC-I r={out['MHC-I_r']:+.3f}  CD47 r={out.get('CD47_r',float('nan')):+.3f}  "
          f"HLA-E r={out.get('HLA-E_r',float('nan')):+.3f}  PD-L1 r={out.get('PD-L1_r',float('nan')):+.3f}")

res = pd.DataFrame(rows)
res.to_csv("results/tables/I7_validation_cohort.csv", index=False)

# ── replication verdict vs discovery (I5 within-TNBC signs: MHC-I+, CD47+, HLA-E+, PD-L1~0) ──
print("\n" + "=" * 64)
print("REPLICATION VERDICT (discovery within-TNBC: MHC-I+, CD47+, HLA-E+, PD-L1~0)")
print("=" * 64)
disc = {"MHC-I": +1, "CD47": +1, "HLA-E": +1}
for prog, sign in disc.items():
    col = f"{prog}_r"
    vals = res[col].dropna()
    same = int((np.sign(vals) == sign).sum())
    print(f"  {prog:6s}: discovery=+  |  validation mean r={vals.mean():+.3f}  "
          f"| same direction in {same}/{len(vals)} TNBC tumours")
pdl1 = res["PD-L1_r"].dropna()
print(f"  PD-L1 : discovery~0 |  validation mean r={pdl1.mean():+.3f} (expected near 0)")

# ── figure ──
fig, ax = plt.subplots(figsize=(9, 5))
progs = ["MHC-I","CD47","HLA-E","PD-L1"]
x = np.arange(len(progs)); w = 0.2
for i, (_, row) in enumerate(res.iterrows()):
    ax.bar(x + (i-1.5)*w, [row.get(f"{p}_r", np.nan) for p in progs], w, label=row["sample"])
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(progs)
ax.set_ylabel("Spearman r (program vs stemness, within malignant cells)")
ax.set_title("Independent-cohort validation (Pal 2021, 4 TNBC)\n"
             "Discovery expectation: MHC-I+, CD47+, HLA-E+, PD-L1~0")
ax.legend(fontsize=8, title="TNBC tumour")
plt.tight_layout()
plt.savefig("results/figures/I7_validation_cohort.png", dpi=130, bbox_inches="tight")
plt.close()
print("\n  Saved: results/tables/I7_validation_cohort.csv")
print("  Saved: results/figures/I7_validation_cohort.png")
