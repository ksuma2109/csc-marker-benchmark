# Phase I2 — Spatial CSC–T-cell exclusion (Visium, Wu et al. 2021)
#
# Follow-up to I1: does CSC state spatially EXCLUDE T-cells? Uses the 6 Visium
# sections of the breast atlas (Zenodo 4739739), each with pathologist spot
# annotations ("Invasive cancer +/- lymphocytes"). Two tests per section:
#   (1) per-spot Spearman(stemness, T-cell score) — negative = anti-colocalized
#   (2) stemness in cancer spots WITH vs WITHOUT lymphocytes (annotation-based)
#
# No wet lab. Output:
#   results/tables/I2_spatial_exclusion.csv
#   results/figures/I2_spatial_exclusion.png

import os, glob, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io, scipy.sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

MTX_DIR  = "data/raw/spatial/filtered_count_matrices"
META_DIR = "data/raw/spatial/meta/metadata"

# CSC / stemness signature (consensus + canonical) and T-cell program
CSC_SIG = ["CD44","SOX9","KLF4","MYC","VIM","EPCAM","FN1","SERPINE2","S100A6",
           "ID3","HMGA1","ALDH1A3","PROM1"]
TCELL   = ["CD3D","CD3E","CD3G","CD2","CD8A","CD8B","TRAC","IL7R","CD247"]
EXHAUST = ["PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","TOX"]

def load_section(sec):
    # NB: the Zenodo files carry a .gz extension but are plain text — read raw.
    p = f"{MTX_DIR}/{sec}_filtered_count_matrix"
    with open(f"{p}/matrix.mtx.gz", "rb") as fh:
        mat = scipy.io.mmread(fh)                       # genes x spots
    X = scipy.sparse.csr_matrix(mat).T.tocsr()          # spots x genes
    feat = pd.read_csv(f"{p}/features.tsv.gz", header=None, sep="\t", compression=None)
    genes = feat.iloc[:, -1].astype(str).values         # last column = gene symbol
    barc = pd.read_csv(f"{p}/barcodes.tsv.gz", header=None, compression=None)[0].astype(str).values
    ad = sc.AnnData(X=X, obs=pd.DataFrame(index=barc), var=pd.DataFrame(index=genes))
    ad.var_names_make_unique()
    meta = pd.read_csv(f"{META_DIR}/{sec}_metadata.csv", index_col=0)
    shared = ad.obs_names.intersection(meta.index)
    ad = ad[shared].copy(); ad.obs = ad.obs.join(meta.loc[shared])
    return ad

def score(ad, genes, name):
    g = [x for x in genes if x in ad.var_names]
    sc.tl.score_genes(ad, g, score_name=name, ctrl_size=50)

print("=" * 64)
print("I2 — SPATIAL CSC–T-CELL EXCLUSION (Visium, 6 sections)")
print("=" * 64)

sections = sorted({os.path.basename(d).replace("_filtered_count_matrix", "")
                   for d in glob.glob(f"{MTX_DIR}/*_filtered_count_matrix")})
rows = []
panels = []
for sec in sections:
    ad = load_section(sec)
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    score(ad, CSC_SIG, "stemness"); score(ad, TCELL, "tcell"); score(ad, EXHAUST, "exhaust")
    subtype = str(ad.obs["subtype"].iloc[0]) if "subtype" in ad.obs else "?"

    # (1) per-spot anti-colocalization
    r, p = spearmanr(ad.obs["stemness"], ad.obs["tcell"])

    # (2) annotation-based: cancer spots WITH vs WITHOUT lymphocytes
    cls = ad.obs["Classification"].astype(str)
    is_cancer = cls.str.contains("cancer", case=False)
    with_lymph = is_cancer & cls.str.contains("lymphocyte", case=False)
    no_lymph   = is_cancer & ~cls.str.contains("lymphocyte", case=False)
    stem_with = ad.obs.loc[with_lymph, "stemness"]
    stem_no   = ad.obs.loc[no_lymph, "stemness"]
    if len(stem_with) >= 20 and len(stem_no) >= 20:
        u, pmw = mannwhitneyu(stem_no, stem_with, alternative="two-sided")
        delta = float(stem_no.mean() - stem_with.mean())   # >0 => lymph-poor cancer more stem-like
    else:
        pmw, delta = np.nan, np.nan

    rows.append({"section": sec, "subtype": subtype, "n_spots": ad.n_obs,
                 "spearman_stem_tcell": round(r, 3), "spearman_p": p,
                 "n_cancer_noLymph": int(no_lymph.sum()), "n_cancer_withLymph": int(with_lymph.sum()),
                 "stemness_delta_noLymph_minus_withLymph": round(delta, 3) if not np.isnan(delta) else np.nan,
                 "mw_p": pmw})
    panels.append((sec, subtype, ad.obs["stemness"].values, ad.obs["tcell"].values, r))
    print(f"  {sec} ({subtype}): n={ad.n_obs:4d}  stem~tcell r={r:+.3f} (p={p:.1e})  "
          f"cancer noLymph-withLymph Δstem={delta:+.3f} (p={pmw:.1e})" if not np.isnan(delta)
          else f"  {sec} ({subtype}): n={ad.n_obs:4d}  stem~tcell r={r:+.3f} (p={p:.1e})  [annotation test n/a]")

res = pd.DataFrame(rows)
res.to_csv("results/tables/I2_spatial_exclusion.csv", index=False)

# meta-summary
sig_excl = (res["spearman_stem_tcell"] < 0) & (res["spearman_p"] < 0.05)
print("\n" + "=" * 64)
print("SUMMARY")
print("=" * 64)
print(f"  Sections with significant NEGATIVE stemness–T-cell correlation "
      f"(spatial exclusion): {sig_excl.sum()}/{len(res)}")
print(f"  Mean stemness~tcell Spearman r across sections: {res['spearman_stem_tcell'].mean():+.3f}")
valid = res.dropna(subset=["stemness_delta_noLymph_minus_withLymph"])
if len(valid):
    print(f"  Cancer spots WITHOUT lymphocytes are more stem-like in "
          f"{(valid['stemness_delta_noLymph_minus_withLymph']>0).sum()}/{len(valid)} sections "
          f"(mean Δstemness={valid['stemness_delta_noLymph_minus_withLymph'].mean():+.3f})")

# ── figure: per-section stemness vs T-cell hexbin ──
n = len(panels); ncol = 3; nrow = int(np.ceil(n/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 5*nrow))
axes = np.array(axes).flatten()
for ax, (sec, st, stem, tc, r) in zip(axes, panels):
    ax.hexbin(stem, tc, gridsize=30, cmap="viridis", mincnt=1)
    ax.set_title(f"{sec} ({st})\nstem~T-cell r={r:+.3f}", fontsize=10)
    ax.set_xlabel("stemness (spot)"); ax.set_ylabel("T-cell score (spot)")
for ax in axes[n:]: ax.axis("off")
plt.suptitle("Figure. Spatial co-localization of CSC state and T-cells (Visium, per spot)", fontsize=13)
plt.tight_layout()
plt.savefig("results/figures/I2_spatial_exclusion.png", dpi=120, bbox_inches="tight")
plt.close()

print("\n  Saved: results/tables/I2_spatial_exclusion.csv")
print("  Saved: results/figures/I2_spatial_exclusion.png")
