# Phase I6 — Spatial CSC–macrophage–CD47 niche (TNBC Visium)
#
# Deepens the I4/I5 finding (TNBC CSCs upregulate CD47, the macrophage
# "don't eat me" ligand) with SPATIAL mechanism: within TNBC Visium sections,
# do CSC-high regions co-localize with macrophages AND CD47 — i.e. sit in
# CD47-high, macrophage-associated (phagocytosis-suppressed) niches?
#
# Per TNBC section: Spearman(stemness, macrophage), (stemness, CD47),
# (macrophage, CD47), and whether high-stemness+high-CD47 spots are
# macrophage-enriched.
#
# Output: results/tables/I6_macrophage_niche.csv
#         results/figures/I6_macrophage_niche.png

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

CSC_SIG   = ["CD44","SOX9","KLF4","MYC","VIM","EPCAM","FN1","SERPINE2","S100A6","ID3","HMGA1"]
MACRO     = ["CD68","CD163","LYZ","CSF1R","C1QA","C1QB","C1QC","ITGAM","MRC1","SIRPA","AIF1"]
TCELL     = ["CD3D","CD3E","CD2","CD8A","TRAC"]

def load_section(sec):
    p = f"{MTX_DIR}/{sec}_filtered_count_matrix"
    with open(f"{p}/matrix.mtx.gz", "rb") as fh:
        mat = scipy.io.mmread(fh)
    X = scipy.sparse.csr_matrix(mat).T.tocsr()
    feat = pd.read_csv(f"{p}/features.tsv.gz", header=None, sep="\t", compression=None)
    genes = feat.iloc[:, -1].astype(str).values
    barc = pd.read_csv(f"{p}/barcodes.tsv.gz", header=None, compression=None)[0].astype(str).values
    ad = sc.AnnData(X=X, obs=pd.DataFrame(index=barc), var=pd.DataFrame(index=genes))
    ad.var_names_make_unique()
    meta = pd.read_csv(f"{META_DIR}/{sec}_metadata.csv", index_col=0)
    shared = ad.obs_names.intersection(meta.index)
    ad = ad[shared].copy(); ad.obs = ad.obs.join(meta.loc[shared])
    return ad

def score(ad, genes, name):
    sc.tl.score_genes(ad, [x for x in genes if x in ad.var_names], score_name=name, ctrl_size=50)

print("=" * 64)
print("I6 — SPATIAL CSC–MACROPHAGE–CD47 NICHE (TNBC Visium)")
print("=" * 64)

sections = sorted({os.path.basename(d).replace("_filtered_count_matrix", "")
                   for d in glob.glob(f"{MTX_DIR}/*_filtered_count_matrix")})
rows = []
panels = []
for sec in sections:
    ad = load_section(sec)
    subtype = str(ad.obs["subtype"].iloc[0]) if "subtype" in ad.obs else "?"
    if subtype != "TNBC":
        continue
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    score(ad, CSC_SIG, "stemness"); score(ad, MACRO, "macro"); score(ad, TCELL, "tcell")
    cd47 = (np.asarray(ad[:, "CD47"].X.todense()).ravel() if "CD47" in ad.var_names
            else np.zeros(ad.n_obs))
    ad.obs["cd47"] = cd47

    r_sm, p_sm = spearmanr(ad.obs["stemness"], ad.obs["macro"])
    r_sc, p_sc = spearmanr(ad.obs["stemness"], ad.obs["cd47"])
    r_mc, p_mc = spearmanr(ad.obs["macro"], ad.obs["cd47"])
    # niche test: are high-stemness+high-CD47 spots macrophage-enriched?
    s_hi = ad.obs["stemness"] >= ad.obs["stemness"].quantile(0.75)
    c_hi = ad.obs["cd47"] >= np.quantile(ad.obs["cd47"], 0.75)
    niche = s_hi & c_hi
    rest = ~niche
    macro_niche = float(ad.obs.loc[niche, "macro"].mean())
    macro_rest  = float(ad.obs.loc[rest, "macro"].mean())
    try:
        _, p_niche = mannwhitneyu(ad.obs.loc[niche,"macro"], ad.obs.loc[rest,"macro"], alternative="greater")
    except Exception:
        p_niche = np.nan

    rows.append({"section": sec, "n_spots": ad.n_obs,
                 "stem~macro_r": round(r_sm,3), "stem~macro_p": p_sm,
                 "stem~CD47_r": round(r_sc,3), "stem~CD47_p": p_sc,
                 "macro~CD47_r": round(r_mc,3),
                 "macro_in_niche": round(macro_niche,3), "macro_rest": round(macro_rest,3),
                 "niche_macro_p": p_niche})
    panels.append((sec, ad.obs["stemness"].values, ad.obs["macro"].values, r_sm))
    print(f"  {sec}: n={ad.n_obs:4d}  stem~macro r={r_sm:+.3f}  stem~CD47 r={r_sc:+.3f}  "
          f"macro~CD47 r={r_mc:+.3f}  |  macro in CSC+CD47 niche {macro_niche:+.3f} vs {macro_rest:+.3f} (p={p_niche:.1e})")

res = pd.DataFrame(rows)
res.to_csv("results/tables/I6_macrophage_niche.csv", index=False)

print("\n" + "=" * 64)
print("SUMMARY (TNBC sections)")
print("=" * 64)
print(f"  stemness~macrophage: positive in {(res['stem~macro_r']>0).sum()}/{len(res)} "
      f"(mean r={res['stem~macro_r'].mean():+.3f})")
print(f"  stemness~CD47:       positive in {(res['stem~CD47_r']>0).sum()}/{len(res)} "
      f"(mean r={res['stem~CD47_r'].mean():+.3f})")
print(f"  CSC+CD47 niche macrophage-enriched in "
      f"{(res['macro_in_niche']>res['macro_rest']).sum()}/{len(res)} sections")

# figure
n = len(panels); fig, axes = plt.subplots(1, max(n,1), figsize=(5*max(n,1), 4.5))
axes = np.atleast_1d(axes)
for ax, (sec, stem, macro, r) in zip(axes, panels):
    ax.hexbin(stem, macro, gridsize=30, cmap="magma", mincnt=1)
    ax.set_title(f"{sec} (TNBC)\nstem~macrophage r={r:+.3f}", fontsize=10)
    ax.set_xlabel("stemness (spot)"); ax.set_ylabel("macrophage score (spot)")
plt.suptitle("Figure. Spatial co-localization of CSC state and macrophages (TNBC Visium)", fontsize=12)
plt.tight_layout()
plt.savefig("results/figures/I6_macrophage_niche.png", dpi=120, bbox_inches="tight")
plt.close()
print("\n  Saved: results/tables/I6_macrophage_niche.csv")
print("  Saved: results/figures/I6_macrophage_niche.png")
