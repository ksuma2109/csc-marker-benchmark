# Phase G6 — Method Reliability Comparison
#
# 4 tests that ask: is Stage 1 (DE) or Stage 2 (Geneformer attention) more reliable?
#   1. Precision against gold-standard CSC markers from literature
#   2. Pathway enrichment — are genes enriched for known stemness pathways?
#   3. Pseudotime correlation — do genes track the stem-to-differentiated axis?
#   4. TCGA-BRCA survival — do signatures predict patient outcome?
#
# Input:  results/tables/A5_csc_markers_DE.csv          (Stage 1 — local)
#         results/tables/G4_geneformer_gene_ranking.csv  (Stage 2 — download from Drive)
#         results/tables/A6_gene_pseudotime_correlation.csv
# Output: results/figures/G6_reliability_*.png
#         results/tables/G6_reliability_summary.csv

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.stats as stats

warnings.filterwarnings("ignore")
os.makedirs("results/figures", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

TOP_N = 200  # compare top N genes from each method

# ─────────────────────────────────────────────────────────────────
# LOAD GENE LISTS
# ─────────────────────────────────────────────────────────────────
stage1_df   = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
stage1_genes = stage1_df["gene_symbol"].dropna().head(TOP_N).tolist()

G4_PATH = "results/tables/G4_geneformer_gene_ranking.csv"
if os.path.exists(G4_PATH):
    stage2_df    = pd.read_csv(G4_PATH).dropna(subset=["gene_symbol"])
    stage2_genes = stage2_df["gene_symbol"].head(TOP_N).tolist()
    print(f"Loaded Stage 2 top {len(stage2_genes)} genes from {G4_PATH}")
else:
    # Hardcoded top 20 from Colab output (fallback — full analysis needs Drive file)
    stage2_genes = [
        "SOX9","KLF4","FZD7","CD44","ALDH1A3","SERPINE2","MAGEA3","TFCP2L1",
        "CSAG1","CAPN6","FOXI1","GABRE","ATXN7","BMPR1B","HOXA4","TOX3",
        "MAGEA12","MME","ZC3H11A","SCNN1B",
    ]
    print(f"G4 CSV not found locally — using hardcoded top {len(stage2_genes)} genes.")
    print("Download G4_geneformer_gene_ranking.csv from Drive for full TOP-200 comparison.")

print(f"Stage 1: {len(stage1_genes)} genes  |  Stage 2: {len(stage2_genes)} genes\n")

# ─────────────────────────────────────────────────────────────────
# ANALYSIS 1 — KNOWN CSC MARKER PRECISION
# ─────────────────────────────────────────────────────────────────
print("=" * 60)
print("ANALYSIS 1 — GOLD-STANDARD CSC MARKER PRECISION")
print("=" * 60)

# Established CSC markers from peer-reviewed literature
# Sources: Bhatt et al. 2013; Prat et al. 2010; Liu et al. 2014;
#          Marjanovic et al. 2013; Gupta et al. 2009; Al-Hajj et al. 2003
GOLD_STANDARD = {
    # Canonical surface markers
    "CD44":   "canonical CSC surface marker (Al-Hajj 2003)",
    "CD24":   "low in CSCs — CD44hi/CD24lo = stem phenotype",
    "EPCAM":  "epithelial CSC marker",
    "CD133":  "PROM1 — CSC marker in multiple cancers",
    "ITGA6":  "alpha-6 integrin — breast CSC marker",
    "MME":    "CD10 — basal CSC marker",
    # Functional markers
    "ALDH1A1":"ALDH isoform — ALDEFLUOR assay basis",
    "ALDH1A3":"ALDH isoform — highest in breast CSCs",
    "ABCG2":  "drug efflux — side-population CSC assay",
    # Stemness TFs
    "SOX2":   "Yamanaka factor — pluripotency",
    "SOX9":   "master CSC TF in breast cancer",
    "KLF4":   "Yamanaka factor — pluripotency",
    "MYC":    "Yamanaka factor + oncogene",
    "NANOG":  "pluripotency master TF",
    "OCT4":   "POU5F1 — pluripotency master TF",
    "BMI1":   "polycomb — CSC self-renewal",
    # EMT / invasion
    "VIM":    "mesenchymal marker — EMT/invasive CSCs",
    "FN1":    "fibronectin — EMT and stemness niche",
    "TWIST1": "EMT TF — promotes CSC state",
    "SNAI1":  "SNAIL — EMT master TF",
    "CDH2":   "N-cadherin — EMT marker",
    # Signalling
    "FZD7":   "Wnt receptor — Wnt/β-catenin CSC self-renewal",
    "NOTCH1": "Notch signalling — CSC maintenance",
    "NOTCH3": "Notch signalling",
    "BMPR1B": "BMP receptor — CSC quiescence",
    # Immune evasion
    "CD274":  "PD-L1 — CSC immune evasion",
    "HLA-E":  "NK-cell evasion by CSCs",
    # Other validated
    "PROM1":  "CD133 alternative symbol",
    "NES":    "nestin — neural CSC marker",
    "HMGA1":  "chromatin remodeler — CSC programs",
    "ID3":    "inhibitor of differentiation — stemness",
}

gold_genes = set(GOLD_STANDARD.keys())

s1_hit = [g for g in stage1_genes if g in gold_genes]
s2_hit = [g for g in stage2_genes if g in gold_genes]
s1_prec = len(s1_hit) / len(stage1_genes) * 100
s2_prec = len(s2_hit) / len(stage2_genes) * 100

print(f"Gold-standard CSC marker set: {len(gold_genes)} genes")
print(f"\nStage 1 (DE)         : {len(s1_hit)}/{len(stage1_genes)} = {s1_prec:.1f}% precision")
for g in s1_hit:
    print(f"   ✓ {g:12s} — {GOLD_STANDARD[g]}")
print(f"\nStage 2 (Geneformer) : {len(s2_hit)}/{len(stage2_genes)} = {s2_prec:.1f}% precision")
for g in s2_hit:
    print(f"   ✓ {g:12s} — {GOLD_STANDARD[g]}")

# ─────────────────────────────────────────────────────────────────
# ANALYSIS 2 — PATHWAY ENRICHMENT
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("ANALYSIS 2 — PATHWAY ENRICHMENT (Enrichr / MSigDB Hallmarks)")
print("=" * 60)

try:
    import gseapy as gp

    gene_sets = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023"]

    results = {}
    for label, gene_list in [("Stage1_DE", stage1_genes), ("Stage2_Geneformer", stage2_genes)]:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism="human",
            outdir=None,
            verbose=False,
        )
        results[label] = enr.results

    # Stemness-relevant keywords to highlight
    STEM_KEYWORDS = [
        "stem", "wnt", "notch", "hedgehog", "tgf", "bmp", "epithelial_mesenchymal",
        "emt", "pluripoten", "self.renew", "sox", "myc", "klf",
    ]

    def is_stemness(term):
        t = term.lower().replace(" ", "_").replace("-", "_")
        return any(kw in t for kw in STEM_KEYWORDS)

    print("\nTop stemness-related enriched terms:")
    summary_rows = []
    for label in ["Stage1_DE", "Stage2_Geneformer"]:
        df = results[label].copy()
        df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(1e-300))
        stem_df = df[df["Term"].apply(is_stemness)].sort_values("-log10(padj)", ascending=False)
        top5 = stem_df.head(5)
        print(f"\n  {label}:")
        for _, row in top5.iterrows():
            print(f"    {row['Term'][:55]:55s}  padj={row['Adjusted P-value']:.2e}  "
                  f"overlap={row['Overlap']}")
            summary_rows.append({
                "method": label, "term": row["Term"],
                "padj": row["Adjusted P-value"],
                "neg_log10_padj": row["-log10(padj)"],
                "overlap": row["Overlap"],
            })

    # Plot: top 10 stemness terms per method
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    colors = {"Stage1_DE": "#d62728", "Stage2_Geneformer": "#1f77b4"}
    for ax, label in zip(axes, ["Stage1_DE", "Stage2_Geneformer"]):
        df  = results[label].copy()
        df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(1e-300))
        top = df[df["Term"].apply(is_stemness)].sort_values("-log10(padj)", ascending=False).head(10)
        if top.empty:
            ax.text(0.5, 0.5, "No stemness terms enriched", ha="center", va="center",
                    transform=ax.transAxes)
        else:
            short = [t[:45] for t in top["Term"][::-1]]
            ax.barh(short, top["-log10(padj)"][::-1], color=colors[label])
            ax.axvline(x=-np.log10(0.05), color="gray", linestyle="--", alpha=0.7, label="padj=0.05")
            ax.set_xlabel("-log10(adj. p-value)")
        ax.set_title(f"{label.replace('_', ' ')} — stemness terms", fontsize=11)
    plt.suptitle("Pathway Enrichment: stemness-related terms", fontsize=13)
    plt.tight_layout()
    plt.savefig("results/figures/G6_pathway_enrichment.png", dpi=120, bbox_inches="tight")
    plt.close()
    print("\n  Saved: results/figures/G6_pathway_enrichment.png")

except Exception as e:
    print(f"  gseapy enrichment failed: {e}")
    print("  (requires internet connection to Enrichr API)")

# ─────────────────────────────────────────────────────────────────
# ANALYSIS 3 — PSEUDOTIME CORRELATION
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("ANALYSIS 3 — PSEUDOTIME CORRELATION")
print("=" * 60)

pt_df = pd.read_csv("results/tables/A6_gene_pseudotime_correlation.csv")
pt_df = pt_df.set_index("gene")

def pseudotime_stats(gene_list, label):
    found   = [g for g in gene_list if g in pt_df.index]
    missing = [g for g in gene_list if g not in pt_df.index]
    corrs   = pt_df.loc[found, "spearman_r"].values
    mean_r  = np.mean(np.abs(corrs))          # mean |r| — magnitude of association
    pos_pct = (corrs > 0.05).mean() * 100     # % positively correlated (more expressed early/stem)
    neg_pct = (corrs < -0.05).mean() * 100    # % negatively correlated (more expressed late/diff)
    print(f"\n  {label}:")
    print(f"    Genes found in pseudotime table : {len(found)}/{len(gene_list)}")
    print(f"    Mean |Spearman r|               : {mean_r:.4f}")
    print(f"    Positive correlation (stem→)    : {pos_pct:.1f}%")
    print(f"    Negative correlation (→diff)    : {neg_pct:.1f}%")
    top5 = pt_df.loc[found].reindex(found).sort_values("spearman_r", ascending=False).head(5)
    print(f"    Top 5 stem-correlated genes:")
    for g, row in top5.iterrows():
        print(f"      {g:12s}  r={row['spearman_r']:+.3f}")
    return {"method": label, "n_found": len(found), "mean_abs_r": mean_r,
            "pct_positive": pos_pct, "pct_negative": neg_pct,
            "corrs": corrs, "found_genes": found}

s1_pt = pseudotime_stats(stage1_genes, "Stage 1 (DE)")
s2_pt = pseudotime_stats(stage2_genes, "Stage 2 (Geneformer)")

# Mann-Whitney test — are the two distributions significantly different?
if len(s1_pt["corrs"]) > 1 and len(s2_pt["corrs"]) > 1:
    stat, pval = stats.mannwhitneyu(
        np.abs(s1_pt["corrs"]), np.abs(s2_pt["corrs"]), alternative="two-sided"
    )
    print(f"\n  Mann-Whitney |r| comparison: U={stat:.0f}, p={pval:.4f}")
    winner = "Stage 2 (Geneformer)" if np.mean(np.abs(s2_pt["corrs"])) > np.mean(np.abs(s1_pt["corrs"])) \
             else "Stage 1 (DE)"
    print(f"  Higher mean |r|: {winner}")

# Plot: distribution of Spearman r values
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
for ax, res, color, label in [
    (axes[0], s1_pt, "#d62728", "Stage 1 (DE)"),
    (axes[1], s2_pt, "#1f77b4", "Stage 2 (Geneformer)"),
]:
    if len(res["corrs"]) == 0:
        ax.text(0.5, 0.5, "No genes found", ha="center", va="center", transform=ax.transAxes)
        continue
    ax.hist(res["corrs"], bins=30, color=color, alpha=0.75, edgecolor="white")
    ax.axvline(x=0, color="black", linewidth=1)
    ax.axvline(x=np.mean(res["corrs"]), color="gold", linewidth=2,
               label=f"mean={np.mean(res['corrs']):.3f}")
    ax.set_xlabel("Spearman r with pseudotime")
    ax.set_ylabel("Gene count")
    ax.set_title(f"{label}\nmean|r|={res['mean_abs_r']:.3f}  n={res['n_found']}")
    ax.legend(fontsize=9)
plt.suptitle("Pseudotime Correlation: does the gene track the stem→differentiated axis?",
             fontsize=12)
plt.tight_layout()
plt.savefig("results/figures/G6_pseudotime_correlation.png", dpi=120, bbox_inches="tight")
plt.close()
print("\n  Saved: results/figures/G6_pseudotime_correlation.png")

# ─────────────────────────────────────────────────────────────────
# ANALYSIS 4 — TCGA-BRCA SURVIVAL
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("ANALYSIS 4 — TCGA-BRCA SURVIVAL (cBioPortal API)")
print("=" * 60)

try:
    import requests
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    CBIO  = "https://www.cbioportal.org/api"
    STUDY = "brca_tcga"
    MPROF = "brca_tcga_rna_seq_v2_mrna"

    # Use top 15 genes from each method (enough for a robust signature score)
    sig_s1 = stage1_genes[:15]
    sig_s2 = stage2_genes[:15]
    all_sig_genes = list(set(sig_s1 + sig_s2))

    print(f"  Fetching expression for {len(all_sig_genes)} genes from cBioPortal...")

    # Resolve gene symbols → Entrez IDs
    import mygene
    mg   = mygene.MyGeneInfo()
    hits = mg.querymany(all_sig_genes, scopes="symbol",
                        fields="entrezgene", species="human", returnall=False)
    sym2entrez = {r["query"]: int(r["entrezgene"])
                  for r in hits if "entrezgene" in r and r["entrezgene"]}
    entrez_ids = list(sym2entrez.values())
    print(f"  Resolved {len(entrez_ids)}/{len(all_sig_genes)} gene symbols to Entrez IDs")

    # Fetch samples
    r = requests.get(f"{CBIO}/studies/{STUDY}/samples",
                     params={"pageSize": 2000, "pageNumber": 0}, timeout=30)
    r.raise_for_status()
    sample_ids = [s["sampleId"] for s in r.json()]
    print(f"  TCGA-BRCA samples: {len(sample_ids)}")

    # Fetch expression (POST with sample IDs + entrez IDs)
    payload = {"sampleIds": sample_ids, "entrezGeneIds": entrez_ids}
    r = requests.post(
        f"{CBIO}/molecular-profiles/{MPROF}/molecular-data/fetch",
        json=payload, params={"projection": "SUMMARY"}, timeout=60,
    )
    r.raise_for_status()
    expr_raw  = pd.DataFrame(r.json())
    # API returns entrezGeneId — pivot then rename columns back to symbols
    entrez2sym = {v: k for k, v in sym2entrez.items()}
    expr_wide  = expr_raw.pivot(index="sampleId", columns="entrezGeneId", values="value")
    expr_wide.columns = [entrez2sym.get(c, str(c)) for c in expr_wide.columns]
    expr_wide  = expr_wide.apply(pd.to_numeric, errors="coerce")
    print(f"  Expression matrix: {expr_wide.shape}")

    # Fetch survival (patient-level) via GET
    r = requests.get(
        f"{CBIO}/studies/{STUDY}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
        timeout=60,
    )
    r.raise_for_status()
    clin_all = pd.DataFrame(r.json())
    clin_all = clin_all[clin_all["clinicalAttributeId"].isin(["OS_MONTHS", "OS_STATUS"])]
    clin = clin_all.pivot(index="patientId", columns="clinicalAttributeId", values="value")
    clin["OS_MONTHS"] = pd.to_numeric(clin["OS_MONTHS"], errors="coerce")
    clin["OS_STATUS"] = clin["OS_STATUS"].str.contains("DECEASED", na=False).astype(int)

    # Map sample → patient: use first 12 chars of sampleId (TCGA barcode convention)
    expr_wide["patientId"] = expr_wide.index.str[:12]
    expr_wide = expr_wide.groupby("patientId").mean()

    # Merge expression + survival
    merged = expr_wide.join(clin, how="inner").dropna(subset=["OS_MONTHS", "OS_STATUS"])
    print(f"  Matched samples with survival: {len(merged)}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    results_4 = {}
    for ax, (label, sig) in zip(axes, [("Stage 1 (DE)", sig_s1),
                                        ("Stage 2 (Geneformer)", sig_s2)]):
        avail = [g for g in sig if g in merged.columns]
        if len(avail) < 3:
            ax.text(0.5, 0.5, f"Insufficient genes ({len(avail)})", ha="center",
                    va="center", transform=ax.transAxes)
            continue

        # Z-score each gene then average → signature score
        scored = merged[avail].apply(stats.zscore, nan_policy="omit")
        merged[f"sig_{label}"] = scored.mean(axis=1)

        # Median split
        median_score = merged[f"sig_{label}"].median()
        high = merged[f"sig_{label}"] >= median_score
        low  = ~high

        t_h = merged.loc[high, "OS_MONTHS"]
        e_h = merged.loc[high, "OS_STATUS"]
        t_l = merged.loc[low,  "OS_MONTHS"]
        e_l = merged.loc[low,  "OS_STATUS"]

        lr  = logrank_test(t_h, t_l, event_observed_A=e_h, event_observed_B=e_l)
        p   = lr.p_value

        kmf = KaplanMeierFitter()
        kmf.fit(t_h, e_h, label=f"High (n={high.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color="red")
        kmf.fit(t_l, e_l, label=f"Low (n={low.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color="blue")

        ax.set_title(f"{label}\nlog-rank p = {p:.2e}  ({len(avail)} genes)", fontsize=11)
        ax.set_xlabel("Months")
        ax.set_ylabel("Overall Survival")
        ax.legend(fontsize=9)
        results_4[label] = {"n_genes": len(avail), "logrank_p": p}

        print(f"\n  {label}: {len(avail)} genes used, log-rank p = {p:.4e}")

    plt.suptitle("TCGA-BRCA Overall Survival: high vs low CSC signature score",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig("results/figures/G6_tcga_survival.png", dpi=120, bbox_inches="tight")
    plt.close()
    print("\n  Saved: results/figures/G6_tcga_survival.png")

except Exception as e:
    print(f"  TCGA analysis failed: {e}")
    print("  Requires internet connection. Check cBioPortal is accessible.")

# ─────────────────────────────────────────────────────────────────
# SUMMARY TABLE
# ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("RELIABILITY COMPARISON SUMMARY")
print("=" * 60)

summary = pd.DataFrame([
    {
        "Analysis": "1. Gold-standard precision",
        "Stage1_DE": f"{len(s1_hit)}/{len(stage1_genes)} ({s1_prec:.1f}%)",
        "Stage2_Geneformer": f"{len(s2_hit)}/{len(stage2_genes)} ({s2_prec:.1f}%)",
        "Winner": "Stage 2" if s2_prec > s1_prec else "Stage 1" if s1_prec > s2_prec else "Tie",
    },
    {
        "Analysis": "3. Mean |pseudotime r|",
        "Stage1_DE": f"{s1_pt['mean_abs_r']:.4f}",
        "Stage2_Geneformer": f"{s2_pt['mean_abs_r']:.4f}",
        "Winner": "Stage 2" if s2_pt["mean_abs_r"] > s1_pt["mean_abs_r"] else "Stage 1",
    },
])
print(summary.to_string(index=False))
summary.to_csv("results/tables/G6_reliability_summary.csv", index=False)

print("\n✓ Phase G6 complete.")
print("  Figures in: results/figures/G6_*.png")
print("  Table:      results/tables/G6_reliability_summary.csv")
