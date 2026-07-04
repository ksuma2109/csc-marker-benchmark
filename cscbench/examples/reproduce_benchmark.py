"""Reproduce the study's multi-cancer functional benchmark using the cscbench API.

Rebuilds a subset of the paper's benchmark (four methods x functional gates
across cancers) from the processed functional matrices, demonstrating the
package on real data. Paths are relative to the project root; run as:

    python cscbench/examples/reproduce_benchmark.py

Requires the processed functional datasets in data/raw/ (see paper Data
Availability) and the method rankings in results/tables/.
"""
import os, gzip, shutil
import numpy as np
import pandas as pd
from cscbench import FunctionalGate, run_benchmark

ROOT = os.environ.get("CSC_PROJECT_ROOT", ".")

# ── method rankings (gene -> score) ──────────────────────────────────────────
def load_rankings():
    t = lambda p: pd.read_csv(os.path.join(ROOT, "results/tables", p))
    s1 = t("A5_csc_markers_DE.csv")
    s2 = t("G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
    lr = t("F4_logreg_ranking.csv"); rf = t("F4_rf_ranking.csv")
    lr = lr[lr["lr_coef"] > 0]; rf = rf[rf["rf_signed"] > 0]
    return {
        "DE":           dict(zip(s1["gene_symbol"], s1["wilcoxon_score"])),
        "Geneformer":   dict(zip(s2["gene_symbol"], s2["attention_score"])),
        "LogReg":       dict(zip(lr["gene_symbol"], lr["lr_coef"])),
        "RandomForest": dict(zip(rf["gene_symbol"], rf["rf_importance"])),
    }

# ── functional gates ─────────────────────────────────────────────────────────
def load_gates():
    raw = lambda p: os.path.join(ROOT, "data/raw", p)
    gates = []

    # Prostate ALDH (GSE270565)
    with gzip.open(raw("GSE270565/data.csv.gz"), "rt") as f:
        d = pd.read_csv(f)
    gates.append(FunctionalGate.from_matrix(
        d, ["LNALPosNo1", "LNALPosNo2", "LNALPosNo3"],
        ["LNALNegNo1", "LNALNegNo2", "LNALNegNo3"],
        name="Prostate ALDH+ (GSE270565)", cancer="Prostate", criterion="ALDH",
        gene_col="Gene.name"))

    # Melanoma ALDH1A3 (GSE243840)
    with gzip.open(raw("GSE243840/tpm.csv.gz"), "rt") as f:
        d = pd.read_csv(f)
    gates.append(FunctionalGate.from_matrix(
        d, ["HIGH1", "HIGH2", "HIGH3"], ["LOW1", "LOW2", "LOW3"],
        name="Melanoma ALDH1A3+ (GSE243840)", cancer="Melanoma", criterion="ALDH",
        gene_col="geneName"))

    # Prostate sphere (GSE228203)
    with gzip.open(raw("GSE228203/fpkm.txt.gz"), "rt") as f:
        d = pd.read_csv(f, sep="\t")
    gates.append(FunctionalGate.from_matrix(
        d, ["sphere_1", "sphere_2", "sphere_3"], ["Adh_1", "Adh_2", "Adh_3"],
        name="Prostate sphere (GSE228203)", cancer="Prostate", criterion="sphere",
        gene_col="gene_name"))

    # Ovarian sphere (GSE232783) — old-format .xls
    with gzip.open(raw("GSE232783/data.xls.gz"), "rb") as fi, open("/tmp/_ov.xls", "wb") as fo:
        shutil.copyfileobj(fi, fo)
    d = pd.read_excel("/tmp/_ov.xls", engine="xlrd")
    gates.append(FunctionalGate.from_matrix(
        d, ["4_FPKM", "5_FPKM", "6_FPKM"], ["1_FPKM", "2_FPKM", "3_FPKM"],
        name="Ovarian sphere (GSE232783)", cancer="Ovarian", criterion="sphere",
        gene_col="gene_ID"))

    # Bladder ALDH (GSE166947) — parse symbol from "ENSG__SYMBOL_strand", exclude WildType
    with gzip.open(raw("GSE166947/cpm.txt.gz"), "rt") as f:
        d = pd.read_csv(f, sep="\t")
    seg = lambda c: c.split(".")[2] if len(c.split(".")) > 2 else ""
    h = [c for c in d.columns if seg(c).startswith("H")]
    l = [c for c in d.columns if seg(c).startswith("L")]
    d = d.copy()
    d["sym"] = d["Gene"].astype(str).apply(
        lambda g: g.split("__")[1].rsplit("_", 1)[0] if "__" in g else g)
    gates.append(FunctionalGate.from_matrix(
        d, h, l, name="Bladder ALDH+ (GSE166947)", cancer="Bladder", criterion="ALDH",
        gene_col="sym"))
    return gates


def main():
    rankings, gates = load_rankings(), load_gates()
    print(f"{len(rankings)} methods x {len(gates)} gates\n")
    table = run_benchmark(rankings, gates, topk=100, n_random=1000, random_state=42)
    pd.set_option("display.width", 160)
    print(table.to_string(index=False))

    print("\nMean across held-out cancers (AUROC — the key comparison):")
    print(table.groupby("method")["auroc"].mean().sort_values(ascending=False).round(3).to_string())
    print("\nMean precision@100:")
    print(table.groupby("method")["precision_at_k"].mean().sort_values(ascending=False).round(3).to_string())


if __name__ == "__main__":
    main()
