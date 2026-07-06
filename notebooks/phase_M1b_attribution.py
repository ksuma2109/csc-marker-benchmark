# Phase M1b — Integrated-gradients gene attribution vs. attention ranking
#
# Paper 1's Geneformer gene ranking uses raw ATTENTION — a contested
# interpretability signal. Here we compute a more principled attribution:
# integrated gradients (IG) of the CSC-high logit w.r.t. the input gene-token
# embeddings, aggregated across CSC-high cells → a gradient-based gene
# importance ranking. We compare it to the attention ranking (G4) and DE.
#
# Question: does the model's causally-load-bearing gene importance (IG) agree
# with the attention ranking, or does the more rigorous method reorder it?
#
# HONEST FRAMING: IG reveals genes the MODEL relies on, not biological cause.
#
# Output: results/tables/M1b_ig_attribution.csv
#         results/figures/M1b_ig_vs_attention.png

import os, pickle, warnings
import numpy as np
import pandas as pd
import torch
from datasets import load_from_disk
from transformers import AutoModelForSequenceClassification
from captum.attr import LayerIntegratedGradients
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

os.makedirs("results/tables", exist_ok=True); os.makedirs("results/figures", exist_ok=True)
CKPT = "results/geneformer_finetuned/checkpoint-918"
MAXLEN, N_CELLS, N_STEPS = 512, 150, 20
torch.manual_seed(0)

print("=" * 64)
print("M1b — INTEGRATED-GRADIENTS ATTRIBUTION vs ATTENTION")
print("=" * 64)

dev = "mps" if torch.backends.mps.is_available() else "cpu"
model = AutoModelForSequenceClassification.from_pretrained(CKPT).eval().to(dev)
pad_id = model.config.pad_token_id or 0

# token_id -> ensembl -> symbol
with open("geneformer_repo/geneformer/token_dictionary_gc104M.pkl", "rb") as f:
    tok2ens = {v: k for k, v in pickle.load(f).items()}     # ensembl->tok inverted
g4 = pd.read_csv("results/tables/G4_geneformer_gene_ranking.csv").dropna(subset=["gene_symbol"])
ens2sym = dict(zip(g4["ensembl_id"], g4["gene_symbol"]))
def tok_to_symbol(t):
    ens = tok2ens.get(t)
    return ens2sym.get(ens)

# CSC-high cells (label 0 = csc_high)
ds = load_from_disk("results/geneformer_finetuned/brca_csc_labeled_test.dataset")
high_idx = [i for i in range(len(ds)) if ds[i]["label"] == 0][:N_CELLS]
print(f"CSC-high cells for IG: {len(high_idx)}  (target logit = CSC-high)")

emb_layer = model.bert.embeddings.word_embeddings
def fwd(input_ids, attention_mask):
    return model(input_ids, attention_mask=attention_mask).logits
lig = LayerIntegratedGradients(fwd, emb_layer)

# accumulate mean |IG| per gene
gene_attr, gene_cnt = {}, {}
print("Running integrated gradients (this takes a few minutes)...")
for k, i in enumerate(high_idx):
    seq = ds[i]["input_ids"][:MAXLEN]
    ids = torch.tensor([seq], device=dev)
    mask = torch.ones_like(ids)
    base = torch.full_like(ids, pad_id)
    attr = lig.attribute(inputs=ids, baselines=base,
                         additional_forward_args=(mask,), target=0, n_steps=N_STEPS)
    tok_attr = attr.sum(dim=-1).squeeze(0).abs().detach().cpu().numpy()   # per-token |IG|
    for pos, t in enumerate(seq):
        sym = tok_to_symbol(t)
        if sym is None: continue
        gene_attr[sym] = gene_attr.get(sym, 0.0) + float(tok_attr[pos])
        gene_cnt[sym]  = gene_cnt.get(sym, 0) + 1
    if (k+1) % 30 == 0: print(f"  {k+1}/{len(high_idx)} cells")

ig = pd.DataFrame([{"gene_symbol": g, "ig_mean": gene_attr[g]/gene_cnt[g], "n_cells": gene_cnt[g]}
                   for g in gene_attr if gene_cnt[g] >= 10])   # genes in >=10 cells
ig = ig.sort_values("ig_mean", ascending=False).reset_index(drop=True)
ig["ig_rank"] = np.arange(1, len(ig)+1)
ig.to_csv("results/tables/M1b_ig_attribution.csv", index=False)
print(f"\nGenes ranked by IG: {len(ig)}")
print("Top 15 IG genes:", list(ig["gene_symbol"].head(15)))

# ── compare to attention (G4) and DE ──
attn = g4[g4["gene_symbol"].notna()].copy()
attn["attn_rank"] = np.arange(1, len(attn)+1)
de = pd.read_csv("results/tables/A5_csc_markers_DE.csv")
de_rank = {g: i+1 for i, g in enumerate(de["gene_symbol"])}
attn_rank = dict(zip(attn["gene_symbol"], attn["attn_rank"]))

merged = ig.merge(attn[["gene_symbol", "attn_rank", "attention_score"]], on="gene_symbol", how="inner")
rho_ig_attn, _ = spearmanr(merged["ig_rank"], merged["attn_rank"])
top_ig = set(ig["gene_symbol"].head(50)); top_attn = set(attn["gene_symbol"].head(50))
overlap = len(top_ig & top_attn)

print("\n" + "=" * 64)
print("IG vs ATTENTION")
print("=" * 64)
print(f"  Spearman(IG rank, attention rank), shared genes: rho={rho_ig_attn:.3f}")
print(f"  Top-50 overlap (IG ∩ attention): {overlap}/50")
print(f"  IG-top genes NOT in attention top-50: {sorted(top_ig - top_attn)[:15]}")
print(f"  Attention-top genes NOT in IG top-50: {sorted(top_attn - top_ig)[:15]}")

# ── figure ──
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
top = ig.head(20)[::-1]
colors = ["#d62728" if g in top_attn else "#ff7f0e" for g in top["gene_symbol"]]
axes[0].barh(top["gene_symbol"], top["ig_mean"], color=colors)
axes[0].set_xlabel("mean |integrated gradient| (CSC-high logit)")
axes[0].set_title("Top 20 genes by IG attribution\n(red=also in attention top-50, orange=IG-unique)")

axes[1].scatter(merged["attn_rank"], merged["ig_rank"], s=12, alpha=0.5, color="#1f77b4")
axes[1].set_xlabel("attention rank (Paper 1)"); axes[1].set_ylabel("IG rank")
axes[1].set_title(f"IG vs attention gene ranking\nSpearman rho={rho_ig_attn:.3f}, top-50 overlap {overlap}/50")
axes[1].plot([0, merged["attn_rank"].max()], [0, merged["attn_rank"].max()], "k--", lw=0.6)
plt.suptitle("M1b — Integrated-gradients attribution vs attention (does the rigorous method reorder Paper 1's ranking?)",
             fontsize=12)
plt.tight_layout(); plt.savefig("results/figures/M1b_ig_vs_attention.png", dpi=130, bbox_inches="tight")
plt.close()
print("\n  Saved: results/tables/M1b_ig_attribution.csv")
print("  Saved: results/figures/M1b_ig_vs_attention.png")
