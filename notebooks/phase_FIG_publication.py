# Composite publication figures — assembles existing analysis panels into
# 4 labeled multi-panel main figures for the manuscript.
#
# Figure 1  CSC identification (Stage 1): annotation, stemness, DE, trajectory
# Figure 2  Method comparison & functional benchmark (the centerpiece)
# Figure 3  Integrated candidate prioritization & pathway context
# Figure 4  Cross-cancer validation & survival
#
# Output: manuscript/figures/Figure{1..4}.png

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

FIG = "results/figures"
OUT = "manuscript/figures"
os.makedirs(OUT, exist_ok=True)

def panel(ax, path, letter, title=None):
    if os.path.exists(path):
        ax.imshow(mpimg.imread(path))
    else:
        ax.text(0.5, 0.5, f"[missing:\n{os.path.basename(path)}]",
                ha="center", va="center", fontsize=8)
    ax.axis("off")
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes,
            fontsize=20, fontweight="bold", va="top", ha="right")
    if title:
        ax.set_title(title, fontsize=10)

# ─────────────────────────────────────────────────────────────────
# FIGURE 1 — CSC identification (Stage 1 pipeline)
# ─────────────────────────────────────────────────────────────────
# NOTE: Figure 1 is generated separately and cleanly from data by
# notebooks/phase_FIG1_clean.py (single-panel regeneration). Not reassembled
# here to avoid overwriting the clean version with a tiled composite.

# ─────────────────────────────────────────────────────────────────
# FIGURE 2 — Method comparison & functional benchmark (centerpiece)
# F1_benchmark.png is already a 3-panel; present it as the main figure
# with the pathway-loading panel from F3 appended.
# ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 15))
gs = fig.add_gridspec(3, 1, hspace=0.14, height_ratios=[1, 1, 1])
panel(fig.add_subplot(gs[0]), f"{FIG}/F4_benchmark_4methods.png", "A")
panel(fig.add_subplot(gs[1]), f"{FIG}/F1_benchmark.png", "B")
panel(fig.add_subplot(gs[2]), f"{FIG}/F3_candidate_landscape.png", "C")
fig.suptitle("Figure 2. Functional benchmarking of CSC marker-ranking methods",
             fontsize=13, y=0.99)
fig.text(0.5, 0.005,
         "(A) Four-method comparison (DE, Geneformer attention, logistic regression, random forest) against "
         "ALDH+, CD44+CD24-, and mammosphere gates: genome-wide AUROC and top-marker precision. "
         "Attention exceeds both supervised baselines on AUROC; DE retains best precision. "
         "(B) DE vs Geneformer detail: precision vs random, AUROC, top functionally-validated genes. "
         "(C) Integrated candidate priority, attention-rank vs functional support, and CSC-pathway loading.",
         ha="center", fontsize=8, wrap=True)
fig.savefig(f"{OUT}/Figure2.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved {OUT}/Figure2.png")

# ─────────────────────────────────────────────────────────────────
# FIGURE 3 — Multi-cancer functional benchmark (Expansion 1)
# ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 6.5))
gs = fig.add_gridspec(1, 1)
panel(fig.add_subplot(gs[0]), f"{FIG}/F5_multicancer_benchmark.png", "")
fig.suptitle("Figure 3. Multi-cancer functional benchmark: rankings are tissue-specific",
             fontsize=13, y=1.02)
fig.text(0.5, -0.04,
         "Left: genome-wide AUROC of the four breast-derived rankings against nine functional gates "
         "spanning six cancers (breast = discovery; prostate/melanoma/ovarian/bladder = held-out). "
         "Right: mean AUROC within breast vs. held-out cancers — attention's within-breast lead does not "
         "transfer, indicating CSC marker rankings are substantially tissue-specific.",
         ha="center", fontsize=8, wrap=True)
fig.savefig(f"{OUT}/Figure3.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved {OUT}/Figure3.png")

# ─────────────────────────────────────────────────────────────────
# FIGURE 4 — Cross-cancer validation
# ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 11))
gs = fig.add_gridspec(2, 1, hspace=0.12)
panel(fig.add_subplot(gs[0]), f"{FIG}/B2_gbm_neoplastic.png", "A")
panel(fig.add_subplot(gs[1]), f"{FIG}/B3_3cancer_heatmap.png", "B")
fig.suptitle("Figure 4. Cross-cancer validation identifies a recurrent CSC program",
             fontsize=13, y=0.99)
fig.text(0.5, 0.005,
         "(A) GBM neoplastic-cell CSC analysis (GSE84465). "
         "(B) Three-cancer evidence matrix (breast, GBM, melanoma): CD44/VIM/FN1/MYC and "
         "ANXA1/CAV1/OSMR recur across tumor types.",
         ha="center", fontsize=8, wrap=True)
fig.savefig(f"{OUT}/Figure4.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved {OUT}/Figure4.png")

# ─────────────────────────────────────────────────────────────────
# FIGURE 5 — Survival (ssGSEA + subtype-stratified Cox)
# ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 6.5))
gs = fig.add_gridspec(1, 2, wspace=0.06, width_ratios=[1.4, 1])
panel(fig.add_subplot(gs[0]), f"{FIG}/G8_ssgsea_km.png",   "A")
panel(fig.add_subplot(gs[1]), f"{FIG}/G8_cox_forest.png",  "B")
fig.suptitle("Figure 5. Subtype-stratified survival association of CSC signatures (METABRIC, n=1,974)",
             fontsize=13, y=1.02)
fig.text(0.5, -0.03,
         "(A) Kaplan-Meier by ssGSEA stemness tertile. (B) Subtype-stratified Cox hazard ratios per SD: "
         "all signatures show the expected adverse-prognosis direction (HR>1); "
         "consensus signature significant for overall survival (HR=1.06, p=0.048).",
         ha="center", fontsize=8, wrap=True)
fig.savefig(f"{OUT}/Figure5.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved {OUT}/Figure5.png")

print("\nAll 5 main figures assembled in", OUT)
