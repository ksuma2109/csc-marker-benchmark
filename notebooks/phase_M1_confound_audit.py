# Phase M1 — Mechanistic interpretability: confound audit of the stemness transformer
#
# THE DECISION GATE. Does the fine-tuned Geneformer detect genuine stemness, or
# is its "stemness" signal confounded with proliferation (cell cycle) or subtype?
#
# We extract the model's per-layer [CLS] representation for held-out CSC-high /
# CSC-low cells and ask, by linear probing at each layer:
#   - is STEMNESS (the model's target) linearly decodable, and how does it build
#     across layers (pretrained input vs. fine-tuning)?
#   - are the CONFOUNDS (cell-cycle "Cancer Cycling"; subtype) decodable?
#   - can stemness be predicted from the confound LABELS alone (the confound
#     ceiling)? If the model's representation only matches this ceiling, its
#     stemness detector is confounded.
#
# HONEST FRAMING: this interprets the MODEL's representation, not biology.
#
# Output: results/tables/M1_confound_audit.csv
#         results/figures/M1_confound_audit.png

import os, warnings
import numpy as np
import pandas as pd
import torch
from datasets import load_from_disk
from transformers import AutoModelForSequenceClassification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

os.makedirs("results/tables", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)
CKPT = "results/geneformer_finetuned/checkpoint-918"
MAXLEN = 512               # matches training model_input_size
torch.manual_seed(0)

print("=" * 64)
print("M1 — CONFOUND AUDIT OF THE STEMNESS TRANSFORMER (decision gate)")
print("=" * 64)

dev = "mps" if torch.backends.mps.is_available() else "cpu"
model = AutoModelForSequenceClassification.from_pretrained(CKPT, output_hidden_states=True).eval().to(dev)
pad_id = model.config.pad_token_id or 0
n_layers = model.config.num_hidden_layers + 1  # + embeddings

ds = load_from_disk("results/geneformer_finetuned/brca_csc_labeled_test.dataset")
print(f"Held-out cells: {len(ds)}  | layers (incl. embeddings): {n_layers}")

# ── labels / confounds ──────────────────────────────────────────
label = np.array(ds["label"])                                   # 0=csc_high, 1=csc_low
y_stem = 1 - label                                              # 1 = CSC-high (intuitive)
minor = np.array(ds["celltype_minor"])
y_cyc  = (minor == "Cancer Cycling").astype(int)                # proliferation confound
subtype = np.array(ds["subtype"])
sub_map = {"ER+": 0, "TNBC": 1, "HER2+": 2}
y_sub = np.array([sub_map[s] for s in subtype])

# data-level confounding: base rates
print("\nData-level confounding (in held-out set):")
print(f"  CSC-high cells: {y_stem.mean()*100:.0f}%  | cycling cells: {y_cyc.mean()*100:.0f}%")
print(f"  P(cycling | CSC-high)={y_cyc[y_stem==1].mean():.2f}  vs  P(cycling | CSC-low)={y_cyc[y_stem==0].mean():.2f}")
for s,i in sub_map.items():
    hi = (y_sub[y_stem==1]==i).mean(); lo = (y_sub[y_stem==0]==i).mean()
    print(f"  {s}: {hi:.2f} of CSC-high vs {lo:.2f} of CSC-low")

# ── extract per-layer [CLS] hidden states ───────────────────────
print("\nExtracting hidden states...")
def batches(n, bs=16):
    for i in range(0, n, bs): yield range(i, min(i+bs, n))

cls_by_layer = [[] for _ in range(n_layers)]
for idx in batches(len(ds)):
    seqs = [ds[i]["input_ids"][:MAXLEN] for i in idx]
    m = max(len(s) for s in seqs)
    ids = torch.full((len(seqs), m), pad_id, dtype=torch.long)
    mask = torch.zeros((len(seqs), m), dtype=torch.long)
    for j, s in enumerate(seqs):
        ids[j, :len(s)] = torch.tensor(s); mask[j, :len(s)] = 1
    with torch.no_grad():
        out = model(ids.to(dev), attention_mask=mask.to(dev))
    for L in range(n_layers):
        cls_by_layer[L].append(out.hidden_states[L][:, 0, :].float().cpu().numpy())  # CLS = pos 0
cls_by_layer = [np.concatenate(x, 0) for x in cls_by_layer]
print(f"  extracted: {n_layers} layers x {cls_by_layer[0].shape}")

# ── per-layer linear probes ─────────────────────────────────────
def probe(X, y, multiclass=False):
    Xs = StandardScaler().fit_transform(X)
    clf = LogisticRegression(max_iter=2000, C=1.0)   # auto-handles multiclass
    return cross_val_score(clf, Xs, y, cv=5, scoring="accuracy").mean()

# confound-only ceiling: predict stemness from [cycling, subtype] labels alone
conf_feat = np.column_stack([y_cyc, (y_sub[:,None] == np.arange(3)).astype(int)])
conf_ceiling = probe(conf_feat, y_stem)
base_rate = max(y_stem.mean(), 1 - y_stem.mean())

print("\nPer-layer linear-probe accuracy:")
rows = []
for L in range(n_layers):
    X = cls_by_layer[L]
    a_stem = probe(X, y_stem)
    a_cyc  = probe(X, y_cyc)
    a_sub  = probe(X, y_sub, multiclass=True)
    rows.append({"layer": L, "stemness_acc": round(a_stem,3),
                 "cellcycle_acc": round(a_cyc,3), "subtype_acc": round(a_sub,3)})
    print(f"  layer {L:2d}: stemness={a_stem:.3f}  cell-cycle={a_cyc:.3f}  subtype={a_sub:.3f}")

res = pd.DataFrame(rows)
res.to_csv("results/tables/M1_confound_audit.csv", index=False)
final_stem = res["stemness_acc"].iloc[-1]

# ── verdict ─────────────────────────────────────────────────────
print("\n" + "=" * 64)
print("DECISION-GATE VERDICT")
print("=" * 64)
print(f"  Base rate (majority class):        {base_rate:.3f}")
print(f"  Stemness from confounds alone:     {conf_ceiling:.3f}  (confound ceiling)")
print(f"  Stemness from model repr (final):  {final_stem:.3f}")
margin = final_stem - conf_ceiling
print(f"  Margin over confound ceiling:      {margin:+.3f}")
if margin > 0.05:
    print("  -> Model captures stemness BEYOND cell-cycle/subtype confounds. Deep MI warranted.")
else:
    print("  -> Model's stemness signal is largely explained by confounds. Important caveat for Papers 1-2.")

# ── figure ─────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(res["layer"], res["stemness_acc"], "-o", color="#d62728", label="stemness (model target)")
ax.plot(res["layer"], res["cellcycle_acc"], "-s", color="#1f77b4", label="cell cycle (confound)")
ax.plot(res["layer"], res["subtype_acc"], "-^", color="#2ca02c", label="subtype (confound)")
ax.axhline(conf_ceiling, color="#d62728", ls="--", lw=1, label=f"stemness-from-confounds ceiling ({conf_ceiling:.2f})")
ax.axhline(base_rate, color="gray", ls=":", lw=1, label=f"base rate ({base_rate:.2f})")
ax.set_xlabel("layer (0 = token embeddings)"); ax.set_ylabel("linear-probe accuracy (5-fold CV)")
ax.set_title("M1 — What the stemness transformer represents, by layer\n"
             "(does stemness exceed cell-cycle / subtype confounds?)")
ax.legend(fontsize=8, loc="lower right"); ax.set_ylim(0.4, 1.0)
plt.tight_layout(); plt.savefig("results/figures/M1_confound_audit.png", dpi=130, bbox_inches="tight")
plt.close()
print("\n  Saved: results/tables/M1_confound_audit.csv")
print("  Saved: results/figures/M1_confound_audit.png")
