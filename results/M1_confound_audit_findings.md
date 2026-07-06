# M1 — Confound Audit of the Stemness Transformer: Findings (MI decision gate)

**Status:** standalone MI analysis (Phase M1 of the interpretability project). NOT in either manuscript.
**Model:** fine-tuned Geneformer (checkpoint-918, BertForSequenceClassification, 12 layers, 768-dim), CSC-high/low classifier. **Data:** 1,225 held-out cancer-epithelial cells.
**Script:** `notebooks/phase_M1_confound_audit.py` · **Outputs:** `results/tables/M1_confound_audit.csv`, `results/figures/M1_confound_audit.png`.

## Question (the decision gate)
Does the fine-tuned Geneformer detect genuine stemness, or is its "stemness" signal confounded with proliferation (cell cycle) or molecular subtype? Method: extract the model's per-layer [CLS] representation and linearly probe it for stemness, cell-cycle ("Cancer Cycling"), and subtype; compare to a confound-only ceiling (predict stemness from cycling + subtype labels alone).

## Honest framing
This interprets what the **model's representation** encodes, not biological causation, and the label is a stemness *signature*, not functional CSCs.

## Key results

### Data-level confounding (held-out set)
- CSC-high vs CSC-low: **69% vs 29% TNBC**, 27% vs 62% ER+ → stemness label is strongly entangled with subtype.
- Cycling cells: P(cycling|high)=0.29 vs P(cycling|low)=0.17 → modest proliferation confounding.

### Per-layer linear-probe accuracy
| Layer | Stemness | Cell-cycle | Subtype |
|---|---|---|---|
| 0 (embeddings) | **0.50 (chance)** | 0.77 | 0.49 |
| final (12) | **0.91** | 0.88 | **0.98** |

- Stemness is **absent from the input token embeddings** (chance) and **built up by fine-tuning** — a constructed classification feature, not a natural pretrained one.
- Cell-cycle **is** present in the pretrained embeddings (0.77); subtype becomes near-perfectly decodable (0.98) from layer 1 on.

### Decision-gate verdict
| | Accuracy |
|---|---|
| Base rate (majority class) | 0.50 |
| Stemness from confounds alone (ceiling) | **0.70** |
| Stemness from model representation | **0.91** |
| **Margin over confound ceiling** | **+0.21** |

## Interpretation — two things are simultaneously true
1. **The model captures genuine, confound-independent stemness signal** (0.91 vs 0.70 ceiling; +0.21 margin). There is a real stemness representation to interpret → **deeper MI (M2/M3) is warranted.**
2. **But the stemness construct is substantially confounded** — ~70% of it is predictable from subtype + cell-cycle alone, and subtype is near-perfectly represented (0.98). A large fraction of the model's "stemness detector" is effectively a **TNBC/proliferation detector**, with genuine stemness layered on top (~21 points).

## Why this matters (connects to Papers 1–2)
This **quantifies** caveats already stated qualitatively in both manuscripts:
- Paper 1's finding that marker rankings are **tissue-specific**, and Paper 2's **subtype-confounding** caveat, both follow from the stemness construct being ~70% explainable by subtype/proliferation.
- Concretely: the CSC marker rankings partly reflect subtype (TNBC) and cell-cycle biology, not pure stemness. M1 puts a number on the split (~70% confound-aligned, ~21 points stemness-specific).
- This is an honest, publishable strengthening of both papers' limitations sections, and it is the kind of finding MI is well-suited to produce.

## Verdict for the MI project
**PROCEED** to M2 (causal localization) / M3 (feature discovery) — the model has a distinct, non-confounded stemness component worth reverse-engineering. But any downstream MI claim must be read against the substantial subtype/proliferation confounding documented here.

## Caveats
- Linear probes only (nonlinear structure not tested); [CLS]-pooled representation.
- Uses checkpoint-918 (epoch 1, F1 92.95%), not the best checkpoint-3676 (on Drive) — results expected to be similar.
- "Stemness" = signature-derived label; confound audit bounds interpretation but does not make it causal.
