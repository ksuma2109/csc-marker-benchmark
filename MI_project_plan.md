# Project Plan — Mechanistic Interpretability of a Cancer-Stemness Transformer

**Working title:** *What does a transformer learn about cancer stemness? A mechanistic interpretation of a fine-tuned single-cell foundation model.*

**Status:** planning (a potential third project / extension to Paper 1). Not started.

---

## 0. The honest scope — read this first

Mechanistic interpretability (MI) reverse-engineers what a neural network *computes*. Applied to our fine-tuned Geneformer, it reveals **which genes, features, and internal circuits the model uses to call a cell stemness-high** — a statement about the *learned function on our data*, **not** about what biologically *causes* stemness.

Two consequences that must frame every claim in this project:
1. **Model-internal ≠ biological cause.** MI outputs are refined hypotheses for perturbation, never mechanisms. Only wet-lab perturbation (CRISPR) certifies cause.
2. **The label is a stemness *signature*, not functional CSCs.** So MI reveals "what the model uses to predict signature-high cells" — doubly removed from real stemness biology. MI on a model trained on correlational labels can *look* mechanistic while being circular. Every deliverable must be phrased as "the model represents/uses X," not "X causes stemness."

The correct goal: **reverse-engineer the model's stemness detector, audit it for confounds, and generate better-prioritized causal hypotheses** — an upgrade over the raw attention ranking Paper 1 rests on.

---

## 1. Research questions

- **Q1 (confound audit):** Does the model detect genuine stemness, or is its "stemness" signal confounded with cell cycle / proliferation / subtype?
- **Q2 (attribution):** Which input genes are *causally load-bearing for the model's prediction* (vs. merely high-attention)? How does this compare to the attention and DE rankings?
- **Q3 (localization):** Which layers, positions, and attention heads carry the stemness computation?
- **Q4 (features/programs):** Does the model use interpretable latent *programs* (e.g., Wnt-high + differentiation-low) rather than single genes, and which gene *interactions* does it exploit — the contextual signal that is the model's claimed advantage over DE?
- **Q5 (grounding):** Do the model's load-bearing genes/programs overlap with the functionally-validated ones (Paper 1 F1 benchmark, DepMap F6)?

---

## 2. Assets we already have

- **Fine-tuned Geneformer** (Geneformer-V2-104M_CLcancer, CSC-high/low classifier; best checkpoint-3676, macro-F1 93.2% on Google Drive; checkpoint-918/919 local).
- **Tokenized CSC-high / CSC-low dataset** (12,246 cancer-epithelial cells, balanced).
- **Attention-based gene ranking** (G4) — the baseline to beat/compare.
- **Functional benchmark + CRISPR results** (F1, F6) — for biological grounding (Q5).
- **Per-cell metadata**: stemness score, subtype, cell-cycle-scoreable genes — needed for the confound probes.

---

## 3. Phased plan

### Phase M1 — Attribution + confound audit  ·  *decision gate*  ·  ~2–4 days
The cheap, high-value entry point; determines whether deeper MI is worthwhile.
- **Integrated gradients** (Captum or manual) attributing the CSC-high logit to input gene tokens. Compare to the attention ranking and DE (rank correlation; which genes appear only under gradient attribution).
- **Linear probes** on hidden states at each layer for: stemness score, **cell-cycle/proliferation** (MKI67, TOP2A, PCNA…), **subtype**, ribosomal/translation load. *Question: is stemness linearly decodable and distinct from these confounds, or is the model's "stemness" mostly proliferation/subtype?*
- **Deliverable / gate:** (a) a gradient-based gene ranking; (b) a verdict on whether the model detects genuine stemness vs. a confound. **If the signal is mostly cell-cycle/subtype → that's an important, honest finding that qualifies Paper 1, and deep MI is lower priority. If stemness is distinct → proceed.**

### Phase M2 — Causal localization (activation patching)  ·  ~1–2 weeks
- **Activation/attention patching** between matched CSC-high and CSC-low cells: patch layer activations / specific gene positions and measure the effect on the stemness logit → identify the load-bearing layers and genes (causal *within the model*).
- **Attention-head ablation**: knock out heads, measure ΔF1 → which heads implement the classification.
- **Deliverable:** a minimal "circuit" (layers, heads, gene positions) computing the stemness call.

### Phase M3 — Feature & interaction discovery  ·  ~weeks (frontier)
- **Sparse autoencoders (SAEs)** on residual-stream activations → interpretable latent features; test whether features map to gene programs (Wnt, EMT, cell cycle).
- **Gene-interaction analysis**: which gene *pairs/contexts* the model relies on (gene A matters only when B is high) — directly probes the contextual advantage over DE.
- **Deliverable:** latent stemness programs + a set of model-used gene interactions.

### Phase M4 — Biological grounding & hypotheses  ·  ~days
- Cross-reference M1–M3 load-bearing genes/programs with the **functional benchmark (F1)** and **CRISPR dependency (F6)**: do the model's mechanistic features overlap functionally validated genes?
- Produce a ranked list of **perturbation hypotheses** (for the wet-lab collaborators already scouted).

---

## 4. Technical requirements & tooling

- **Model access:** fine-tuned checkpoint + tokenizer + the tokenized dataset. (Retrieve checkpoint-3676 from Drive for the definitive model.)
- **Hooks:** HuggingFace `transformers` forward hooks for activations/attention (Geneformer is a BERT-style `BertForSequenceClassification` — standard hooks work). Note: TransformerLens targets GPT-style models, so patching/SAE code will need light custom adaptation for BERT-style Geneformer.
- **Attribution:** Captum (integrated gradients, layer attribution).
- **Probing:** scikit-learn linear models on extracted hidden states.
- **SAEs:** a small SAE training loop (PyTorch) on cached activations.
- **Compute:** M1–M2 feasible on CPU/MPS (inference-only, slow but fine); M3 (SAEs) benefits from GPU/Colab.

---

## 5. Risks & failure modes (name them up front)

- **Circularity / label problem** (Section 0.2) — the dominant risk; frame all outputs as model-internal.
- **Over-interpretation** — MI results *look* mechanistic; guard with ablations and the biological-grounding step (M4).
- **Confound reveal** — M1 may show the model rides cell-cycle/subtype. Deflating but honest and publishable; it directly informs Papers 1–2.
- **Tooling immaturity** — single-cell foundation-model MI is nascent; budget time for custom patching/SAE code.
- **Attention-baseline weakness** — since Paper 1 already used attention, the novelty must come from the *more rigorous* methods (gradients, patching, SAEs) and the *confound audit*, not from re-doing attention.

---

## 6. Potential outputs

- **Best case (third paper):** *"Mechanistic interpretation of a fine-tuned single-cell foundation model for cancer stemness"* — novel because MI of single-cell foundation models is underexplored; contributes a confound-audited, causally-localized, interaction-aware account of the model's stemness representation, grounded against functional data.
- **Minimum useful:** a **confound-audit result** (M1) that strengthens or qualifies Paper 1, even if deeper MI is not pursued.
- **Not a claimed output:** "we found what causes stemness." (It can't be — see Section 0.)

---

## 7. Honest bottom line

This is a genuine, somewhat novel direction and a natural fit (we already have the model). But two honest flags:
1. It **cannot answer "what causes stemness"** — only "what the model uses." Framed otherwise it's overclaiming.
2. It is **a third project atop two papers not yet submitted**, and it's the most technically demanding part of the whole effort.

**Recommended sequencing:** submit the two papers first (imminent), then run **Phase M1 only** as a scoped, ~few-day decision gate. M1 is high-value regardless (it audits the existing papers) and cheap. Only commit to M2–M3 (the real research program) if M1 shows the model has learned a distinct, non-confounded, interpretable stemness representation.

**First concrete step when ready:** retrieve checkpoint-3676, extract hidden states + integrated-gradients attributions on the CSC-high/low cells, and run the layer-wise probes for stemness vs. cell-cycle/subtype.
