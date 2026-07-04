---
name: project-csc-research
description: Bioinformatics project — Cancer Stem Cell Identification and Characterization using scRNA-seq, trajectory analysis, and cross-cancer marker discovery
metadata:
  type: project
---

# Project: Cancer Stem Cell (CSC) Identification & Characterization

---

## Goal

Build a computational pipeline to **identify, characterize, and compare Cancer Stem Cells (CSCs)** across cancer types using publicly available single-cell RNA-seq data. The project will:

1. Identify CSC subpopulations within tumors using transcriptomic and stemness signatures
2. Characterize the molecular programs (gene expression, signaling pathways, trajectory) that define CSC identity
3. Discover shared and cancer-type-specific CSC surface markers across multiple cancers
4. Understand how CSCs relate to normal stem cell programs through cross-dataset comparison

**End deliverable:** A reproducible, well-documented bioinformatics workflow that can be applied to new tumor scRNA-seq datasets to identify and profile CSC subpopulations — suitable as a portfolio project, research manuscript, or lab tool.

---

## Why This Matters

- CSCs drive tumor recurrence, metastasis, and therapy resistance
- They comprise only 0.01–2% of tumor cells — invisible to bulk sequencing, detectable only by scRNA-seq
- Current CSC markers (CD44, CD133, ALDH1A1) are not universal — only 32% of ovarian cancer patients express ALDH1A1 (PMC9916537)
- CSCs are now understood as **dynamic states**, not fixed populations (Huang et al. 2025) — requiring trajectory-based computational approaches
- A cross-cancer CSC atlas does not yet exist at the single-cell level at scale

---

## Execution Strategy — Two Stages

The project runs in two parallel stages on the same breast cancer dataset (GSE176078), then compares results:

```
GSE176078 (100,064 cells — Wu et al. 2021 Breast Cancer)
         │
         ├──► STAGE 1: Standard scRNA-seq pipeline (Baseline)
         │    Clustering → Marker genes via differential expression
         │    Output: CSC marker list (literature-anchored)
         │
         └──► STAGE 2: Geneformer transformer (Novel)
              Attention weights → Data-driven gene program
              Output: CSC marker list (co-expression-anchored)
                                │
                                ▼
                    COMPARISON: Overlap + Novel discoveries
                    - Genes in both  → high-confidence CSC markers
                    - Stage 2 only   → new candidates missed by clustering
                    - Stage 1 only   → literature markers not captured by model
```

**Why two stages?**
Stage 1 (clustering) identifies CSC markers by comparing cells against each other using differential expression — it finds genes that are statistically enriched in CSC clusters relative to other cells. It is anchored to known biology through the stemness scoring step.

Stage 2 (Geneformer) identifies CSC markers by learning which gene co-expression patterns distinguish CSC states — without using any prior gene list. Attention weights reveal which genes the model relies on to make the CSC/non-CSC distinction.

Comparing both exposes what clustering misses and what the transformer cannot find without context.

---

## Project Scope

Three interconnected sub-projects, ordered by complexity:

| Sub-project | Focus | Difficulty |
|---|---|---|
| A | scRNA-seq trajectory analysis of CSC differentiation states (both stages) | Beginner–Intermediate |
| B | CSC surface marker discovery across cancer types | Intermediate |
| C | Cross-cancer CSC transcriptomic atlas | Intermediate–Advanced |

Sub-project A uses both stages. Sub-projects B and C extend the winning or combined approach.

---

## Implementation Plan

### Sub-project A: scRNA-seq Trajectory Analysis of CSC States

**Goal:** For a single cancer type, identify CSC subpopulations and map how they relate to more differentiated tumor cells using pseudotime and RNA velocity.

---

#### Phase A1 — Data Acquisition (Week 1–2)

**Task:** Download a well-characterized scRNA-seq tumor dataset.

**Recommended starting datasets:**
- Breast cancer: GSE176078 (Wu et al. 2021, Nature Genetics — 100,064 cells)
- HNSCC: GSE234933 (Zhang et al. 2024, Nature Communications — 120,952 cells)
- Glioblastoma: GSE84465 or GSE131928

**Steps:**
1. Create a GEO account at ncbi.nlm.nih.gov/geo
2. Search for your chosen dataset by GSE accession number
3. Download the count matrix (usually a `.h5`, `.h5ad`, or `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` set)
4. Download the associated metadata (cell barcodes, sample info, clinical annotations)
5. Organize your project directory:
```
project/
├── data/
│   ├── raw/          ← downloaded files go here
│   └── processed/    ← outputs from your pipeline
├── scripts/          ← R or Python scripts
├── results/
│   ├── figures/
│   └── tables/
└── README.md
```

**Tools needed:** `wget` or `curl` (command line), GEO browser

---

#### Phase A2 — Quality Control & Preprocessing (Week 2–3)

**Goal:** Remove low-quality cells and normalize data so cells are comparable.

**Steps:**
1. Load the count matrix into Seurat (R) or scanpy (Python)
2. Calculate per-cell QC metrics:
   - Number of genes detected per cell
   - Total UMI count per cell
   - Percentage of mitochondrial gene reads (high % = dying cell)
3. Plot violin plots of all three metrics — inspect distributions before setting thresholds
4. Filter cells (thresholds are dataset-specific — do NOT use fixed numbers):
   - Remove cells with very few genes (empty droplets)
   - Remove cells with very high gene counts (doublets — two cells as one)
   - Remove cells with high mitochondrial % (dying cells, typically >10–20% depending on tissue)
5. Run DoubletFinder or scDblFinder to detect and remove doublets (~5% of barcodes)
6. Normalize counts (SCTransform in Seurat, or log1p normalization in scanpy)

**Key output:** A clean, normalized count matrix ready for analysis

**Common mistake to avoid:** Do not use fixed QC thresholds copied from a paper. Always plot your data first.

---

#### Phase A3 — Dimensionality Reduction & Clustering (Week 3–4)

**Goal:** Compress 20,000-gene space into a 2D representation and group similar cells.

**Steps:**
1. Find highly variable genes (top 2,000–3,000 genes that vary most across cells)
2. Run PCA (Principal Component Analysis) — reduce to top 30–50 PCs
3. Correct for batch effects if data comes from multiple patients/samples:
   - Use **Harmony** (batch correction tool — NOT pseudotime)
   - Or **scVI** for deep learning-based integration
4. Build a nearest-neighbor graph on the corrected embeddings
5. Run UMAP for visualization (2D layout)
6. Run Leiden or Louvain clustering to group cells
7. Experiment with different resolutions (0.1–1.5) — lower = fewer, bigger clusters

**Key output:** UMAP plot with color-coded clusters

---

#### Phase A4 — Cell Type Annotation (Week 4–5)

**Goal:** Assign a biological identity to each cluster (cancer cell, T cell, fibroblast, endothelial, CSC, etc.).

**Steps:**
1. Find marker genes for each cluster:
   - Seurat: `FindAllMarkers()` — compares each cluster vs. all others
   - Look for genes with high log2FC and low adjusted p-value
2. Annotate clusters using known markers:
   - Cancer epithelial: EPCAM, KRT8, KRT18
   - T cells: CD3D, CD3E, CD8A, CD4
   - B cells: MS4A1 (CD20), CD79A
   - Macrophages: CD68, CSF1R
   - Fibroblasts: VIM, FAP, ACTA2
   - CSC-like: CD44, ALDH1A1, SOX2, OCT4, NANOG, CD24-low
3. Use automated annotation as a starting point:
   - `SingleR` (R) or `CellTypist` (Python) — reference-based annotation
   - Cross-check with manual marker inspection
4. Subset the epithelial/malignant cells for CSC analysis

**Key output:** Annotated UMAP — each cluster labeled with cell type

---

#### Phase A5 — CSC Identification (Week 5–6)

**Goal:** Within the malignant cell compartment, identify which cells have stem-like transcriptional programs.

**Steps:**
1. **Score stemness using published gene signatures:**
   - Use AUCell or UCell in R to score each cell against:
     - Malta et al. 2018 mRNAsi gene list (pluripotency-derived stemness signature)
     - CancerSEA stemness gene set (14-state functional atlas)
     - Cancer-type-specific CSC signature from literature
   - Cells with high scores are CSC-like
2. **Identify CSC marker expression:**
   - Plot expression of CD44, ALDH1A1, SOX2, NANOG, OCT4 on UMAP
   - Use FeaturePlot, VlnPlot, DotPlot
3. **Compare stemness-high vs. stemness-low clusters:**
   - Differential expression between high-score and low-score cells
   - Pathway enrichment: run fgsea or GSEA on the DE results
   - Look for enrichment in: Wnt, Notch, Hedgehog, TGF-β gene sets (from MSigDB)
4. **CytoTRACE2** (optional, newer tool):
   - Deep learning framework that predicts developmental potential from scRNA-seq
   - Ranks each cell from most stem-like to most differentiated
   - Apply to malignant cells to identify the most dedifferentiated (CSC-like) subpopulation

**Key output:** CSC subpopulation identified, scored, and labeled on UMAP

---

#### Phase A6 — Trajectory & Pseudotime Analysis (Week 6–8)

**Goal:** Map the path from CSC state → differentiated cancer cell, and find which genes change along that path.

**Steps:**
1. **Monocle 3 pseudotime:**
   - Input: the malignant cell subset from Phase A4
   - Monocle 3 learns the trajectory structure directly from data (principal graph)
   - Set the root cells = your CSC cluster (the starting point)
   - Output: each cell gets a pseudotime value (0 = most stem-like, high = most differentiated)
   - Plot pseudotime on UMAP
2. **Find genes that change along pseudotime:**
   - Monocle 3: `graph_test()` — identifies genes with expression that varies significantly across the trajectory
   - Plot top genes as heatmaps across pseudotime
   - These are candidate genes driving CSC → differentiation transition
3. **RNA velocity with scVelo:**
   - Requires re-processing raw BAM files to get spliced/unspliced counts (use `velocyto` CLI or `STARsolo`)
   - Run scVelo dynamic model (more accurate than steady-state velocyto)
   - Overlay velocity arrows on UMAP — arrows show where each cell is heading
   - This confirms the direction of the trajectory found by Monocle
4. **Validate trajectory against known biology:**
   - CSC markers should be high at root (pseudotime 0)
   - Differentiation markers should increase at later pseudotime
   - EMT markers (VIM, SNAI1, TWIST1) often peak in the middle

**Key output:** Pseudotime trajectory + RNA velocity plot + list of trajectory-associated genes

---

### Stage 2: Geneformer — Attention-Based CSC Gene Program Discovery

**Goal:** Use the pre-trained Geneformer transformer model (trained on 30M single cells) to discover which gene co-expression patterns define CSC identity — without using any prior gene list. Compare results against Stage 1.

**Dataset:** Same GSE176078 breast cancer data. Focus on Cancer Epithelial cells (24,489 cells).

**Why Geneformer works here:** Geneformer represents each cell as a ranked sequence of its expressed genes (highest → lowest). Its transformer attention mechanism learns which genes matter *in the context of which other genes are co-expressed* — capturing the dynamic, co-occurrence-based nature of CSC states that static gene scoring cannot.

---

#### Phase G1 — Generate CSC Pseudo-Labels (builds on A5)

**Goal:** Create CSC-high / CSC-low labels for Cancer Epithelial cells to supervise fine-tuning.

**Steps:**
1. Load `data/processed/brca_A2_preprocessed.h5ad`
2. Subset to Cancer Epithelial cells only (`celltype_major == "Cancer Epithelial"`)
3. Score each cell for stemness using `decoupler` with the Malta et al. mRNAsi gene signature
4. Rank cells by stemness score
5. Label top 25% as `csc_high`, bottom 25% as `csc_low`, discard middle 50%
   - This creates a clean binary contrast for fine-tuning
6. Save the labelled subset

**Key output:** Labelled AnnData with 6,122 cells (Cancer Epithelial, top/bottom quartile)

---

#### Phase G2 — Tokenize for Geneformer

**Goal:** Convert scRNA-seq counts into Geneformer's input format.

**How Geneformer tokenizes:**
- For each cell, genes are sorted from highest → lowest expression
- Each gene is looked up in Geneformer's Ensembl token dictionary
- The result is a ranked list of token IDs (like a sentence where words are genes ordered by activity)

**Steps:**
1. Convert gene symbols → Ensembl IDs using `mygene`
   - e.g., `TP53` → `ENSG00000141510`
   - Genes without Ensembl IDs are dropped
2. Add Ensembl IDs as a column in `adata.var`
3. Write `.loom` file (Geneformer's preferred input format)
4. Run `TranscriptomeTokenizer` on the loom file
   - Output: HuggingFace `Dataset` with one tokenized cell per row
5. Attach CSC labels to the tokenized dataset

**Key output:** Tokenized HuggingFace dataset ready for fine-tuning

---

#### Phase G3 — Fine-Tune Geneformer

**Goal:** Train Geneformer to distinguish CSC-high from CSC-low Cancer Epithelial cells.

**What happens during fine-tuning:**
- Geneformer's pre-trained weights (learned from 30M cells) already encode general gene relationships
- Fine-tuning adjusts these weights specifically for the CSC-high vs CSC-low distinction
- The attention mechanism learns to up-weight genes that are most predictive of CSC identity

**Steps:**
1. Load pre-trained Geneformer model (6-layer BERT architecture)
2. Add a classification head (2 classes: csc_high, csc_low)
3. Fine-tune for 5–10 epochs with early stopping
4. Evaluate on held-out 20% test set
5. Save the fine-tuned model

**Hardware note:** Runs on Apple M-series MPS (Metal) or CPU. Expect ~30–60 min on CPU for 6,000 cells.

**Key output:** Fine-tuned Geneformer model with classification accuracy + F1 score

---

#### Phase G4 — Extract Attention Weights → Gene Importance

**Goal:** After fine-tuning, ask the model: "which genes did you pay most attention to when classifying a cell as CSC?"

**How attention extraction works:**
- Each transformer layer has multiple attention heads
- Each head computes attention scores between every pair of genes in the sequence
- A gene with high attention is one the model considers contextually important
- We aggregate attention across all heads and layers to get a per-gene importance score

**Steps:**
1. Forward-pass all CSC-high cells through the fine-tuned model with `output_attentions=True`
2. For each cell, extract attention matrix `[layers × heads × genes × genes]`
3. Average across layers and heads → `[genes × genes]` attention matrix per cell
4. Sum each gene's column → scalar importance score per gene per cell
5. Average across all CSC-high cells → final gene importance ranking
6. Map token IDs back to gene symbols
7. Save ranked gene list as `results/tables/geneformer_csc_gene_ranking.csv`

**Key output:** Ranked list of genes by attention importance — the data-driven CSC gene program

---

#### Phase G5 — Comparison: Stage 1 vs Stage 2

**Goal:** Systematically compare CSC markers found by clustering (Stage 1) vs attention (Stage 2).

**Steps:**
1. Load Stage 1 marker genes (from Phase A5 differential expression)
2. Load Stage 2 attention-ranked genes (top 200 by attention score)
3. Compute overlap (Jaccard index, hypergeometric test for significance)
4. Categorize each gene:
   - **Shared** — appears in both lists → high-confidence CSC marker
   - **Stage 2 only** → novel candidates missed by clustering
   - **Stage 1 only** → classical markers the model de-emphasizes
5. For Stage-2-only genes: check against CSC literature (PubMed, CSCTT database)
   - Known but missed by clustering → validates attention approach
   - Unknown → genuinely novel candidates for experimental follow-up
6. Visualize:
   - Venn diagram of overlap
   - Scatter plot: attention score (x) vs log2FC from DE (y)
   - Heatmap of top 50 shared genes across CSC-high cells

**Key output:** Comparison table, figures, list of novel CSC candidates

---

### Sub-project B: CSC Marker Discovery Across Cancer Types

**Goal:** Apply the Phase A5 stemness scoring to multiple cancer types and identify which markers are shared vs. cancer-specific.

**Steps:**
1. Download 3–5 scRNA-seq datasets from different cancer types (GEO or CellxGene)
2. Run the Phase A1–A5 pipeline independently for each cancer type
3. For each cancer, produce:
   - A ranked list of CSC marker genes (from FindMarkers on CSC vs. non-CSC clusters)
   - Stemness score distribution across clusters
4. Compare marker lists across cancers:
   - Find the intersection (universal CSC markers)
   - Find cancer-specific markers (present in only one or two types)
5. Validate candidates against:
   - CancerSCEM 2.0 (1,466 datasets, 74 cancer types)
   - CSCTT database (135 experimentally validated CSC targets)
6. Prioritize surface proteins (for potential therapeutic targeting) using Uniprot/Human Protein Atlas to filter for membrane-localized proteins

**Key output:** Table of pan-cancer CSC markers + cancer-specific markers, with evidence scores

---

### Sub-project C: Cross-Cancer CSC Atlas

**Goal:** Integrate datasets from multiple cancer types into a unified atlas of CSC transcriptional states.

**Steps:**
1. Collect 5–10 scRNA-seq datasets from CellxGene (pre-standardized, easier to integrate)
2. Extract malignant cells from each dataset (using copy number variation inference or known markers)
3. Integrate across datasets using:
   - **scVI** (recommended for large cross-cancer integration)
   - Or **Seurat v5 integration** pipeline
4. Cluster the integrated malignant cells
5. Score all cells for stemness (mRNAsi, AUCell)
6. Identify CSC clusters — label them as "CSC-like"
7. Characterize each CSC cluster:
   - Which cancer types contribute to it? (pan-cancer vs. type-specific)
   - What are its top marker genes?
   - What pathways are enriched?
8. Build a reference UMAP showing all cancer types together with CSC states highlighted
9. Compare your CSC clusters against CancerSEA's 14 functional states

**Key output:** Cross-cancer CSC atlas — a reference embedding showing how CSC states relate across tumor types

---

## Tools & Software Stack

### Stage 1 — Standard Pipeline (all installed in venv)
```
scanpy          — core scRNA-seq analysis
scvelo          — RNA velocity / pseudotime
decoupler       — stemness scoring (mRNAsi, AUCell-equivalent)
celltypist      — automated cell type annotation
harmonypy       — batch correction
leidenalg       — clustering
umap-learn      — dimensionality reduction
pyscenic        — gene regulatory network analysis
scikit-misc     — required for seurat_v3 HVG selection
```

### Stage 2 — Geneformer Pipeline (all installed in venv)
```
geneformer      — transformer model (installed from HuggingFace repo)
torch           — PyTorch backend
transformers    — pinned to ==4.46 (Geneformer requirement)
datasets        — HuggingFace datasets for tokenized input
peft            — parameter-efficient fine-tuning
mygene          — gene symbol → Ensembl ID conversion
tensorboard     — training monitoring
```

### Key Learning Resources
- scanpy tutorials: scanpy.readthedocs.io
- OSCA book (free): bioconductor.org/books/release/OSCA/
- scVelo docs: scvelo.readthedocs.io
- Geneformer paper: Theodoris et al. 2023, Nature
- StatQuest (YouTube) — statistics concepts explained clearly

---

## Data Sources

| Source | URL | What to Use It For |
|---|---|---|
| GEO | ncbi.nlm.nih.gov/geo | Download raw scRNA-seq datasets |
| CellxGene | cellxgene.cziscience.com | Pre-standardized, ready-to-use datasets |
| CancerSCEM 2.0 | (Nucleic Acids Research 2025) | 1,466 cancer scRNA-seq datasets, 74 types |
| CancerSEA | (Yuan et al. 2019) | 14 CSC functional state gene sets |
| MSigDB | gsea-msigdb.org | Pathway gene sets for enrichment analysis |
| TCGA | portal.gdc.cancer.gov | Bulk RNA-seq + clinical data for validation |
| CSCTT | (MDPI 2025) | 135 validated CSC therapeutic targets |

---

## Milestones

### Stage 1 — Standard Pipeline (Baseline)

| Milestone | Deliverable | Status |
|---|---|---|
| M1 | Raw data downloaded (GSE176078, 100,064 cells) | ✅ Done |
| M2 | QC & preprocessing complete — `brca_A2_preprocessed.h5ad` saved | ✅ Done |
| M3 | Dimensionality reduction + clustering (UMAP, Leiden) | Phase A3 |
| M4 | Cell type annotation — Cancer Epithelial cells isolated | Phase A4 |
| M5 | CSC subpopulation identified via mRNAsi stemness scoring | Phase A5 |
| M6 | Differential expression → Stage 1 CSC marker gene list | Phase A5 |
| M7 | Pseudotime trajectory + RNA velocity (Monocle / scVelo) | Phase A6 |

### Stage 2 — Geneformer Pipeline

| Milestone | Deliverable | Status |
|---|---|---|
| G1 | CSC pseudo-labels generated (top/bottom quartile stemness) | Phase G1 |
| G2 | Data tokenized in Geneformer format (gene symbol → Ensembl → tokens) | Phase G2 |
| G3 | Geneformer fine-tuned on CSC-high vs CSC-low | Phase G3 |
| G4 | Attention weights extracted → Stage 2 gene importance ranking | Phase G4 |
| G5 | Stage 1 vs Stage 2 comparison complete — novel candidates identified | Phase G5 |

### Sub-projects B & C

| Milestone | Deliverable | Target |
|---|---|---|
| B1 | CSC marker comparison across 3+ cancer types (best method from A) | Week 14 |
| C1 | Cross-cancer CSC atlas (integrated, clustered, annotated) | Week 20 |

---

## Key Verified Facts (Research-Backed)

- CSCs are 0.01–2% of total tumor cells — too rare for bulk sequencing
- ALDH1A1 is only expressed in 32% of ovarian cancer patients — markers are not universal (PMC9916537, 3-vote verified)
- ALDH1A1+ ovarian cancer patients: median OS 34 months vs 58 months for ALDH1A1- (p=0.00022)
- Notch + gamma-secretase inhibitors reduce GBM CSC markers: CD133, NESTIN, BMI1, OLIG2 (PMC11856057, 3-vote verified)
- ~5% of cell barcodes in scRNA-seq are doublets — always run doublet removal (PMC10189648, 3-vote verified)
- scVelo uses likelihood-based dynamic modeling, improving on velocyto's steady-state assumption (2-vote verified)
- **Harmony is batch correction, NOT a pseudotime tool** — a common misconception (3-vote refuted)
- QC thresholds are not universal — 86% of papers use 5–10% mitochondrial cutoff, not 20% (refuted claim)
- CSCs are dynamic states cells move in and out of — not fixed marker-defined populations (Huang et al. 2025)

---

## Related Files

- [CSC_Bioinformatics_Beginner_Guide.md](CSC_Bioinformatics_Beginner_Guide.md) — full background knowledge guide
- [project_csc_background_knowledge.md](project_csc_background_knowledge.md) — condensed knowledge map
- [MEMORY.md](MEMORY.md) — project memory index
