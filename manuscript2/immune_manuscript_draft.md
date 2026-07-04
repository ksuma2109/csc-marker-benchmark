<!--
SECOND MANUSCRIPT (separate from the CSC marker-benchmark paper).
Working draft. Analyses: notebooks/phase_I1-I5. Findings: results/B_immune_microenvironment_summary.md.
Target venue candidates: Cancer Immunology Research; J. for ImmunoTherapy of Cancer (JITC);
npj Breast Cancer; Breast Cancer Research. Computational-only; correlational.
STATUS: draft scaffold — needs independent-cohort validation before submission (see Limitations).
-->

# Cancer stem cells in triple-negative breast cancer retain antigen presentation and engage innate, not adaptive, immune evasion

**Author:** Suma Kasa (Independent Researcher) · ksuma2109@gmail.com

**Running title:** Subtype-specific innate immune evasion by breast cancer stem cells

## Abstract

**Background.** Cancer stem cells (CSCs) drive recurrence and are widely presumed to evade immunity, but the mechanism is unclear and often assumed to be loss of MHC-I antigen presentation, an immune-cold niche, or PD-1/PD-L1 engagement. We tested these assumptions directly across transcriptomic, spatial, and protein modalities in human breast cancer.

**Methods.** Using the Wu et al. breast atlas (GSE176078; 100,064 cells, 26 patients), Visium spatial transcriptomics (6 sections), and CITE-seq surface protein, we scored malignant-cell stemness and related it to antigen-presentation machinery, an expanded panel of tumour-intrinsic immune-evasion ligands, T-cell infiltration/exhaustion, and spatial T-cell co-localization. All cell-level associations were re-tested within molecular subtype.

**Results.** CSC state did not fit any classic evasion model. CSC-high malignant cells **retained** (even elevated) classical MHC-I antigen presentation; tumour stemness was uncorrelated with T-cell infiltration; and spatial analysis showed no consistent CSC–T-cell exclusion (mean per-spot r≈+0.06 across sections). PD-L1 was essentially unchanged in CSCs. Instead, a per-ligand analysis nominated **innate-immune** checkpoints — **CD47** (anti-macrophage) and **HLA-E** (anti-NK) — as the CSC-associated axis. Subtype stratification showed this phenotype is **specific to triple-negative (and, weakly, HER2+) breast cancer**: within TNBC, CSC-high cells upregulated CD47 (r=+0.23) and HLA-E (r=+0.27) and retained MHC-I (r=+0.35), whereas in luminal (ER+) tumours MHC-I showed no association and CD47 reversed. An RNA exhaustion signature was cross-validated at the protein level (CITE-seq). Critically, the two core findings — **CD47 upregulation and MHC-I retention** in TNBC CSCs — **replicated in two independent TNBC cohorts** (Pal et al. 2021, GSE161529: CD47 4/4, MHC-I 3/4 tumours; Gao et al. 2020, GSE148673: CD47 4/5, MHC-I 3/5), while PD-L1 remained near zero in both. HLA-E replicated inconsistently across cohorts and is treated as secondary.

**Conclusions.** Breast CSCs are not defeated by antigen-presentation loss and are not immune-cold; rather, TNBC CSCs upregulate the innate-immune checkpoint **CD47** while remaining T-cell-visible (MHC-I retained) — a replicated, subtype-specific phenotype. This nominates **anti-CD47 (magrolimab)**, potentially with concurrent T-cell engagement, over anti-PD-1 for CSC-directed immunotherapy in triple-negative breast cancer.

**Keywords:** cancer stem cells, triple-negative breast cancer, immune evasion, CD47, HLA-E, antigen presentation, single-cell RNA-seq, spatial transcriptomics

## 1. Introduction
*(To draft: CSCs drive relapse and therapy resistance; immunotherapy often spares CSCs; the presumed mechanisms (MHC-I loss, immune-cold, PD-1) are largely untested at single-cell resolution in CSCs; innate checkpoints (CD47, HLA-E) as an emerging axis; subtype heterogeneity of breast cancer immunobiology. Gap: a direct, multi-modal, subtype-resolved test of how CSCs relate to immunity.)*

## 2. Results
### 2.1 Breast CSCs do not fit classic immune-evasion models (I1, I2)
CSC-high malignant cells *retain* (even elevate) MHC-I antigen presentation; tumour stemness is uncorrelated with T-cell infiltration; and Visium spatial analysis shows no consistent CSC–T-cell exclusion (mean per-spot r≈+0.06). CSCs are neither antigen-presentation-null nor immune-cold.

### 2.2 The CSC-associated evasion ligand is CD47, not PD-L1 (I4)
Across an expanded panel of tumour-intrinsic evasion ligands, CSC-high cells upregulate the innate-immune "don't-eat-me" ligand **CD47** (and, more weakly, HLA-E), whereas PD-L1/PD-L2 are essentially unchanged. CD24 — itself a phagocytosis-checkpoint ligand — is *down* in CSCs (as required by the CD44⁺CD24⁻ definition), a sanity check on the labelling.

### 2.3 The phenotype is specific to triple-negative breast cancer (I5, Figure 1A)
Subtype stratification shows CD47 upregulation and MHC-I retention hold within TNBC (and, weakly, HER2⁺) but vanish or reverse in luminal (ER⁺) tumours — a TNBC-specific CSC phenotype.

### 2.4 CD47 upregulation and MHC-I retention replicate in two independent cohorts (I7, I9, Figure 1B) — KEY RESULT
In two independent TNBC scRNA-seq cohorts (Pal et al. 2021; Gao et al. 2020, malignant cells identified by copyKAT), CSC CD47 upregulation replicated (4/4 and 4/5 tumours) as did MHC-I retention (3/4 and 3/5), while PD-L1 stayed near zero (mean r=+0.07 and −0.01). HLA-E replicated inconsistently (2/4, 4/5) and is secondary.

### 2.5 Protein-level validation of the exhaustion signature (I3)
The RNA T-cell-exhaustion signature is confirmed against CITE-seq surface protein (PD-1/CTLA-4/TIGIT/LAG-3/TIM-3) in the two well-powered samples (r=+0.20, +0.41).

*(Numbers/figures in results/ — I1–I9; main figure manuscript2/figures/Figure1.png.)*

## 3. Discussion
*(To draft: reframes CSC immune evasion from adaptive/antigen-loss to innate/subtype-specific; therapeutic implication favouring anti-CD47/anti-NKG2A + the retained-MHC-I argument for concurrent T-cell engagement; relationship to the companion CSC surface-marker work for CAR-T; why ER+ differs.)*

## 4. Limitations
- Single cohort; **independent-cohort validation is required** (an external TNBC scRNA-seq dataset) before firm conclusions.
- Correlational; module scores depend on gene-set choice; HLA-E is itself a CSC-associated gene (its evasion interpretation is partly confounded).
- Patient-level (n=20) and spatial (n=6 sections) are underpowered; Visium spots are ~10-cell mixtures (regional, not single-cell contact).
- Therapeutic nominations are computational hypotheses, not experimental results.

## Data & Code Availability
scRNA-seq GSE176078; Visium Zenodo 4739739; CITE-seq SCP1039. Analysis code: `notebooks/phase_I1–I5.py` (repository shared with the companion study).

---
## Development checklist before submission
- [x] Independent-cohort validation #1 (Pal 2021 / GSE161529): CD47 4/4, MHC-I 3/4, HLA-E 2/4 (I7).
- [x] Independent-cohort validation #2 (Gao 2020 / GSE148673, copyKAT): CD47 4/5, MHC-I 3/5, PD-L1 ~0 (I9).
- [x] Lead framing on **CD47** (replicated in 2 cohorts); HLA-E demoted; PD-L1 softened.
- [x] Main figure (Figure 1): subtype specificity + two-cohort replication.
- [x] Macrophage-niche analysis (I6): NOT robustly supported — reported as a limitation.
- [ ] Full Introduction, Discussion, Methods prose; references; target-journal formatting.
- [ ] Optional: functional phagocytosis assay for CD47 (wet lab); a third validation cohort (Bassez 2021 BIOKEY).
