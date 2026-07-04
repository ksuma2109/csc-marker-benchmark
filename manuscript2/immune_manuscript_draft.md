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

**Results.** CSC state did not fit any classic evasion model. CSC-high malignant cells **retained** (even elevated) classical MHC-I antigen presentation; tumour stemness was uncorrelated with T-cell infiltration; and spatial analysis showed no consistent CSC–T-cell exclusion (mean per-spot r≈+0.06 across sections). PD-L1 was essentially unchanged in CSCs. Instead, a per-ligand analysis nominated **innate-immune** checkpoints — **CD47** (anti-macrophage) and **HLA-E** (anti-NK) — as the CSC-associated axis. Subtype stratification showed this phenotype is **specific to triple-negative (and, weakly, HER2+) breast cancer**: within TNBC, CSC-high cells upregulated CD47 (r=+0.23) and HLA-E (r=+0.27) and retained MHC-I (r=+0.35), whereas in luminal (ER+) tumours MHC-I showed no association and CD47 reversed. An RNA exhaustion signature was cross-validated at the protein level (CITE-seq).

**Conclusions.** Breast CSCs are not defeated by antigen-presentation loss and are not immune-cold; rather, TNBC CSCs engage innate-immune evasion (CD47/HLA-E) while remaining T-cell-visible. This nominates innate-immune combination partners — anti-CD47 (magrolimab) and anti-NKG2A (monalizumab) — over anti-PD-1 for CSC-directed immunotherapy in triple-negative breast cancer.

**Keywords:** cancer stem cells, triple-negative breast cancer, immune evasion, CD47, HLA-E, antigen presentation, single-cell RNA-seq, spatial transcriptomics

## 1. Introduction
*(To draft: CSCs drive relapse and therapy resistance; immunotherapy often spares CSCs; the presumed mechanisms (MHC-I loss, immune-cold, PD-1) are largely untested at single-cell resolution in CSCs; innate checkpoints (CD47, HLA-E) as an emerging axis; subtype heterogeneity of breast cancer immunobiology. Gap: a direct, multi-modal, subtype-resolved test of how CSCs relate to immunity.)*

## 2. Results
### 2.1 CSC-high malignant cells retain antigen presentation (I1)
### 2.2 Tumour stemness is uncorrelated with T-cell infiltration; no spatial exclusion (I1, I2)
### 2.3 The CSC-associated evasion axis is innate (CD47/HLA-E), not PD-1 (I4)
### 2.4 The phenotype is specific to triple-negative breast cancer (I5)
### 2.5 Protein-level validation of the exhaustion signature (I3)

*(Numbers and figures are in results/ — I1_immune_evasion, I2_spatial_exclusion, I3_protein_exhaustion, I4_checkpoint_axis, I5_subtype_stratified.)*

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
- [ ] Independent-cohort validation of the TNBC CSC → CD47/HLA-E finding (external TNBC scRNA-seq).
- [ ] Macrophage-niche analysis: are CSC-high TNBC regions macrophage-rich / phagocytosis-suppressed (spatial CD47–SIRPA / macrophage co-localization)?
- [ ] Explicit subtype-matched figures (main figure = TNBC).
- [ ] References; full Methods; format for target journal.
