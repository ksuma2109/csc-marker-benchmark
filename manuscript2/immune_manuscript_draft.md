<!--
SECOND MANUSCRIPT — CSC innate immune evasion (CD47) in TNBC.
Target venue candidates: Cancer Immunology Research; J. ImmunoTherapy of Cancer (JITC);
npj Breast Cancer; Breast Cancer Research. Computational; two-cohort replicated.
NOTE: reference details (volume/pages/DOIs) must be verified in a reference manager before submission.
Figures: manuscript2/figures/. Analyses: notebooks/phase_I1-I9.
-->

# Triple-negative breast cancer stem cells retain antigen presentation and upregulate the innate-immune checkpoint CD47

**Author:** Suma Kasa¹

**Affiliation:** ¹Independent Researcher

**ORCID:** 0009-0001-3253-6972

**Correspondence:** ksuma2109@gmail.com

**Running title:** CD47 innate immune evasion by TNBC cancer stem cells

---

## Abstract

**Background.** Cancer stem cells (CSCs) drive recurrence and are widely presumed to evade immunity, but the mechanism is unclear and often assumed to be loss of MHC-I antigen presentation, an immune-cold niche, or PD-1/PD-L1 engagement. We tested these assumptions directly across transcriptomic, spatial, and protein modalities in human breast cancer, and sought independent replication.

**Methods.** Using the Wu et al. breast atlas (GSE176078; 100,064 cells, 26 patients) with Visium spatial transcriptomics (6 sections) and CITE-seq surface protein, we scored malignant-cell stemness and related it to antigen-presentation machinery, an expanded panel of tumour-intrinsic immune-evasion ligands, T-cell infiltration and exhaustion, and spatial T-cell co-localization. All cell-level associations were re-tested within molecular subtype and replicated in two independent TNBC scRNA-seq cohorts (Pal et al. 2021, GSE161529; Gao et al. 2020, GSE148673, malignant cells identified by copyKAT).

**Results.** CSC state fit no classic evasion model: CSC-high malignant cells *retained* MHC-I antigen presentation, tumour stemness was uncorrelated with T-cell infiltration, and spatial analysis showed no consistent CSC–T-cell exclusion. Among evasion ligands, CSC-high cells upregulated the innate-immune "don't-eat-me" ligand **CD47** whereas PD-L1 was unchanged; CD24 (a phagocytosis checkpoint) was concordantly *down*, as required by the CD44⁺CD24⁻ definition. Subtype stratification localized the phenotype to **triple-negative (and, weakly, HER2⁺) breast cancer**, with no effect or reversal in luminal tumours. The two core findings — CD47 upregulation and MHC-I retention — **replicated in both independent TNBC cohorts** (CD47 4/4 and 4/5 tumours; MHC-I 3/4 and 3/5), while PD-L1 remained near zero.

**Conclusions.** TNBC CSCs are not defeated by antigen-presentation loss and are not immune-cold; instead they upregulate the innate-immune checkpoint CD47 while remaining T-cell-visible. This replicated, subtype-specific finding nominates anti-CD47 (magrolimab), potentially with concurrent T-cell engagement, over anti-PD-1 for CSC-directed immunotherapy in triple-negative breast cancer.

**Keywords:** cancer stem cells, triple-negative breast cancer, immune evasion, CD47, antigen presentation, single-cell RNA-seq, spatial transcriptomics, tumour immunology

---

## 1. Introduction

Cancer stem cells (CSCs) are a minority tumour subpopulation with self-renewal and tumour-initiating capacity that is thought to underlie relapse, metastasis, and therapy resistance (Reya et al., 2001; Batlle and Clevers, 2017). Because CSCs frequently persist after cytotoxic and targeted therapy, there is strong interest in whether they can be eliminated immunologically. A common assumption is that CSCs are intrinsically immune-evasive — variously attributed to downregulation of major histocompatibility complex class I (MHC-I) antigen presentation, formation of an immune-excluded ("cold") niche, or engagement of the PD-1/PD-L1 checkpoint. These mechanisms, however, have rarely been tested directly in CSCs at single-cell resolution, and breast cancer is strongly heterogeneous in immunobiology across molecular subtypes, with triple-negative breast cancer (TNBC) being both the most stem-like and the most immune-infiltrated subtype.

An orthogonal and increasingly important axis of tumour immune evasion is the innate-immune "don't-eat-me" checkpoint CD47, which engages SIRPα on macrophages to suppress phagocytosis and is upregulated across many cancers (Jaiswal et al., 2009; Willingham et al., 2012); anti-CD47 agents such as magrolimab have entered clinical testing (Advani et al., 2018). A second, related phagocytosis checkpoint is CD24–Siglec-10 (Barkal et al., 2019) — notable here because the canonical breast-CSC definition is CD44⁺CD24⁻/low, i.e. CSCs are defined by *low* CD24. Whether CSCs rely on innate checkpoints such as CD47 rather than the adaptive PD-1 axis has not been systematically examined.

Here we test how CSC state relates to tumour immunity in breast cancer across three data modalities — single-cell transcriptomics, spatial transcriptomics, and CITE-seq surface protein — and, critically, seek independent replication. We find that CSCs do not fit any classic evasion model but instead, specifically in TNBC, retain antigen presentation while upregulating CD47, a phenotype that replicates in two independent cohorts and nominates a distinct therapeutic combination.

---

## 2. Materials and Methods

### 2.1 Datasets
Discovery used the human breast cancer atlas of Wu et al. (2021) (GEO GSE176078; 100,064 cells, 26 patients; 24,489 cancer-epithelial and 35,214 T-cells), together with its Visium spatial transcriptomics (6 sections; Zenodo 4739739) and CITE-seq surface-protein (ADT) data (Single Cell Portal SCP1039). Independent validation used two TNBC scRNA-seq cohorts: Pal et al. (2021) (GSE161529; 4 treatment-naive sporadic TNBC tumours) and Gao et al. (2021) (GSE148673; 5 TNBC tumours), in which malignant cells were identified by the authors' copyKAT aneuploidy calls.

### 2.2 Stemness and immune scoring
Single-cell data were processed in Scanpy (Wolf et al., 2018): library-size normalization to 10⁴ and log1p transformation. A CSC/stemness score and gene-program scores (module scores relative to matched control gene sets) were computed with `score_genes`. Immune programs: classical MHC-I antigen presentation (HLA-A/B/C, TAP1/2, NLRC5; B2M and HLA-E were deliberately excluded because they are CSC-associated genes, to avoid a shared-gene confound), an expanded panel of tumour-intrinsic evasion ligands (CD274/PD-L1, PDCD1LG2/PD-L2, CD47, CD276/B7-H3, VTCN1/B7-H4, LGALS9, NT5E/CD73, ENTPD1/CD39, HLA-E, HLA-G, CD24), CD8 T-cell exhaustion (PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, TOX), and cytotoxicity (GZMB, GZMA, PRF1, IFNG, NKG7). Associations were quantified as per-cell Spearman correlations with stemness and as CSC-high versus CSC-low contrasts (top/bottom stemness quartiles; Mann–Whitney U).

### 2.3 Subtype stratification and validation
All cell-level associations were recomputed within each molecular subtype (ER+, TNBC, HER2+). For validation, malignant cells from each independent TNBC tumour were scored identically and the per-cell Spearman correlation of each immune program with stemness computed per tumour; replication was assessed by direction concordance and magnitude relative to discovery.

### 2.4 Spatial and protein analyses
For each Visium section, spots were scored for stemness, T-cell, macrophage, and CD47 content; spatial co-localization was tested by per-spot Spearman correlation and by pathologist spot annotations (cancer spots with vs. without lymphocytes). For protein validation, CITE-seq ADT exhaustion markers (PD-1, CTLA-4, TIGIT, LAG-3, TIM-3; CLR-normalized) were correlated per T-cell with the matched RNA exhaustion score.

### 2.5 Code and reproducibility
All analyses are scripted (`notebooks/phase_I1–I9.py`); tables and figures are in the repository (Data & Code Availability).

---

## 3. Results

### 3.1 Breast CSCs do not fit classic immune-evasion models
Contrary to the immune-cold / antigen-loss model, CSC-high malignant cells *retained* classical MHC-I antigen presentation (per-cell Spearman r = +0.26 with stemness; CSC-high +0.66 vs. CSC-low +0.19; excluding B2M/HLA-E from the score, so not a shared-gene artifact). At the patient level (20 patients), tumour stemness was uncorrelated with T-cell infiltration (r = −0.07, ns); weak positive trends for CD8 exhaustion (r = +0.32) and cytotoxicity (r = +0.44, p = 0.056) were consistent with the higher immunogenicity of stem-like (largely TNBC) tumours rather than immune coldness. Spatial transcriptomics (6 sections) showed no consistent CSC–T-cell exclusion: the per-spot stemness–T-cell correlation was negative in 2 sections, positive in 2, and null in 2 (mean r ≈ +0.06). CSCs are therefore neither antigen-presentation-null nor immune-excluded.

### 3.2 The CSC-associated evasion ligand is CD47, not PD-L1
Across an expanded panel of tumour-intrinsic evasion ligands, CSC-high cells most upregulated the innate-immune "don't-eat-me" ligand **CD47** (CSC-high − low +0.13) and HLA-E (flagged, as a CSC-associated gene), followed weakly by CD73 and B7-H3, whereas **PD-L1 and PD-L2 were essentially unchanged** (+0.009). As an internal control, CD24 — itself a Siglec-10 phagocytosis checkpoint — was *down* in CSC-high cells (−0.29), exactly as required by the CD44⁺CD24⁻/low CSC definition, validating the stemness labelling and indicating that CSCs shed one phagocytosis checkpoint (CD24) while gaining another (CD47).

### 3.3 The phenotype is specific to triple-negative breast cancer (Figure 1A)
Because TNBC is both the most stem-like and most immunogenic subtype, all cell-level associations were re-tested within subtype. CD47 upregulation and MHC-I retention held within TNBC (r = +0.23 and +0.35) and, more weakly, HER2⁺ (r = +0.11 and +0.24), but were absent or reversed in luminal ER⁺ tumours (CD47 r = −0.05; MHC-I r ≈ 0). The CSC innate-evasion phenotype is thus subtype-specific to TNBC rather than a general property of breast CSCs, and is not a cross-subtype artifact.

### 3.4 CD47 upregulation and MHC-I retention replicate in two independent cohorts (Figure 1B)
In two independent TNBC scRNA-seq cohorts, the two core findings replicated. In Pal et al. (2021; 4 tumours), CSC CD47 upregulation was positive in 4/4 tumours (mean r = +0.11) and MHC-I retention in 3/4 (mean r = +0.14); in Gao et al. (2021; 5 tumours, malignant cells by copyKAT), CD47 was positive in 4/5 (mean r = +0.08) and MHC-I in 3/5 (mean r = +0.08). PD-L1 remained near zero in both (mean r = +0.07 and −0.01), confirming it is not the CSC-associated axis. HLA-E replicated inconsistently (2/4 and 4/5) and is treated as a secondary, confounded signal. Thus CD47 upregulation and antigen-presentation retention in TNBC CSCs are reproducible across three cohorts and two independent malignant-cell-calling strategies.

### 3.5 Protein-level validation of the exhaustion signature
To confirm that the RNA-based T-cell readouts reflect protein, CITE-seq surface-protein exhaustion markers (PD-1, CTLA-4, TIGIT, LAG-3, TIM-3) were correlated per T-cell with the matched RNA exhaustion score; the two agreed in the two well-powered samples (Spearman r = +0.20 and +0.41), validating the signature.

### 3.6 No spatial macrophage niche detected
Given the CD47 finding, we tested whether CSC-high TNBC regions occupy macrophage-rich niches spatially; this was not robustly supported (CSC–macrophage co-localization in only 1 of 4 TNBC sections), indicating that the CD47 upregulation is a robust cell-intrinsic property whose macrophage-interaction consequence is not resolvable at Visium (~55 µm) resolution and would require functional assay.

---

## 4. Discussion

Across transcriptomic, spatial, and protein modalities, and in three cohorts, cancer stem cells in breast cancer do not conform to the textbook picture of immune evasion. They are not immune-cold, they are not spatially excluded from T-cells, they do not lose MHC-I antigen presentation, and they do not rely on the PD-1/PD-L1 axis. Instead, specifically in triple-negative breast cancer, CSC-high malignant cells upregulate the innate-immune checkpoint CD47 while retaining antigen presentation — a replicated, subtype-specific phenotype.

This reframes CSC immune evasion in TNBC from an adaptive/antigen-loss model to an innate/phagocytosis-checkpoint model, with three therapeutic implications. First, because CSCs retain MHC-I, they remain in principle visible to T-cells and T-cell-receptor-based immunotherapy; antigen-presentation restoration is not the required intervention. Second, the innate axis nominates **anti-CD47 (magrolimab)** — rather than anti-PD-1/PD-L1 — as the rational combination partner for CSC-directed immunotherapy in TNBC. Third, the concomitant loss of the CD24 phagocytosis checkpoint in CSCs suggests CSCs are poised toward, not away from, macrophage-mediated clearance if CD47 is blocked. The retained-MHC-I and gained-CD47 phenotype together argue for combining innate (anti-CD47) with adaptive (T-cell) engagement.

The finding is distinct from, and complementary to, the identification of CSC-specific surface antigens for chimeric antigen receptor (CAR) T-cell therapy, which is independent of MHC-I; the present analysis speaks to endogenous and innate immunity rather than to CAR design.

### 4.1 Limitations
The analyses are computational and correlational. Effect sizes are modest (per-cell r ≈ 0.1–0.35), as expected for single-cell correlations; patient-level (n = 20) and spatial (6 sections) analyses are underpowered, and several immune associations are trends rather than robust effects. HLA-E, although prominent, is itself a CSC-associated gene and replicated inconsistently, so its NKG2A-axis interpretation is not secure. The CD47 finding is cell-intrinsic and transcriptomic; its functional consequence — resistance to macrophage phagocytosis and its reversal by anti-CD47 — was not tested here and requires a functional phagocytosis assay. Validation, while spanning two independent cohorts and two malignant-cell-calling methods, used marker-based and copyKAT-based epithelial gating rather than matched CNV inference in every dataset.

### 4.2 Conclusion
TNBC cancer stem cells retain antigen presentation and upregulate the innate-immune checkpoint CD47 — a replicated, subtype-specific finding that nominates anti-CD47, potentially with concurrent T-cell engagement, as a rational strategy for CSC-directed immunotherapy in triple-negative breast cancer, and reframes CSC immune evasion as innate rather than adaptive.

---

## Data & Code Availability
scRNA-seq: GSE176078 (discovery), GSE161529 and GSE148673 (validation). Spatial: Zenodo 4739739. CITE-seq: Single Cell Portal SCP1039. Analysis code (`notebooks/phase_I1–I9.py`), result tables, and figures are available at https://github.com/ksuma2109/csc-marker-benchmark (private during peer review; public upon publication).

## Author Contributions
S.K. designed the study, performed all analyses, and wrote the manuscript.

## Conflict of Interest
The author declares no competing interests.

## Funding
This research received no external funding.

## Author Biography
**Suma Kasa** is an independent researcher in computational biology and bioinformatics, with interests in single-cell genomics, cancer stem cell biology, and tumour immunology.

---

## References
*(Author–date; verify volume/pages/DOIs in a reference manager before submission.)*

Advani, R., Flinn, I., Popplewell, L., Forero, A., Bartlett, N. L., Ghosh, N., et al. (2018). CD47 blockade by Hu5F9-G4 and rituximab in non-Hodgkin's lymphoma. *N. Engl. J. Med.* 379, 1711–1721.

Al-Hajj, M., Wicha, M. S., Benito-Hernandez, A., Morrison, S. J., and Clarke, M. F. (2003). Prospective identification of tumorigenic breast cancer cells. *Proc. Natl. Acad. Sci. USA* 100, 3983–3988.

Barkal, A. A., Brewer, R. E., Markovic, M., Kowarsky, M., Barkal, S. A., Zaro, B. W., et al. (2019). CD24 signalling through macrophage Siglec-10 is a target for cancer immunotherapy. *Nature* 572, 392–396.

Batlle, E., and Clevers, H. (2017). Cancer stem cells revisited. *Nat. Med.* 23, 1124–1134.

Gao, R., Bai, S., Henderson, Y. C., Lin, Y., Schalck, A., Yan, Y., et al. (2021). Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. *Nat. Biotechnol.* 39, 599–608.

Jaiswal, S., Jamieson, C. H. M., Pang, W. W., Park, C. Y., Chao, M. P., Majeti, R., et al. (2009). CD47 is upregulated on circulating hematopoietic stem cells and leukemia cells to avoid phagocytosis. *Cell* 138, 271–285.

Pal, B., Chen, Y., Vaillant, F., Capaldo, B. D., Joyce, R., Song, X., et al. (2021). A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast. *EMBO J.* 40, e107333.

Reya, T., Morrison, S. J., Clarke, M. F., and Weissman, I. L. (2001). Stem cells, cancer, and cancer stem cells. *Nature* 414, 105–111.

Willingham, S. B., Volkmer, J.-P., Gentles, A. J., Sahoo, D., Dalerba, P., Mitra, S. S., et al. (2012). The CD47-signal regulatory protein alpha (SIRPa) interaction is a therapeutic target for human solid tumors. *Proc. Natl. Acad. Sci. USA* 109, 6662–6667.

Wolf, F. A., Angerer, P., and Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* 19, 15.

Wu, S. Z., Al-Eryani, G., Roden, D. L., Junankar, S., Harvey, K., Andersson, A., et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nat. Genet.* 53, 1334–1347.
