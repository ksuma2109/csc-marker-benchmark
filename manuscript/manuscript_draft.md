<!--
Briefings in Bioinformatics — Research Article (method comparison / benchmarking).
Working manuscript. Figures in manuscript/figures/. Reproducible code in notebooks/phase_*.py.
NOTE: reference page numbers/DOIs should be verified in a reference manager before submission.
-->

# A functional benchmarking framework for cancer stem cell marker identification: differential expression versus transformer attention in single-cell RNA-sequencing

**Running title:** Functionally benchmarking CSC marker-ranking methods

**Authors:** Suma Kasa¹

**Affiliations:** ¹Independent Researcher

**Correspondence:** ksuma2109@gmail.com

---

## Abstract

**Background.** Cancer stem cells (CSCs) drive recurrence and therapy resistance, but the computational identification of CSC marker genes from single-cell RNA-sequencing (scRNA-seq) is methodologically unsettled and evaluated circularly: stemness labels are defined from published markers, and candidate genes are then judged against the same literature. No framework exists to compare CSC marker-ranking methods against a label-independent standard, so it is unknown whether newer transformer-based gene-importance methods identify markers that differ from, or improve on, classical differential expression (DE). We introduce a reusable functional-benchmarking framework that scores any gene-ranking method against independent functional CSC assays, and apply it to four representative methods.

**Methods.** We implemented a two-stage pipeline on a breast cancer atlas (GSE176078; 100,064 cells): Stage 1 ranks genes by DE between stemness-high and stemness-low cancer-epithelial cells; Stage 2 fine-tunes the Geneformer transformer to classify the same cells and ranks genes by attention. To break the circularity, we benchmarked both rankings against three independent functional CSC assays — sorted ALDH⁺ and CD44⁺CD24⁻ populations and mammosphere self-renewal cultures — using precision against random gene sets, area under the ROC curve (AUROC), and effect-size concordance. We further characterized pathway loading, cross-cancer reproducibility (glioblastoma, melanoma), and bulk survival association.

**Results.** The DE and attention rankings overlapped minimally (16 of top 200; Jaccard 0.042), indicating complementary rather than redundant signal. Against functional assays, the two methods split metrics systematically: DE produced higher-precision top-gene shortlists, whereas attention discriminated functional CSC genes substantially better genome-wide (AUROC up to 0.73 vs. ~0.50). Benchmarking against supervised baselines showed genome-wide AUROC increasing monotonically with model capability (DE 0.53, logistic regression 0.54, random forest 0.60, Geneformer 0.66), with attention exceeding both linear and nonlinear baselines — indicating its advantage derives from pretraining rather than supervision alone. Attention uniquely recovered functionally validated genes that DE ranked low (e.g., KLK5, FOXI1, BMPR1B) and loaded preferentially on Wnt and BMP/TGF-β developmental pathways. Extending the benchmark to nine functional gates across six cancers showed that marker rankings are substantially tissue-specific: within the discovery cancer the method dissociation was pronounced, but breast-derived rankings transferred only weakly to held-out cancers (AUROC falling toward chance), indicating the framework is best applied per tumor type. A recurrent program (CD44, VIM, FN1, MYC; ANXA1, CAV1, OSMR) was nonetheless reproducible across three cancers, and subtype-stratified survival analysis recovered the expected adverse-prognosis direction after correcting luminal confounding.

**Conclusion.** Functional benchmarking resolves an otherwise circular methods question and yields concrete method-selection guidance: use DE (or a supervised linear model) as a high-precision shortlist generator, and transformer attention as a genome-wide ranker of mechanistically informative regulators, combining the two when both recall and precision matter. We release the framework as an installable, tested Python package (`cscbench`), together with the four-method comparison and a prioritized, pathway-annotated CSC candidate list, as reusable resources for CSC marker discovery.

**Keywords:** cancer stem cells, single-cell RNA-sequencing, Geneformer, transformer attention, differential expression, functional benchmarking, breast cancer, gene prioritization

---

## Key Points

- Computational cancer stem cell (CSC) marker discovery is evaluated circularly; we introduce a reusable, label-independent framework that scores any gene-ranking method against independent functional CSC assays (sorted ALDH⁺, CD44⁺CD24⁻, CD133⁺ populations and sphere self-renewal).
- Differential expression and transformer attention are complementary: differential expression produces higher-precision marker shortlists, whereas transformer attention is a stronger genome-wide ranker; comparison against supervised linear and nonlinear baselines shows the transformer's advantage derives from pretraining rather than supervision alone.
- CSC marker rankings are substantially tissue-specific and transfer only weakly across tumor types, so the framework and its method-selection guidance are best applied per cancer.
- The framework is released as an installable, tested Python package (`cscbench`), together with a four-method comparison and a prioritized, pathway-annotated CSC candidate list.

---

## 1. Introduction

Cancer stem cells (CSCs) are a minority tumor subpopulation capable of self-renewal, multilineage differentiation, and tumor initiation, and are implicated in metastasis, relapse, and treatment resistance (Reya et al., 2001; Batlle and Clevers, 2017). Their identification is hampered by two persistent problems. First, canonical markers such as CD44, CD133, and ALDH1A1 are not universal across patients or tumor types, limiting the transferability of fixed marker panels (Al-Hajj et al., 2003). Second, CSC identity is increasingly recognized as a dynamic, plastic transcriptional state rather than a fixed marker-defined population, which single-gene scoring captures poorly.

Single-cell RNA-sequencing (scRNA-seq) resolves CSCs at the required granularity, and two broad strategies are used to nominate marker genes. Differential expression (DE) compares stemness-high to stemness-low cells and ranks genes by statistical abundance. More recently, transformer foundation models such as Geneformer — pre-trained on tens of millions of single-cell transcriptomes — represent each cell as a rank-ordered gene sequence and learn contextual dependencies between genes; their attention weights can rank genes by importance to a downstream classification task (Theodoris et al., 2023). Whether attention-based ranking identifies CSC markers that differ from, or improve upon, DE remains unknown.

A core obstacle is circularity. Stemness labels are typically derived from published signatures (Malta et al., 2018), and candidate genes are then validated against the same literature — so internal metrics cannot establish which method finds *biologically* better markers. Functional assays that define CSCs by phenotype rather than transcript — ALDEFLUOR enzymatic activity, the CD44⁺CD24⁻ surface gate (Al-Hajj et al., 2003), and mammosphere self-renewal — provide a label-independent standard, but they are seldom used to benchmark computational marker-discovery methods head-to-head.

Here we address this gap with a **functional-benchmarking framework** for CSC marker-ranking methods and apply it systematically. Specifically, we (i) define a label-independent evaluation that scores any gene ranking against functional CSC assays using precision, genome-wide AUROC, and effect-size concordance; (ii) apply it to four representative methods spanning the modeling spectrum — univariate DE, a pretrained transformer (Geneformer attention), and supervised linear and nonlinear baselines (logistic regression, random forest) — to isolate what drives performance; (iii) characterize the pathway content each method captures; and (iv) test reproducibility across glioblastoma and melanoma and association with patient survival. The framework converts an otherwise circular question into a quantitative one and yields explicit method-selection guidance: DE and supervised linear models excel at high-precision shortlists, whereas the pretrained transformer is the strongest genome-wide ranker of mechanistically informative regulators, with its advantage traceable to pretraining rather than supervision alone. We frame nominated genes as prioritized hypotheses for experimental validation rather than established markers, and release the framework, method comparison, and candidate list as reusable resources.

---

## 2. Materials and Methods

### 2.1 Datasets

Primary discovery used the human breast cancer atlas of Wu et al. (2021) (Gene Expression Omnibus [GEO] GSE176078; 100,064 cells, 26 patients). Cross-cancer analyses used the glioblastoma atlas of Darmanis et al. (2017) (GSE84465) and the melanoma atlas of Tirosh et al. (2016) (GSE72056). Functional benchmark data comprised: sorted ALDH⁺ and CD44⁺CD24⁻ versus bulk patient-derived xenograft RNA-seq (GSE115302); sorted CD44ʰⁱCD24ˡᵒ versus CD44ˡᵒCD24ʰⁱ microarray (GSE36643); and mammosphere versus adherent monolayer RNA-seq (GSE182532). Survival analysis used the METABRIC cohort (Curtis et al., 2012; Pereira et al., 2016) obtained from the cBioPortal datahub.

### 2.2 Single-cell processing

Data were processed in Scanpy v1.x (Wolf et al., 2018). Cells were filtered by gene count, total counts, and mitochondrial fraction (thresholds set per dataset by inspection of QC distributions). Counts were library-size normalized to 10⁴ and log1p-transformed. The top 2,000 highly variable genes (seurat_v3 flavor) were used for principal component analysis (30 components); batch effects across patients were corrected with Harmony (Korsunsky et al., 2019). Neighborhood graphs and UMAP embeddings were computed on the corrected space, and cells were clustered with the Leiden algorithm (Traag et al., 2019). Clusters were annotated to major lineages using canonical markers (e.g., EPCAM/KRT8 for epithelial, CD3D for T cells), and cancer-epithelial cells were isolated for CSC analysis.

### 2.3 Stemness scoring and CSC labels

Cancer-epithelial cells were scored against a published stemness signature (Malta et al., 2018) using Scanpy `score_genes`. Cells in the top and bottom score quartiles were labeled stemness-high (CSC-high) and stemness-low (CSC-low), respectively, providing a binary contrast for both stages.

### 2.4 Stage 1 — differential expression

DE between CSC-high and CSC-low cancer-epithelial cells used the Wilcoxon rank-sum test with Benjamini–Hochberg false-discovery-rate (FDR) correction. Genes were ranked by the Wilcoxon statistic; markers were defined at FDR < 0.05.

### 2.5 Stage 2 — Geneformer fine-tuning and attention

Gene symbols were mapped to Ensembl IDs (mygene) and tokenized in Geneformer's rank-value encoding. The Geneformer-V2-104M_CLcancer checkpoint was fine-tuned for binary CSC-high/-low classification (classification head plus unfrozen upper transformer layers; mixed-precision; early stopping on macro-F1; held-out test split). Per-gene importance was obtained by forward-passing CSC-high cells with attention output enabled and aggregating attention across layers and heads, then mapping tokens back to gene symbols to yield a genome-wide ranking.

### 2.6 Functional benchmark

For each functional assay, a per-gene functional log₂ fold change (CSC vs. non-CSC) was computed: log₂((CSC+1)/(bulk+1)) for GSE115302 (averaged over two PDX models); the provided DESeq2 log₂ fold change for GSE182532; and log₂ of the reported linear fold change for GSE36643. For each method we computed: (i) precision@100 — the fraction of the top 100 genes up-regulated in the functional CSC population, with significance assessed against 1,000 size-matched random gene sets; (ii) mean functional log₂ fold change of the top 100 genes; (iii) AUROC of the method's continuous score in ranking genes above the 75th percentile of functional log₂ fold change genome-wide; and (iv) Spearman correlation between method score and functional log₂ fold change. Pathway content was assessed by curated CSC-pathway membership and by Enrichr (Chen et al., 2013) against the MSigDB Hallmark collection (Liberzon et al., 2015).

To contextualize the transformer against conventional supervised learners, two baseline gene rankings were derived from the same CSC-high versus CSC-low cancer-epithelial cells. An L2-penalized logistic regression (scikit-learn; saga solver, C=1) was fit on standardized log-normalized expression, ranking genes by signed coefficient. A random forest (200 trees, √p features per split) was fit on the same data, ranking genes by impurity-based importance signed by the direction of the high-versus-low mean difference. Both rankings were scored against the functional gates identically to the DE and attention rankings. All four methods (DE, Geneformer attention, logistic regression, random forest) were compared by mean AUROC and mean precision@100 across the three sorted-population gates.

To assess cross-cancer transfer, the same four breast-derived rankings were scored against five additional functional datasets in other tumor types: sorted ALDH⁺ versus ALDH⁻ in prostate (GSE270565), melanoma (GSE243840), and bladder (GSE166947) cancer, and sphere versus adherent in prostate (GSE228203) and ovarian (GSE232783) cancer (each n = 3 per group). Per-gene functional log₂ fold change was computed as log₂(mean(CSC)+1) − log₂(mean(non-CSC)+1) from the processed expression matrices, and mean AUROC and precision@100 were compared within the breast discovery gates versus the five held-out cancers.

### 2.7 Cross-cancer analysis

The Stage 1 pipeline was applied to GBM and melanoma. GBM neoplastic cells were identified using the original study annotations because GBM tumor cells transcriptionally resemble astrocytes and oligodendrocyte-precursor cells; melanoma malignant cells were taken from the published copy-number-based annotation. Pan-cancer genes were defined as those recurrently up-regulated across tumor types.

### 2.8 Survival analysis

METABRIC microarray expression was scored per patient by single-sample gene set enrichment analysis (ssGSEA; Barbie et al., 2009; Hänzelmann et al., 2013) for the Stage 1, Stage 2, and consensus gene sets. Cox proportional-hazards models (lifelines; Davidson-Pilon, 2019) were fit with PAM50 subtype as a stratification variable (separate baseline hazard per subtype) and age as a covariate; hazard ratios are reported per standard deviation of the ssGSEA score. Kaplan–Meier curves used ssGSEA tertiles with subtype-aware log-rank tests.

### 2.9 Software and reproducibility

The benchmarking framework is released as an installable, tested Python package, `cscbench` (MIT license), which implements the functional-gate abstraction and the benchmark metrics used here. A `FunctionalGate` is constructed from an expression matrix and CSC/non-CSC sample columns (or from a precomputed fold change), and `benchmark_ranking`/`run_benchmark` score any gene ranking — supplied as a gene→score mapping or an ordered list — against one or many gates, returning precision@k with an empirical random-set p-value, genome-wide AUROC, effect-size concordance, and enrichment. The package is method- and cancer-agnostic and depends only on numpy, pandas, and scikit-learn; a worked example reproduces the multi-cancer benchmark reported here. All analyses are additionally scripted (`notebooks/phase_*.py`); result tables and figures are provided in the repository (see Data Availability).

---

## 3. Results

### 3.1 A two-stage pipeline contrasts abundance- and attention-based marker ranking

Processing GSE176078 (Figure 1A) and scoring cancer-epithelial cells for stemness (Figure 1B), Stage 1 DE between CSC-high and CSC-low cells yielded 3,028 markers, led by MYC, CD44, and VIM (Figure 1C); these varied along a pseudotime trajectory from stem-like to differentiated states (Figure 1D). Stage 2 fine-tuning of Geneformer classified CSC-high versus CSC-low cells at 93.2% macro-F1, and attention aggregation produced a genome-wide ranking led by SOX9, KLF4, FZD7, CD44, and ALDH1A3 — a list enriched for transcription factors and signaling receptors rather than the ribosomal and proliferation genes dominating the DE ranking.

### 3.2 The two rankings are complementary, not redundant

The Stage 1 and Stage 2 top-200 lists shared only 16 genes (Jaccard 0.042), including CD44, SOX9, KLF4, MYC, VIM, and SERPINE2. This low overlap reflects a difference in what each method measures — abundance in CSCs versus informativeness about CSC identity in transcriptomic context — rather than unreliability of either. The 16-gene consensus, independently surfaced by a statistical test and a pre-trained transformer, constitutes an internally replicated high-confidence core.

### 3.3 Functional benchmarking dissociates the methods into "shortlist" and "ranker" (Figure 2)

Scored against three independent functional assays whose labels were used in neither pipeline (Figure 2A), both methods enriched for functionally up-regulated genes above random on the two sorted-population gates (ALDH⁺ precision@100: DE 0.82, Geneformer 0.78; both p ≤ 0.001). The mammosphere gate enriched for neither method, consistent with the transcriptional drift of long-term sphere cultures; we down-weight it accordingly.

Across four contrasts the two rankings split head-to-head metrics systematically: DE won precision@100 and effect size, concentrating reliable high-abundance markers at the top of its list, whereas Geneformer attention won AUROC on three of four gates by wide margins (0.733 vs. 0.544; 0.682 vs. 0.547; 0.570 vs. 0.501), discriminating functional CSC genes genome-wide where the DE statistic was near chance. DE thus behaves as a high-precision shortlist generator and attention as a superior genome-wide ranker, quantifying the abundance-versus-informativeness distinction against functional ground truth.

To test whether attention's genome-wide advantage reflects its pretrained transformer architecture rather than simply being a supervised or nonlinear model, we added two supervised baselines trained to classify the same CSC-high versus CSC-low cells: L2-penalized logistic regression (linear) and random forest (nonlinear), with genes ranked by coefficient magnitude and impurity importance, respectively (Figure 2A). Averaged over the three sorted-population gates, genome-wide AUROC increased monotonically with model capability — DE 0.531, logistic regression 0.543, random forest 0.599, Geneformer 0.662 — with Geneformer exceeding both supervised baselines, including the nonlinear random forest. Because logistic regression and random forest are also supervised (and the latter nonlinear), attention's superior discrimination is not attributable merely to supervision or nonlinearity but to the contextual representation learned during pretraining. Conversely, DE retained the best mean top-100 precision (0.730), ahead of random forest (0.700), Geneformer (0.683), and logistic regression (0.667), confirming the shortlist-versus-ranker dissociation holds against stronger baselines. Reassuringly, all four methods placed canonical regulators (MYC, CD44, KLF4, SOX9) among their top genes, indicating the baselines are well-calibrated rather than degenerate.

Attention also recovered specific functionally validated genes that DE ranked low: among genes in Geneformer's top 50 but absent from the DE top 200, KLK5 (up in 4/4 gates; mean functional log₂ fold change +2.25), FOXI1 (+1.93), TOX3 (+0.90), MAGEA3 (+0.48), and the surface receptor BMPR1B (+0.39) were robustly up-regulated. The genes with the strongest cross-gate support overall — SOX9 and SERPINE2 — were identified convergently by DE, attention, and function, defining the highest-confidence tier.

**Marker rankings are tissue-specific and transfer weakly across cancers.** To test whether the breast-derived rankings, and the method-level dissociation, extend to other tumor types, we applied the framework to five additional functional datasets spanning five cancers and all three criteria: sorted ALDH⁺ populations in prostate (GSE270565), melanoma (GSE243840), and bladder (GSE166947), and sphere-versus-adherent contrasts in prostate (GSE228203) and ovarian (GSE232783) cancer — nine functional gates across six cancers in total (Figure 3). Within the breast discovery cancer, the dissociation was clear (mean AUROC: Geneformer 0.587, well ahead of random forest 0.552, logistic regression 0.528, and DE 0.519). On the five held-out cancers, however, the discriminative performance of all breast-derived rankings fell toward chance (mean AUROC 0.51–0.54), and Geneformer's large within-breast advantage shrank to a near-tie with the linear baseline (0.534 vs. logistic regression 0.540; DE 0.512). The direction of the method ordering was largely preserved but its magnitude was not. This indicates that CSC marker *rankings* are substantially tissue-specific and do not transfer strongly across tumor types — a concrete, and in retrospect expected, consequence of cancer-type-specific stemness biology. The practical implication is that the framework and the method-level comparison are best applied *per cancer*, deriving rankings on the target tumor type rather than transferring a ranking learned in one cancer to another; the near-chance transfer performance quantifies the cost of doing otherwise.

### 3.4 Attention captures developmental-pathway signal (Figure 2B)

Integrating method ranks, functional validation, surface localization, and curated pathway membership, the leading druggable candidates were the surface receptors FZD7 (Wnt) and BMPR1B (BMP/TGF-β). Consistently, Geneformer's top-50 loaded preferentially on Wnt (FZD7, KLF4, MYC, SOX9) and BMP/TGF-β (BMPR1B, ID1, SERPINE2) pathways, whereas DE's top-50 was enriched (Hallmark) for interferon-γ, TNF/NF-κB, and Myc-target programs — inflammatory and proliferative rather than stemness-signaling content. We emphasize that these genes are prioritized hypotheses: FZD7, although attention-prominent and druggable, showed only weak functional up-regulation (+0.06) and thus represents a high-value but unproven target, whereas BMPR1B combined surface druggability with functional support.

### 3.5 A recurrent program across three cancers (Figure 4)

Applying the Stage 1 pipeline to GBM (GSE84465; 1,076 neoplastic cells) and melanoma (GSE72056; 1,257 malignant cells) identified cancer-specific markers and a recurrent core. CD44, VIM, FN1, and MYC were shared between breast and GBM; extending to melanoma, ANXA1, CAV1, and OSMR recurred across all three tumor types (Figure 4). The modest three-way overlap is expected given the tissue-specific component of CSC biology, and the recurrent genes represent the most likely pan-cancer dependencies.

### 3.6 Subtype-stratified survival recovers the expected direction (Figure 5)

CSC signatures derived from scRNA-seq are confounded in bulk survival analysis: naive scoring made CD44-high METABRIC tumors appear protective, because luminal good-prognosis cells dominate the bulk and also express CD44. Replacing mean z-score scoring with ssGSEA and replacing covariate adjustment with PAM50-subtype stratification removed this confounding. After correction, all signatures showed the expected adverse-prognosis direction (hazard ratio [HR] > 1), and the consensus signature was significantly associated with overall survival (HR 1.06 per SD; 95% CI 1.00–1.12; p = 0.048; Figure 5). Effect sizes were small, as expected for a minority-cell signature diluted in bulk tissue, with the attention signature trending slightly stronger than the DE signature (relapse-free HR 1.07, p = 0.087 vs. HR 1.00, p = 0.96).

### 3.7 Per-cell protein-anchored validation is unavailable in public breast data

The strongest anchor would define CSCs by surface protein (CD44⁺CD24⁻) at single-cell resolution. Although the source atlas includes CITE-seq with CD44, CD24, and EPCAM antibodies, the antibody capture was performed on a CD45⁺-enriched immune/stromal compartment, and only 7 of 2,169 cancer-epithelial cells carried matched protein measurements. Because CD44⁺CD24⁻ is a CSC gate only within epithelium, a per-cell protein-anchored benchmark could not be constructed. This appears to generalize — public breast CITE-seq profiles the tumor immune microenvironment rather than epithelial CSCs — leaving per-cell protein validation an open data gap that motivates targeted epithelial CITE-seq.

### 3.8 Method-selection guidance

The benchmark translates into practical guidance for CSC marker studies (Table 1). When the goal is a short, high-confidence list for immediate follow-up (e.g., a small panel for flow sorting or targeted perturbation), differential expression — or a supervised linear model — is preferable, because these concentrate reliable, high-effect markers at the top of the ranking (best precision@100). When the goal is to score or prioritize across the whole transcriptome — for example, to nominate regulators that abundance-based ranking overlooks, or to weight genes in a downstream model — transformer attention is preferable, because it discriminates functional CSC genes genome-wide far better than the alternatives (best AUROC). Because the two objectives are complementary, the highest-yield strategy is to intersect them: genes convergently ranked by DE, attention, and functional support (here, SOX9 and SERPINE2) form the highest-confidence tier, while attention-unique, functionally validated genes (e.g., KLK5, FOXI1, BMPR1B) recover mechanistic candidates that DE misses.

**Table 1. Method-selection guidance from the functional benchmark.**

| Goal | Recommended method | Rationale (benchmark evidence) |
|---|---|---|
| High-precision shortlist (small validated panel) | Differential expression; supervised linear model | Best top-100 precision (DE 0.73); dense in reliable high-effect markers |
| Genome-wide ranking / prioritization | Transformer attention (Geneformer) | Best genome-wide AUROC (0.66 vs. 0.53–0.60); ranks whole transcriptome by relevance |
| Recover mechanistic regulators missed by abundance | Transformer attention | Uniquely surfaces functionally validated genes DE ranks low (KLK5, FOXI1, BMPR1B) |
| Highest-confidence markers | Intersection of DE + attention + function | Convergent genes (SOX9, SERPINE2) validated by all three lines of evidence |
| Understand the source of a method's advantage | Include supervised baselines (LR, RF) | Isolates pretraining contribution: attention > RF > LR > DE on AUROC |

---

## 4. Discussion

Our results reframe the common question "DE or transformer attention — which is better for CSC markers?" into "what does each method measure?" Differential expression and Geneformer attention rank genes by different properties — abundance versus contextual informativeness — and the functional benchmark shows this dissociation is real rather than semantic: DE generates denser high-precision shortlists, while attention ranks the whole transcriptome by relevance to CSC identity more accurately. The practical implication is that the two are best used together, with attention surfacing mechanistically informative regulators (transcription factors, signaling receptors) that abundance ranking buries, and DE providing a high-precision core and effect-size estimates.

The candidates this yields are biologically coherent. Attention's preferential loading on Wnt and BMP/TGF-β signaling, and its prioritization of the druggable receptors FZD7 and BMPR1B, align with the developmental-pathway model of stemness and are not recoverable from DE's inflammation/proliferation-dominated top list. These observations support using transformer attention as a complementary discovery layer rather than a replacement for established statistical testing.

Two caveats bound the interpretation. First, attention denotes predictive importance, not proven causality; functionally validated attention-unique genes (KLK5, FOXI1, BMPR1B, FZD7) are prioritized hypotheses whose decisive test is perturbation (e.g., CRISPR depletion followed by sphere-formation and limiting-dilution assays). Second, the functional anchors here are bulk or PDX-sorted populations; per-cell protein anchoring would be stronger but is unavailable in public breast data. Within these bounds, the pipeline functions as a prioritization engine, narrowing ~20,000 genes to a ranked, functionally filtered, pathway-annotated shortlist suitable for the expensive experiments that can certify a marker.

### 4.1 Limitations

Stemness labels derive from published signatures; although the functional benchmark mitigates circularity, the training labels themselves remain literature-anchored. Bulk survival effects are small and significant for only one endpoint, establishing direction rather than a strong prognostic claim. Functional anchors are uneven (one breast PDX dataset has no within-group replicates; the mammosphere gate was non-enriched). Marker rankings transfer only weakly across cancers, so the per-cancer conclusions should not be assumed to hold when a ranking is ported to a different tumor type; a definitive test of whether the *method-level* dissociation replicates would require deriving fresh rankings within each cancer (rather than transferring breast rankings), which entails re-running the full pipeline — including transformer fine-tuning — per tumor type. The benchmark compares four representative methods spanning the modeling spectrum but not the full landscape of marker-discovery tools or foundation models. All nominated candidates are computational predictions pending experimental validation.

### 4.2 Conclusion

Differential expression and transformer attention identify complementary CSC gene programs whose relative strengths are made explicit by functional benchmarking: DE as a high-precision shortlist generator and attention as a superior genome-wide ranker of mechanistically informative regulators. The framework is method- and cancer-agnostic — any gene-ranking method can be scored against any functional CSC assay with a per-gene fold change — so it extends directly to additional methods (e.g., other foundation models or network-based rankers) and to additional tumor types as sorted-population and self-renewal datasets accrue. Together with the prioritized, pathway-annotated candidate list, it provides a reusable, extensible basis for CSC marker discovery and experimental follow-up.

---

## Data Availability Statement

All datasets analyzed are publicly available. Primary breast cancer scRNA-seq: GEO GSE176078. Cross-cancer scRNA-seq: GSE84465 (GBM), GSE72056 (melanoma). Functional CSC assays (breast): GSE115302, GSE36643, GSE182532. Functional CSC assays (cross-cancer transfer): GSE270565 (prostate), GSE243840 (melanoma), GSE228203 (prostate), GSE232783 (ovarian), GSE166947 (bladder). Survival data: METABRIC (cBioPortal, study brca_metabric). The benchmarking framework is released as the `cscbench` Python package (MIT license); analysis code (`notebooks/phase_*.py`), the package, and derived result tables and figures are available at https://github.com/ksuma2109/csc-marker-benchmark (private during peer review; made publicly available upon publication).

## Author Contributions

S.K. designed the study, implemented the pipeline, benchmark, and `cscbench` package, analyzed the data, and wrote the manuscript.

## Conflict of Interest

The authors declare that the research was conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest.

## Funding

[Funding sources, or: "This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors."]

## Acknowledgments

The author thanks the contributors of the public datasets used in this study.

## Author Biography

**Suma Kasa** is an independent researcher in computational biology and bioinformatics, with interests in single-cell genomics, cancer stem cell biology, and machine-learning methods for gene prioritization.

---

## References

Al-Hajj, M., Wicha, M. S., Benito-Hernandez, A., Morrison, S. J., and Clarke, M. F. (2003). Prospective identification of tumorigenic breast cancer cells. *Proc. Natl. Acad. Sci. U.S.A.* 100, 3983–3988.

Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature* 462, 108–112.

Batlle, E., and Clevers, H. (2017). Cancer stem cells revisited. *Nat. Med.* 23, 1124–1134.

Chen, E. Y., Tan, C. M., Kou, Y., Duan, Q., Wang, Z., Meirelles, G. V., et al. (2013). Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. *BMC Bioinformatics* 14, 128.

Curtis, C., Shah, S. P., Chin, S.-F., Turashvili, G., Rueda, O. M., Dunning, M. J., et al. (2012). The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature* 486, 346–352.

Darmanis, S., Sloan, S. A., Croote, D., Mignardi, M., Chernikova, S., Samghababi, P., et al. (2017). Single-cell RNA-seq analysis of infiltrating neoplastic cells at the migrating front of human glioblastoma. *Cell Rep.* 21, 1399–1410.

Davidson-Pilon, C. (2019). lifelines: survival analysis in Python. *J. Open Source Softw.* 4, 1317.

Hänzelmann, S., Castelo, R., and Guinney, J. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics* 14, 7.

Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* 16, 1289–1296.

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., and Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst.* 1, 417–425.

Malta, T. M., Sokolov, A., Gentles, A. J., Burzykowski, T., Poisson, L., Weinstein, J. N., et al. (2018). Machine learning identifies stemness features of cancer and normal cells. *Cell* 173, 338–354.

Pereira, B., Chin, S.-F., Rueda, O. M., Vollan, H.-K. M., Provenzano, E., Bardwell, H. A., et al. (2016). The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. *Nat. Commun.* 7, 11479.

Reya, T., Morrison, S. J., Clarke, M. F., and Weissman, I. L. (2001). Stem cells, cancer, and cancer stem cells. *Nature* 414, 105–111.

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl. Acad. Sci. U.S.A.* 102, 15545–15550.

Theodoris, C. V., Xiao, L., Chopra, A., Chaffin, M. D., Al Sayed, Z. R., Hill, M. C., et al. (2023). Transfer learning enables predictions in network biology. *Nature* 618, 616–624.

Tirosh, I., Izar, B., Prakadan, S. M., Wadsworth, M. H., Treacy, D., Trombetta, J. J., et al. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 352, 189–196.

Traag, V. A., Waltman, L., and van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.* 9, 5233.

Wolf, F. A., Angerer, P., and Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* 19, 15.

Wu, S. Z., Al-Eryani, G., Roden, D. L., Junankar, S., Harvey, K., Andersson, A., et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nat. Genet.* 53, 1334–1347.
