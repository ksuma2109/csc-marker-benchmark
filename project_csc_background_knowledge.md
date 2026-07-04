---
name: project-csc-background-knowledge
description: "Background knowledge map for CSC bioinformatics — what the user needs to learn, organized by domain"
metadata: 
  node_type: memory
  type: project
  originSessionId: eebb9393-a09e-4a1e-97ba-44e1a57275d0
---

Background knowledge required for CSC Identification & Characterization, as discussed with user. Organized by domain.

## Biology Foundations
- Stemness properties: self-renewal, multipotency, quiescence
- Normal stem cell hierarchies (hematopoietic system as model)
- Pluripotency TFs: OCT4, SOX2, NANOG, KLF4, MYC
- Core stemness pathways: Wnt, Notch, Hedgehog, TGF-β
- Clonal evolution vs. CSC model of tumor heterogeneity
- EMT (Epithelial-Mesenchymal Transition) — linked to stemness
- Classic CSC markers: Breast (CD44+/CD24−, ALDH1+), CRC (CD133+, LGR5+), AML (CD34+/CD38−), GBM (CD133+, SOX2+)

## Epigenomics
- DNA methylation (CpG islands), histone modifications (H3K4me3, H3K27me3, H3K27ac)
- Bivalent domains as stemness hallmark
- Chromatin accessibility

## Computational Skills
- R (primary), Python, Bash/Unix
- Statistics: FDR, PCA, UMAP, t-SNE
- scRNA-seq pipeline: alignment → QC → normalization → clustering → marker ID → trajectory

## Key Tools
- Seurat, scanpy — scRNA-seq analysis
- Monocle3, Palantir, scVelo — trajectory/pseudotime
- DESeq2, edgeR — differential expression
- GSEA, fgsea — pathway enrichment
- Harmony, scVI — batch correction
- AUCell, UCell — gene signature scoring
- SingleR, CellTypist — cell type annotation

## Key Databases
- GEO, TCGA, CellxGene, DepMap, MSigDB

## Landmark Papers
- Lapidot 1994 (Nature) — first CSC demonstration in AML
- Al-Hajj 2003 (PNAS) — breast CSCs, CD44+/CD24−
- Reya 2001 (Nature) — conceptual CSC overview
- Malta 2018 (Cell) — ML stemness index (mRNAsi)
- Puram 2017 (Cell) — scRNA-seq of HNSCC, partial EMT in CSCs
- Tirosh 2016 (Science) — scRNA-seq of melanoma heterogeneity

## Suggested Learning Sequence
Month 1-2: Biology (stem cell textbook + foundational papers)
Month 2-3: Programming (R for Data Science, StatQuest)
Month 3-4: Bioinformatics methods (Seurat vignettes, OSCA book)
Month 4-6: Apply to real CSC dataset from GEO/CellxGene
