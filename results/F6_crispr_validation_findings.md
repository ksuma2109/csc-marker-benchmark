# F6 — CRISPR-Dependency Validation (DepMap 24Q2): Findings

**Status:** standalone analysis note — NOT yet incorporated into the manuscript.
**Data:** DepMap 24Q2 genome-wide CRISPR-knockout fitness screens (Chronos gene effect), 51 breast cancer cell lines with CRISPR data (of 91 breast lines), 18,443 genes.
**Script:** `notebooks/phase_F6_crispr_validation.py` · **Outputs:** `results/tables/F6_crispr_validation.csv`, `results/tables/F6_candidate_dependencies.csv`, `results/figures/F6_crispr_validation.png`.

## Question
Do the CSC marker predictions (DE, Geneformer attention, LogReg, RandomForest) hold up under independent CRISPR perturbation — i.e., are predicted markers enriched for CRISPR *dependencies* in breast cancer cell lines, and does the Geneformer ranking capture more dependency signal than DE?

## Key result 1 — Only MYC is a fitness dependency; the CSC markers are not

| Gene | mean CRISPR gene effect (breast) | dependency (<−0.5)? | DE rank | Geneformer rank |
|---|---|---|---|---|
| **MYC** | −1.76 | **YES** | 1 | 26 |
| EGFR | −0.26 | no | 560 | 1999 |
| SOX9 | −0.22 | no | 14 | 1 |
| EPCAM | −0.20 | no | 69 | 22 |
| CD44 | −0.07 | no | 2 | 4 |
| KLK5 | −0.05 | no | 474 | 34 |
| KLF4 | −0.01 | no | 52 | 2 |
| BMPR1B | 0.00 | no | – | 14 |
| VIM | 0.01 | no | 3 | 21 |
| OSMR | 0.04 | no | 1104 | 71 |
| FN1 | 0.04 | no | 78 | 135 |
| FOXI1 | 0.08 | no | 2089 | 11 |
| FZD7 | 0.10 | no | 851 | 3 |
| ALDH1A3 | 0.16 | no | 1400 | 5 |
| SERPINE2 | 0.17 | no | 46 | 6 |

## Key result 2 — The DE-vs-Geneformer dissociation reproduces on the perturbation axis

| Method | dependency-fraction@100 (precision) | AUROC (essentiality, genome-wide) |
|---|---|---|
| Stage1_DE | **0.33** (3.5× random, p<0.001) | 0.50 (chance) |
| RandomForest | 0.17 (1.8×, p=0.008) | 0.58 |
| Stage2_Geneformer | 0.05 (below random, p=0.96) | **0.63** (best) |
| LogReg | 0.02 (below random) | 0.47 |

Same pattern as the functional-sorting benchmark: DE wins precision, Geneformer wins genome-wide AUROC — now on a third, independent (CRISPR perturbation) data type.

## Honest interpretation — NOT a clean "validated" result

1. **CSC markers are mostly NOT CRISPR fitness dependencies — and that is expected and correct.** DepMap measures 2D proliferation/fitness, not stemness/self-renewal. Genes like SOX9, FZD7, KLF4 are **cell-state regulators**: knocking them out does not stop growth in standard culture; it changes identity/self-renewal, which DepMap does not assay. Their near-zero gene effects are exactly what a stemness regulator (not an essential gene) should show.

2. **DepMap is therefore the *wrong* perturbation assay to validate stemness regulators.** It confirms MYC (a proliferation driver) but is blind to state regulators — most of the interesting CSC biology. The real lesson: the definitive perturbation test needs **stemness-context CRISPRi with a self-renewal / sphere readout** (e.g., the glioma-CSC CRISPRi platform of S. John Liu, UCSF), not bulk-fitness screens.

3. **DE's high dependency-precision (0.33) is partly a confound.** DE ranks abundant genes high, and abundant genes (ribosomal proteins, MYC) are disproportionately common-essential — so DE's "essentiality enrichment" partly reflects housekeeping essentiality, not CSC biology.

4. **Geneformer's low dependency-precision (0.05) is consistent with its markers being state regulators, not fitness genes** — the biologically coherent outcome.

## What this validation establishes
- ✅ Reproduces the method dissociation on an independent perturbation axis.
- ✅ Confirms MYC as a genuine dependency.
- ✅ Makes a defensible scientific point: CSC state regulators are non-essential in bulk fitness screens; validating them requires stemness-specific perturbation (self-renewal readouts), not DepMap-style screens.
- ❌ Does NOT show the markers are essential genes — nor should it, because they are state regulators.

## Caveats
- Only 51/91 breast lines had CRISPR data.
- Bulk-fitness CRISPR is a coarse proxy; a truer computational version would use published CSC/sphere-context CRISPR screens (self-renewal readout), which remains a possible follow-up.

## Decision
Kept as a standalone note for now (not folded into the manuscript). Its clearest use, if added later, is a short perturbation-validation subsection framed around "why stemness-context validation is needed," which also motivates the wet-lab collaboration (Liu, glioma-CSC CRISPRi).
