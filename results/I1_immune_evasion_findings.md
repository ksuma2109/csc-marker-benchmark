# B / I1 — CSC State and the Tumor Immune Microenvironment: Findings

**Status:** standalone analysis note — NOT in the manuscript.
**Data:** GSE176078 breast atlas (100,064 cells; 35,214 T-cells incl. CD8+/CD4+ subsets; 24,489 cancer-epithelial; 26 patients).
**Script:** `notebooks/phase_I1_immune_evasion.py` · **Outputs:** `results/tables/I1_immune_evasion.csv`, `results/tables/I1_patient_level.csv`, `results/figures/I1_immune_evasion.png`.

## Question
Is CSC state associated with immune evasion — do CSC-high malignant cells downregulate antigen presentation and/or upregulate checkpoint ligands, and do high-CSC tumors have less T-cell infiltration / more T-cell exhaustion?

## Headline: the simple "CSCs are immune-cold / evade by MHC-I loss" model is NOT supported here

### I1a — cell-level (cancer-epithelial cells; well-powered, ~24k cells)
| Program | Spearman r vs stemness | CSC-high vs CSC-low | Direction |
|---|---|---|---|
| Antigen presentation (classical MHC-I: HLA-A/B/C, TAP1/2, NLRC5) | **+0.257** (p≈0) | +0.66 vs +0.19 | CSC-high express **MORE** MHC-I |
| Checkpoint ligands (PD-L1, PD-L2, CD47) | +0.124 (p≈0) | −0.18 vs −0.23 | CSC-high weakly **higher** |

- CSC-high malignant cells express **higher** classical antigen-presentation machinery, not lower — the *opposite* of the textbook MHC-I-loss evasion model. (Robust: B2M and HLA-E were **excluded** from the MHC-I score because they are themselves CSC consensus genes, so this is not a shared-gene artifact.)
- CSC-high cells show a **weak** upregulation of checkpoint ligands (PD-L1/PD-L2/CD47) — a small effect (r=+0.12) that is the only evasion-consistent signal.

### I1b — patient-level (20 patients with sufficient cells)
| Association (stemness vs) | Spearman r | p |
|---|---|---|
| T-cell infiltration | −0.068 | 0.78 (null) |
| CD8 exhaustion | +0.319 | 0.17 (trend) |
| CD8 cytotoxicity | +0.435 | 0.056 (borderline) |

- **No association** between tumor stemness and T-cell infiltration — high-CSC tumors are not immune-cold.
- Weak/borderline **positive** trends for CD8 exhaustion and cytotoxicity — i.e., high-CSC tumors trend toward *more* immune activity, not less.

## Interpretation (honest, with confounds)

1. **Against simple immune evasion.** In this dataset CSC state does not fit the classic immune-cold / antigen-presentation-loss phenotype. CSC-high cells retain (even elevate) MHC-I, and T-cell infiltration is uncorrelated with stemness. This pushes back on an oversimplified narrative.
2. **The evasion signal, if any, is on the checkpoint axis.** The only evasion-consistent finding is a weak upregulation of checkpoint ligands (PD-L1/PD-L2/CD47) in CSC-high cells — suggesting that in breast cancer, CSC immune modulation (if present) acts through co-inhibitory checkpoints rather than antigen-presentation loss.
3. **Subtype is a major confound.** TNBC is both the most stem-like (mean stemness: TNBC +0.30 vs ER+ −0.20, HER2+ −0.15) and the most immune-infiltrated breast subtype. The borderline positive patient-level trends (exhaustion, cytotoxicity) are therefore likely subtype-driven rather than a direct CSC effect. With only 20 patients, these cannot be disentangled here.

## Translational relevance (connects to the CAR-T / immunotherapy angle)
- Because CSC-high cells **retain antigen presentation**, they are in principle **still visible to T-cells / CAR-T** — CSCs do not obviously escape by becoming MHC-I-null.
- The weak **checkpoint-ligand** upregulation suggests that combining CSC-directed immunotherapy with **checkpoint blockade** is the more plausible evasion axis to address, rather than antigen-presentation restoration.

## Caveats
- Patient-level analysis is **underpowered (n=20)** and **subtype-confounded**; correlations are trends, not robust associations.
- Cell-level results are well-powered but correlational; module scores depend on gene-set choice.
- This is bulk-of-atlas scRNA-seq; spatial proximity of CSCs to T-cells (which matters for evasion) is not assessed.

## Possible next steps
- Control for subtype explicitly (within-TNBC analysis, or partial correlation) — limited by n.
- Use the immune **CITE-seq** protein data (already downloaded) to score T-cell exhaustion at the protein level.
- Spatial transcriptomics (GSE176078 includes Visium) to test CSC–T-cell spatial exclusion directly.
