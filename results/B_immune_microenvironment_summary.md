# Extension B — CSC State and the Tumor Immune Microenvironment: Consolidated Findings

**Status:** standalone analysis (I1–I3) — NOT in the manuscript.
**Data:** GSE176078 breast atlas — scRNA-seq (100k cells, 26 patients), Visium spatial (6 sections, Zenodo 4739739), CITE-seq immune protein (SCP1039).
**Scripts:** `notebooks/phase_I1_immune_evasion.py`, `phase_I2_spatial_exclusion.py`, `phase_I3_protein_exhaustion.py`.

## Overall conclusion
**Across three independent analyses (transcriptomic, spatial, and protein), CSC state is NOT associated with a consistent immune-evasion phenotype in this breast cancer cohort.** The naive models — CSCs are immune-cold, downregulate MHC-I, or spatially exclude T-cells — are not supported. The only evasion-consistent signal is a weak upregulation of checkpoint ligands. This is a coherent, honest negative/nuanced result.

---

## I1 — Transcriptomic (scRNA-seq)
**Cell-level (24k cancer-epithelial cells, well-powered):**
- CSC-high cells express **more** classical MHC-I antigen-presentation machinery (Spearman r=+0.26; robust after excluding B2M/HLA-E which are themselves CSC genes) — the *opposite* of MHC-I-loss evasion.
- CSC-high cells show a **weak** upregulation of checkpoint ligands (PD-L1/PD-L2/CD47; r=+0.12) — the only evasion-consistent signal.

**Patient-level (20 patients):**
- No association between stemness and T-cell infiltration (r=−0.07, ns).
- Borderline *positive* trends for CD8 exhaustion (r=+0.32, p=0.17) and cytotoxicity (r=+0.44, p=0.056) — likely **subtype-confounded** (TNBC is both most stem-like, +0.30 vs ER+ −0.20, and most immune-infiltrated).

## I2 — Spatial (Visium, 6 sections)
Per-spot correlation of stemness vs T-cell score, and stemness in cancer spots with vs without lymphocytes:
- **No consistent spatial exclusion.** 2/6 sections negative (exclusion), 2/6 positive (co-localization), 2/6 null; mean Spearman r≈+0.06 (essentially zero).
- Notably CID4535 (ER) shows strong **co-localization** (r=+0.50); CID4290 (ER) cancer spots *with* lymphocytes are *more* stem-like — opposite of exclusion.
- Conclusion: CSC-high regions do not consistently exclude T-cells; the pattern is patient-specific and averages to null, reinforcing I1.
- Caveat: Visium spots (~55 µm) are ~10-cell mixtures, so this measures regional co-occurrence, not single-cell contact.

## I4 — Which immune-evasion axis do CSCs upregulate? (per-ligand)
Breaking the I1 combined checkpoint score into individual tumour-intrinsic evasion ligands (CSC-high vs CSC-low cancer cells), ranked by effect:

| Ligand | Axis / target | CSC-high − low | % expr high/low |
|---|---|---|---|
| HLA-E | NKG2A (NK/CD8 inhibition) → anti-NKG2A (monalizumab) | **+0.56** ⚠ | 64% / 40% |
| **CD47** | macrophage "don't eat me" → **anti-CD47 (magrolimab)** | **+0.13** | 56% / 43% |
| NT5E (CD73) | adenosine → anti-CD73 | +0.11 | 12% / 0% |
| CD276 (B7-H3) | B7 family → anti-B7-H3 | +0.06 | 27% / 14% |
| **CD274 (PD-L1)** | **PD-1 axis** | **+0.009** (negligible) | 3% / 1% |
| CD24 (sanity) | Siglec-10 "don't eat me" | **−0.29** (expected ↓) | 80% / 89% |

**Key findings:**
- **The dominant CSC-associated axis is NOT the classic PD-1/PD-L1 pathway** — PD-L1 is essentially unchanged (+0.009). This argues *against* anti-PD-1/PD-L1 as the natural CSC combination partner.
- The actionable CSC-upregulated evasion ligands are **innate-immune**: **CD47** (macrophage phagocytosis checkpoint → magrolimab), **HLA-E** (NK/CD8 inhibition via NKG2A → monalizumab), and **CD73** (adenosine). HLA-E is flagged (⚠) because it is itself a CSC-consensus gene, so its association is partly circular — but the NKG2A-axis interpretation remains plausible.
- **Sanity check passes:** CD24 is *down* in CSC-high cells (−0.29), as required by the CD44⁺CD24⁻/low definition — validating the stemness labelling. Note the nuance: losing CD24 (itself a "don't-eat-me" Siglec-10 ligand) would make CSCs *more* phagocytosable, partly opposing their CD47 upregulation.

**Interpretation:** in breast cancer, CSC immune modulation operates through **innate-immune evasion (CD47/macrophage, HLA-E/NK)** rather than the T-cell PD-1 axis. Effect sizes are small (except HLA-E), correlational, and subtype-confounded.

## I3 — Protein cross-validation (CITE-seq)
Validates the I1 RNA exhaustion signature at the protein level in T-cells (ADT: PD-1, CTLA-4, TIGIT, LAG-3, TIM-3):
- RNA vs protein exhaustion correlate positively in 3/4 tumour samples; **significant in both well-powered samples** (CID3838 n=1351, r=+0.20; CID4515 n=419, r=+0.41).
- Confirms the RNA exhaustion score is a meaningful proxy — the modest r (~0.2–0.4) is expected for RNA-vs-surface-protein in sparse scRNA-seq.

---

## Translational relevance — two therapeutic hypotheses, correctly scoped

Critically, "T-cell targeting" splits into two mechanisms with **opposite** dependence on MHC-I, and the findings apply differently to each:

| Approach | Needs MHC-I? | What the findings say |
|---|---|---|
| Endogenous T-cells / checkpoint / TCR-T / vaccines | **Yes** (TCR reads peptide on MHC-I) | ✅ Encouraging — CSCs retain MHC-I, so they stay recognizable |
| **CAR-T** | **No** (CAR binds surface antigen, bypassing MHC) | ➖ MHC-I irrelevant — CAR-T needs a CSC surface antigen (the FZD7/BMPR1B surface-marker work), a separate asset |

**Hypothesis 1 — endogenous/TCR-based immunotherapy of CSCs is not defeated by antigen-loss.** CSC-high cells retain (even elevate) MHC-I antigen-presentation machinery, so the relapse-driving cells are *not* immunologically invisible. This is mildly contrarian to the common "CSCs evade via MHC-I loss" narrative and, if validated, is genuinely encouraging for T-cell/checkpoint approaches. Necessary but not sufficient — MHC-I presence does not guarantee killing.

**Hypothesis 2 — the CSC evasion axis is INNATE, not PD-1.** The per-ligand analysis (I4) shows CSCs do **not** rely on PD-L1; instead they upregulate **CD47** (anti-macrophage) and **HLA-E** (anti-NK), with **CD73** (adenosine) a secondary signal. The nominated combination partners are therefore **anti-CD47 (magrolimab)** and/or **anti-NKG2A (monalizumab)** — engaging macrophages and NK cells — rather than the default anti-PD-1/PD-L1. The lost CD24 "don't-eat-me" signal in CSCs (I4) further points to a **macrophage-phagocytosis** vulnerability.

**For CAR-T specifically:** MHC-I and checkpoint biology are moot; the relevant asset is a CSC-specific surface antigen from the surface-marker prioritization (FZD7, BMPR1B), independent of this immune analysis.

**Spatial:** the absence of consistent T-cell exclusion argues against a physical "immune-excluded CSC niche" in this cohort.

*(All of the above are computational hypotheses from one correlational, subtype-confounded cohort — testable, not established.)*

## Honest limitations
- Patient-level (n=20) and spatial (n=6 sections) are underpowered and **subtype-confounded**; several signals are trends, not robust associations.
- All analyses are correlational; module scores depend on gene-set choice.
- CITE-seq cannot test the CSC–exhaustion link directly (only CID4515 has both tumour and immune-protein data) — I3 is a signature validation, not a CSC test.

## Possible next steps
- Explicit subtype control (within-TNBC, or partial correlation) — limited by n.
- Neighborhood/distance analysis using Visium coordinates (already downloaded) for true spatial-contact exclusion.
- Larger spatial cohorts (e.g., newer breast spatial atlases) to power the CSC–immune spatial question.
