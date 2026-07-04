# Cover Letter — Paper 1 (Briefings in Bioinformatics)

Dear Editors,

Please consider our manuscript, **"A functional benchmarking framework for cancer stem cell marker identification: differential expression versus transformer attention in single-cell RNA-sequencing,"** for publication as a Research Article in *Briefings in Bioinformatics*.

Computational identification of cancer stem cell (CSC) marker genes from single-cell RNA-sequencing is evaluated circularly: stemness labels are defined from published markers, and candidate genes are then judged against the same literature. No framework has existed to compare CSC marker-ranking methods against a label-independent standard. We address this gap with a reusable **functional-benchmarking framework** that scores any gene-ranking method against independent functional CSC assays (sorted ALDH⁺, CD44⁺CD24⁻, CD133⁺ populations and sphere self-renewal), and apply it systematically to four representative methods spanning the modelling spectrum.

Our study offers three contributions of interest to the journal's readership:

1. **A label-independent benchmarking framework**, released as an installable, tested Python package (`cscbench`), that is method- and cancer-agnostic and reproduces our reported results.
2. **A quantitative dissociation** between differential expression (a high-precision shortlist generator) and transformer attention (a superior genome-wide ranker), with supervised baselines showing the transformer's advantage derives from pretraining rather than supervision alone.
3. **A multi-cancer evaluation** across nine functional gates and six cancers demonstrating that CSC marker rankings are substantially tissue-specific, with explicit method-selection guidance.

This work is a comparative/benchmarking contribution that fits the scope of *Briefings in Bioinformatics*. All data are public, and code and the framework are openly available. The manuscript has not been published elsewhere and is not under consideration by another journal. The author declares no competing interests.

Thank you for your consideration.

Sincerely,
Suma Kasa (Independent Researcher)
ksuma2109@gmail.com
