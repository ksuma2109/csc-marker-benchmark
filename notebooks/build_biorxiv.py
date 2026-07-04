# Build clean single-column preprint PDFs (bioRxiv-appropriate) for both manuscripts.
# Single column, figures + legends embedded, all required statements, "Preprint" banner.
import os, re, base64
import markdown
from weasyprint import HTML

ROOT = os.path.abspath(".")

PAPERS = {
    "paper1": {
        "md": "manuscript/manuscript_draft.md",
        "out": "manuscript/CSC_marker_benchmark_bioRxiv.pdf",
        "figdir": "manuscript/figures",
        "figures": [
            ("Figure1.png", "Figure 1.", "Single-cell identification of breast cancer stem cells (GSE176078). "
             "(A) Cell-type annotated UMAP. (B) Stemness score on cancer-epithelial cells. "
             "(C) Stage 1 differential expression. (D) Pseudotime trajectory."),
            ("Figure2.png", "Figure 2.", "Functional benchmarking of CSC marker-ranking methods. "
             "(A) Four-method comparison against ALDH+, CD44+CD24-, and mammosphere gates (genome-wide AUROC and "
             "top-marker precision). (B) DE vs Geneformer detail. (C) Candidate prioritization and CSC-pathway loading."),
            ("Figure3.png", "Figure 3.", "Multi-cancer functional benchmark: marker rankings are tissue-specific. "
             "Genome-wide AUROC of breast-derived rankings across nine functional gates / six cancers, and mean AUROC "
             "within breast (discovery) vs. five held-out cancers."),
            ("Figure4.png", "Figure 4.", "Cross-cancer validation identifies a recurrent CSC program. "
             "(A) GBM neoplastic-cell CSC analysis (GSE84465). (B) Three-cancer evidence matrix."),
            ("Figure5.png", "Figure 5.", "Subtype-stratified survival association of CSC signatures (METABRIC, n=1,974). "
             "(A) Kaplan-Meier by ssGSEA stemness tertile. (B) Subtype-stratified Cox hazard ratios per SD."),
        ],
    },
    "paper2": {
        "md": "manuscript2/immune_manuscript_draft.md",
        "out": "manuscript2/CSC_immune_bioRxiv.pdf",
        "figdir": "manuscript2/figures",
        "figures": [
            ("Figure1.png", "Figure 1.", "TNBC cancer stem cells upregulate CD47 and retain antigen presentation "
             "(subtype-specific, independently replicated). (A) Subtype specificity (discovery, Wu 2021). "
             "(B) Replication across two independent TNBC cohorts (Pal 2021; Gao 2020); points = individual tumours."),
        ],
    },
}

CSS = """
@page { size: A4; margin: 2cm 2.2cm;
        @top-center { content: "Preprint — not peer reviewed"; font-family: Arial, sans-serif;
                      font-size: 7.5pt; color:#999; }
        @bottom-center { content: counter(page); font-family: Georgia, serif; font-size: 9pt; color:#666; } }
body { font-family: Georgia, 'Times New Roman', serif; font-size: 10.5pt;
       line-height: 1.55; color:#111; text-align: justify; }
h1 { font-size: 17pt; text-align:center; line-height:1.3; margin:0 0 6pt 0;
     font-family:'Helvetica Neue', Arial, sans-serif; }
.titleblock { text-align:center; } .titleblock p { text-align:center; margin:2pt 0; }
h2 { font-size: 12.5pt; margin:15pt 0 4pt; border-bottom:1px solid #ccc; padding-bottom:2pt;
     font-family:'Helvetica Neue', Arial, sans-serif; }
h3 { font-size: 11pt; margin:11pt 0 3pt; font-family:'Helvetica Neue', Arial, sans-serif; color:#222; }
p { margin:0 0 7pt; } strong { color:#000; } em { color:#333; }
hr { border:none; border-top:1px solid #ddd; margin:9pt 0; }
table { border-collapse:collapse; width:100%; font-size:9.5pt; margin:8pt 0; }
th,td { border:1px solid #bbb; padding:3pt 5pt; text-align:left; } th { background:#eef; }
code { font-family:'Courier New',monospace; font-size:9pt; background:#f5f5f5; padding:0 2pt; }
.figs { page-break-before: always; }
figure { margin:12pt 0; text-align:center; page-break-inside:avoid; }
figure img { max-width:100%; border:1px solid #ddd; }
figcaption { font-size:9pt; text-align:justify; margin-top:5pt; color:#222; line-height:1.4; }
"""

def data_uri(p):
    with open(p, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

for name, cfg in PAPERS.items():
    text = open(cfg["md"]).read()
    text = re.sub(r"^<!--.*?-->\s*", "", text, count=1, flags=re.DOTALL)
    body = markdown.markdown(text, extensions=["tables", "sane_lists", "attr_list"])
    body = re.sub(r"(</h1>)(.*?)(<hr\s*/?>)", r'\1<div class="titleblock">\2</div>\3',
                  body, count=1, flags=re.DOTALL)
    figs = ['<h2 class="figs">Figures</h2>'] if cfg["figures"] else []
    for fn, num, cap in cfg["figures"]:
        p = os.path.join(cfg["figdir"], fn)
        if os.path.exists(p):
            figs.append(f'<figure><img src="{data_uri(p)}"/>'
                        f'<figcaption><strong>{num}</strong> {cap}</figcaption></figure>')
    html = (f'<!doctype html><html><head><meta charset="utf-8"><style>{CSS}</style></head>'
            f'<body>{body}{"".join(figs)}</body></html>')
    HTML(string=html, base_url=ROOT).write_pdf(cfg["out"])
    print(f"Saved {cfg['out']}  ({os.path.getsize(cfg['out'])/1024:.0f} KB)")
