# Build a journal-styled PDF of the manuscript with figures embedded.
# markdown -> HTML (+CSS) -> PDF via weasyprint.

import os, re, base64
import markdown
from weasyprint import HTML

ROOT = os.path.abspath(".")
MD   = "manuscript/manuscript_draft.md"
OUT  = "manuscript/CSC_manuscript.pdf"
FIGDIR = "manuscript/figures"

text = open(MD).read()
# strip the leading HTML comment block (editor note)
text = re.sub(r"^<!--.*?-->\s*", "", text, count=1, flags=re.DOTALL)

# convert body markdown -> html
body = markdown.markdown(text, extensions=["tables", "sane_lists", "attr_list"])

# figure captions (match the composite figures produced)
FIGURES = [
    ("Figure1.png", "Figure 1.", "Single-cell identification of breast cancer stem cells (GSE176078). "
        "(A) Cell-type annotated UMAP of all 100,064 cells. (B) Stemness score on cancer-epithelial cells. "
        "(C) Stage 1 differential expression (CSC-high vs CSC-low). (D) Pseudotime trajectory."),
    ("Figure2.png", "Figure 2.", "Functional benchmarking of CSC marker-ranking methods. "
        "(A) Four-method comparison (DE, Geneformer attention, logistic regression, random forest) against "
        "ALDH+, CD44+CD24-, and mammosphere gates: genome-wide AUROC and top-marker precision; attention "
        "exceeds both supervised baselines on AUROC while DE retains best precision. (B) DE vs Geneformer detail. "
        "(C) Integrated candidate prioritization and CSC-pathway loading by method."),
    ("Figure3.png", "Figure 3.", "Multi-cancer functional benchmark: marker rankings are tissue-specific. "
        "Genome-wide AUROC of the four breast-derived rankings against nine functional gates across six cancers, "
        "and mean AUROC within breast (discovery) versus five held-out cancers — attention's within-breast lead "
        "does not transfer, indicating CSC marker rankings are substantially tissue-specific."),
    ("Figure4.png", "Figure 4.", "Cross-cancer validation identifies a recurrent CSC program. "
        "(A) GBM neoplastic-cell CSC analysis (GSE84465). (B) Three-cancer evidence matrix (breast, GBM, melanoma)."),
    ("Figure5.png", "Figure 5.", "Subtype-stratified survival association of CSC signatures (METABRIC, n=1,974). "
        "(A) Kaplan-Meier by ssGSEA stemness tertile. (B) Subtype-stratified Cox hazard ratios per SD."),
]

def img_data_uri(path):
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode()
    return f"data:image/png;base64,{b64}"

fig_html = ['<h2 class="figures-header">Figures</h2>']
for fname, num, cap in FIGURES:
    p = os.path.join(FIGDIR, fname)
    if not os.path.exists(p):
        continue
    fig_html.append(
        f'<figure><img src="{img_data_uri(p)}"/>'
        f'<figcaption><strong>{num}</strong> {cap}</figcaption></figure>')
fig_html = "\n".join(fig_html)

CSS = """
@page { size: A4; margin: 2cm 2.2cm; @bottom-center { content: counter(page);
        font-family: Georgia, serif; font-size: 9pt; color:#666; } }
body { font-family: Georgia, 'Times New Roman', serif; font-size: 10.5pt;
       line-height: 1.5; color:#1a1a1a; text-align: justify; }
h1 { font-size: 17pt; text-align:center; line-height:1.3; margin: 0 0 4pt 0;
     font-family: 'Helvetica Neue', Arial, sans-serif; }
h2 { font-size: 12.5pt; margin: 16pt 0 4pt 0; border-bottom:1px solid #ccc;
     padding-bottom:2pt; font-family:'Helvetica Neue', Arial, sans-serif; }
h3 { font-size: 11pt; margin: 12pt 0 3pt 0; color:#222;
     font-family:'Helvetica Neue', Arial, sans-serif; }
p { margin: 0 0 7pt 0; }
strong { color:#000; }
hr { border:none; border-top:1px solid #ddd; margin:10pt 0; }
em { color:#333; }
table { border-collapse: collapse; width:100%; font-size:9.5pt; margin:8pt 0; }
th,td { border:1px solid #bbb; padding:3pt 5pt; text-align:left; }
th { background:#f0f0f0; }
code { font-family: 'Courier New', monospace; font-size:9pt; background:#f5f5f5; padding:0 2pt; }
figure { margin: 14pt 0; page-break-inside: avoid; text-align:center; }
figure img { max-width: 100%; height:auto; border:1px solid #e0e0e0; }
figcaption { font-size: 9pt; text-align:justify; margin-top:5pt; color:#333; line-height:1.4; }
.figures-header { page-break-before: always; }
"""

html = f"""<!doctype html><html><head><meta charset="utf-8">
<style>{CSS}</style></head><body>{body}{fig_html}</body></html>"""

HTML(string=html, base_url=ROOT).write_pdf(OUT)
sz = os.path.getsize(OUT) / 1024
print(f"Saved {OUT}  ({sz:.0f} KB)")
