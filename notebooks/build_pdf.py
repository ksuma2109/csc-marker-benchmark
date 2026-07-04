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

# Split into full-width front matter (title..Key Points) and two-column body
# (Introduction onward), Oxford/Briefings-in-Bioinformatics style.
MD_EXT = ["tables", "sane_lists", "attr_list"]
split_marker = "## 1. Introduction"
if split_marker in text:
    front_md, body_md = text.split(split_marker, 1)
    body_md = split_marker + body_md
else:
    front_md, body_md = text, ""
front = markdown.markdown(front_md, extensions=MD_EXT)
body  = markdown.markdown(body_md,  extensions=MD_EXT)

# Center only the title block (running title / authors / affiliation / correspondence):
# everything between the </h1> and the first horizontal rule.
front = re.sub(r"(</h1>)(.*?)(<hr\s*/?>)",
               r'\1<div class="titleblock">\2</div>\3', front, count=1, flags=re.DOTALL)
# Box the Key Points list (a Briefings in Bioinformatics requirement)
front = re.sub(r"(<h2>Key Points</h2>\s*<ul>.*?</ul>)",
               r'<div class="keypoints">\1</div>', front, flags=re.DOTALL)

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
@page { size: A4; margin: 1.6cm 1.5cm;
        @bottom-center { content: counter(page);
            font-family: 'Times New Roman', serif; font-size: 8.5pt; color:#666; } }
body { font-family: 'Times New Roman', Georgia, serif; font-size: 9.5pt;
       line-height: 1.42; color:#111; text-align: justify; }

/* ── Full-width front matter ── */
.front { }
h1 { font-size: 18pt; text-align:center; line-height:1.28; margin: 0 0 6pt 0;
     font-family: 'Helvetica Neue', Arial, sans-serif; font-weight:700; }
.front p { margin: 3pt 0; text-align:justify; }   /* abstract/keywords justified */
.titleblock { text-align:center; margin: 2pt 0 0; }
.titleblock p { text-align:center; margin: 2pt 0; }
.front h2 { font-size: 11pt; text-transform:uppercase; letter-spacing:0.4pt;
     margin: 12pt 0 3pt 0; font-family:'Helvetica Neue', Arial, sans-serif;
     color:#333; text-align:left; }
.front hr { border:none; border-top:1px solid #ccc; margin:8pt 0; }

/* Key Points box (Briefings in Bioinformatics requirement) */
.keypoints { border:1px solid #4a6785; background:#eef3f8; border-radius:4px;
     padding:8pt 12pt 6pt; margin:10pt 0; break-inside: avoid; }
.keypoints h2 { margin-top:0; color:#274763; border:none; }
.keypoints ul { margin:0; padding-left:16pt; }
.keypoints li { margin-bottom:4pt; text-align:justify; }

/* ── Two-column body (Oxford / BiB article style) ── */
.body { column-count: 2; column-gap: 0.7cm; column-fill: auto; }
.body h2 { font-size: 10.5pt; margin: 11pt 0 3pt 0; break-after: avoid;
     font-family:'Helvetica Neue', Arial, sans-serif; color:#1a3a5c; }
.body h3 { font-size: 9.8pt; margin: 8pt 0 2pt 0; break-after: avoid;
     font-family:'Helvetica Neue', Arial, sans-serif; color:#222; font-style:italic; }
.body p { margin: 0 0 5pt 0; }
.body hr { display:none; }

/* Tables span both columns */
.body table, table { border-collapse: collapse; width:100%; font-size:8.3pt;
     margin:6pt 0; column-span: all; break-inside: avoid; }
th,td { border:1px solid #bbb; padding:2.5pt 4pt; text-align:left; }
th { background:#e8eef4; }

strong { color:#000; }
em { color:#222; }
code { font-family: 'Courier New', monospace; font-size:8.3pt;
       background:#f4f4f4; padding:0 2pt; }

/* ── Figures: full page width after the body ── */
.figures-header { break-before: page; font-size: 12pt; text-transform:uppercase;
     letter-spacing:0.5pt; font-family:'Helvetica Neue', Arial, sans-serif; }
figure { margin: 12pt 0; break-inside: avoid; text-align:center; }
figure img { max-width: 100%; height:auto; border:1px solid #ddd; }
figcaption { font-size: 8.5pt; text-align:justify; margin-top:4pt; color:#222;
     line-height:1.35; }
"""

html = f"""<!doctype html><html><head><meta charset="utf-8">
<style>{CSS}</style></head><body>
<div class="front">{front}</div>
<div class="body">{body}</div>
{fig_html}</body></html>"""

HTML(string=html, base_url=ROOT).write_pdf(OUT)
sz = os.path.getsize(OUT) / 1024
print(f"Saved {OUT}  ({sz:.0f} KB)")
