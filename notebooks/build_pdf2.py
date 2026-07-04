# Build a clean single-column journal-style PDF of the second (immune) manuscript.
import os, re, base64
import markdown
from weasyprint import HTML

ROOT = os.path.abspath(".")
MD   = "manuscript2/immune_manuscript_draft.md"
OUT  = "manuscript2/CSC_immune_manuscript.pdf"
FIGDIR = "manuscript2/figures"

text = open(MD).read()
text = re.sub(r"^<!--.*?-->\s*", "", text, count=1, flags=re.DOTALL)
body = markdown.markdown(text, extensions=["tables", "sane_lists", "attr_list"])

# center the title block (title..correspondence, up to first hr)
body = re.sub(r"(</h1>)(.*?)(<hr\s*/?>)", r'\1<div class="titleblock">\2</div>\3',
              body, count=1, flags=re.DOTALL)

FIGURES = [
    ("Figure1.png", "Figure 1.",
     "TNBC cancer stem cells upregulate CD47 and retain antigen presentation "
     "(subtype-specific, independently replicated). "
     "(A) Subtype specificity (discovery, Wu 2021): the CSC innate-evasion phenotype "
     "(CD47 up, MHC-I retained) holds within TNBC and HER2+ but is absent/reversed in ER+. "
     "(B) Replication across two independent TNBC cohorts (Pal 2021; Gao 2020): CD47 upregulation "
     "and MHC-I retention replicate; PD-L1 stays near zero. Points = individual tumours."),
]

def data_uri(p):
    with open(p, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

figs = ['<h2 class="figs">Figure</h2>']
for fn, num, cap in FIGURES:
    p = os.path.join(FIGDIR, fn)
    if os.path.exists(p):
        figs.append(f'<figure><img src="{data_uri(p)}"/>'
                    f'<figcaption><strong>{num}</strong> {cap}</figcaption></figure>')
figs = "\n".join(figs)

CSS = """
@page { size: A4; margin: 2cm 2.2cm; @bottom-center { content: counter(page);
        font-family: Georgia, serif; font-size: 9pt; color:#666; } }
body { font-family: Georgia, 'Times New Roman', serif; font-size: 10.5pt;
       line-height: 1.5; color:#111; text-align: justify; }
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

html = f'<!doctype html><html><head><meta charset="utf-8"><style>{CSS}</style></head><body>{body}{figs}</body></html>'
HTML(string=html, base_url=ROOT).write_pdf(OUT)
print(f"Saved {OUT}  ({os.path.getsize(OUT)/1024:.0f} KB)")
