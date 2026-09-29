#!/usr/bin/env python3
"""论文笔记 Markdown（带 YAML frontmatter）→ 排版好的 PDF。用法：md2pdf_note.py in.md out.pdf"""
import sys, os, re, subprocess, yaml, markdown, html as H

src, out = sys.argv[1], sys.argv[2]
raw = open(src, encoding="utf-8").read()
m = re.match(r"^---\n(.*?)\n---\n(.*)$", raw, re.S)
meta, body = (yaml.safe_load(m.group(1)), m.group(2)) if m else ({}, raw)
body = re.sub(r"^# .*\n", "", body.lstrip(), count=1)            # 标题单独排
body = re.sub(r"\[\[(.+?)\]\]", r'<span class="wl">\1</span>', body)
body = re.sub(r"([A-Za-z])\u0304", r'<span class="bar">\1</span>', body)   # 组合长音符在 PDF 字体里会丢，改用 CSS 上划线
core = markdown.markdown(body, extensions=["tables", "sane_lists", "attr_list"])
core = core.replace("<li>[ ] ", '<li class="todo">')

def row(k, v):
    if v is None: return ""
    v = "、".join(str(x) for x in v) if isinstance(v, list) else str(v)
    return f"<tr><th>{H.escape(k)}</th><td>{H.escape(v)}</td></tr>"

LBL = {"venue": "发表于", "year": "年份", "paper_type": "类型", "institute": "机构",
       "doi": "DOI", "status": "状态", "topics": "主题", "Annotator": "笔记作者", "date_created": "记于"}
rows = "".join(row(LBL.get(k, k), meta.get(k)) for k in
               ["venue", "year", "paper_type", "doi", "institute", "topics", "Annotator", "date_created"])

CSS = """
@page { size: A4; margin: 18mm 17mm 16mm 17mm;
  @bottom-center { content: counter(page); } }
:root { --ink:#111; --ink2:#454443; --muted:#7c7a76; --rule:#e2e0db; --accent:#8a3324; --surf:#fbfaf8; }
* { box-sizing: border-box; }
body { font-family:"Noto Serif CJK SC","Noto Serif",serif; color:var(--ink);
  font-size:10.2pt; line-height:1.66; margin:0; -webkit-font-smoothing:antialiased; }
h1 { font-family:"Noto Sans CJK SC",sans-serif; font-size:17pt; line-height:1.34; margin:0 0 2mm;
  letter-spacing:-.01em; }
.sub { font-family:"Noto Sans CJK SC",sans-serif; color:var(--muted); font-size:9.4pt; margin:0 0 5mm; }
table.meta { border-collapse:collapse; width:100%; margin:0 0 8mm; font-size:9.1pt;
  font-family:"Noto Sans CJK SC",sans-serif; background:var(--surf);
  border-top:1.4px solid var(--ink); border-bottom:1.4px solid var(--ink); }
table.meta th { text-align:left; color:var(--muted); font-weight:500; white-space:nowrap;
  width:22mm; padding:1.5mm 3mm 1.5mm 3mm; vertical-align:top; }
table.meta td { padding:1.5mm 3mm; color:var(--ink2); vertical-align:top; }
table.meta tr + tr th, table.meta tr + tr td { border-top:1px solid var(--rule); }
h2 { font-family:"Noto Sans CJK SC",sans-serif; font-size:12.2pt; font-weight:600; color:var(--accent);
  margin:7mm 0 2.2mm; padding-bottom:1.4mm; border-bottom:1px solid var(--rule);
  break-after:avoid; page-break-after:avoid; }
h2:first-of-type { margin-top:2mm; }
h3 { font-family:"Noto Sans CJK SC",sans-serif; font-size:10.4pt; font-weight:600; color:var(--ink);
  margin:5.5mm 0 1.5mm; break-after:avoid; page-break-after:avoid; }
p { margin:0 0 2.6mm; text-align:justify; }
ul,ol { margin:0 0 3mm; padding-left:5.6mm; }
li { margin:0 0 1.5mm; }
li > p { margin:0 0 1.2mm; }
li::marker { color:var(--muted); }
strong { font-weight:700; }
em { font-style:normal; color:var(--accent); }
code { font-family:"DejaVu Sans Mono",monospace; font-size:8.9pt; background:var(--surf);
  border:1px solid var(--rule); border-radius:2px; padding:0 1mm; }
.wl { font-family:"Noto Sans CJK SC",sans-serif; font-size:9.3pt; color:var(--accent);
  background:#f6f1ee; border-radius:2px; padding:0 1.2mm; }
.bar { text-decoration:overline; text-decoration-thickness:.7px; text-underline-offset:0; }
li.todo { list-style:none; margin-left:-4mm; }
li.todo::before { content:"☐"; color:var(--muted); margin-right:2mm; }
h2, h3, li, tr { break-inside:avoid; page-break-inside:avoid; }
p { orphans:2; widows:2; }
"""

doc = f"""<!doctype html><html lang="zh"><meta charset="utf-8">
<title>{H.escape(str(meta.get('title','')))}</title><style>{CSS}</style>
<body>
<h1>{H.escape(str(meta.get('title','')))}</h1>
<p class="sub">论文笔记 · Paper Note</p>
<table class="meta">{rows}</table>
{core}
</body></html>"""

tmp = os.path.splitext(out)[0] + ".html"
open(tmp, "w", encoding="utf-8").write(doc)
import glob
cands = sorted(glob.glob("/opt/pw-browsers/chromium*/chrome-linux*/chrome")) + \
        [p for p in ("/usr/bin/chromium", "/usr/bin/google-chrome") if os.path.exists(p)]
chrome = next(p for p in cands if os.path.isfile(p) and os.access(p, os.X_OK))
subprocess.run([chrome, "--headless", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
                f"--print-to-pdf={out}", "file://" + os.path.abspath(tmp)], check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
os.remove(tmp)
print("→", out)
