"""build_pdf.py - DRAFT_EN.md -> HTML (marked + KaTeX) -> PDF via headless Edge/Chrome.

Usage: python -I build_pdf.py DRAFT_EN.md DRAFT_EN.pdf
Needs network access to cdn.jsdelivr.net at render time (marked, KaTeX).
HTML comments in the markdown (source notes) are not rendered.
"""
import json
import os
import subprocess
import sys
import tempfile

BROWSERS = [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
            r"C:\Program Files\Google\Chrome\Application\chrome.exe"]

TEMPLATE = """<!doctype html><html><head><meta charset="utf-8">
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script src="https://cdn.jsdelivr.net/npm/marked@12.0.2/marked.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"></script>
<style>
@page { size: A4; margin: 22mm 20mm; }
body { font-family: "Times New Roman", serif; font-size: 11pt; line-height: 1.45; color: #000; }
h1 { font-size: 17pt; text-align: center; } h2 { font-size: 13.5pt; margin-top: 1.4em; } h3 { font-size: 12pt; }
table { border-collapse: collapse; font-size: 9pt; margin: 0.6em 0; width: 100%; }
th, td { border: 1px solid #888; padding: 3px 5px; vertical-align: top; }
code { font-family: Consolas, monospace; font-size: 8.5pt; word-break: break-all; }
pre { white-space: pre-wrap; font-size: 8.5pt; }
blockquote { margin: 0.5em 1em; color: #222; }
</style></head><body><div id="c"></div>
<script>
const src = __SRC__;
// protect math from markdown processing
const store = [];
let s = src.replace(/\\$\\$[\\s\\S]+?\\$\\$|\\$[^$\\n]+?\\$/g, m => { store.push(m); return "@@MATH" + (store.length - 1) + "@@"; });
let html = marked.parse(s);
html = html.replace(/@@MATH(\\d+)@@/g, (_, i) => store[+i].replace(/&/g, "&amp;").replace(/</g, "&lt;"));
document.getElementById("c").innerHTML = html;
renderMathInElement(document.body, {delimiters: [{left: "$$", right: "$$", display: true}, {left: "$", right: "$", display: false}], throwOnError: false});
document.title = "done";
</script></body></html>"""


def main(md_path, pdf_path):
    with open(md_path, encoding="utf-8") as f:
        src = f.read()
    import re
    # Only absolute http(s) links are links; "[x](y)" in formulas must stay literal text.
    src = re.sub(r"\]\((?!https?://)", "]\\(", src)
    html = TEMPLATE.replace("__SRC__", json.dumps(src))
    tmp = tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8")
    tmp.write(html)
    tmp.close()
    browser = next(b for b in BROWSERS if os.path.exists(b))
    out = os.path.abspath(pdf_path)
    prof = tempfile.mkdtemp(prefix="pdfprof_")
    cmd = [browser, "--headless=new", "--disable-gpu", "--user-data-dir=" + prof, "--no-first-run", "--no-pdf-header-footer",
           "--virtual-time-budget=20000", "--print-to-pdf=" + out, "file:///" + tmp.name.replace("\\", "/")]
    subprocess.run(cmd, check=True, timeout=180)
    import time
    for _ in range(60):
        if os.path.exists(out) and os.path.getsize(out) > 0:
            break
        time.sleep(1)
    os.unlink(tmp.name)
    print("wrote", out, os.path.getsize(out), "bytes")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
