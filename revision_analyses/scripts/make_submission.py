#!/usr/bin/env python3
"""Assemble SUBMISSION/ - the files that go to the journal, under plain numbered names.

Everything here is a copy; figures/ stays the working output directory and analysis/ the
reproducible source. Run from build/.
"""
import os, shutil, subprocess

R = ".."
S = os.path.join(R, "SUBMISSION")
os.makedirs(S, exist_ok=True)

DOCS = [
    ("01_Response_to_Reviewers.docx", "01_Response_to_Reviewers.docx"),
    ("02_Revised_Manuscript_clean.docx", "02_Revised_Manuscript_clean.docx"),
    ("03_Revised_Manuscript_tracked.docx", "03_Revised_Manuscript_tracked.docx"),
    ("04_Supplementary_Methods.docx", "04_Supplementary_Methods.docx"),
    ("supplementary/Supplementary_Tables_S1-S8.xlsx", "05_Supplementary_Tables.xlsx"),
]
FIGS = [
    ("figures/Figure1_composite.pdf", "Figure1.pdf"),
    ("figures/Figure_2_manual-merge.pdf", "Figure2.pdf"),
    ("figures/Figure_3_manual-merge.pdf", "Figure3.pdf"),
    ("figures/SupplementaryFigure1_epitope.pdf", "SupplementaryFigure1.pdf"),
    ("figures/SupplementaryFigure2_RF2_metrics.pdf", "SupplementaryFigure2.pdf"),
    ("figures/SupplementaryFigure3_mAb_panel.pdf", "SupplementaryFigure3.pdf"),
    ("figures/SupplementaryFigure4_crossreactivity.pdf", "SupplementaryFigure4.pdf"),
    ("figures/SupplementaryFigure5_normal_tissue.pdf", "SupplementaryFigure5.pdf"),
]

# Figures whose artboard is far larger than the page get their page box scaled to the
# journal's 170 mm double-column width. This is a lossless transform of the whole page:
# vector stays vector and the raster pixels are untouched, so the effective resolution
# rises by the same factor. It fixes the trim size, not the relative size of any type -
# an element drawn too small for its panel is still too small afterwards.
SCALE_TO_MM = {"Figure2.pdf": 170.0, "Figure3.pdf": 170.0, "SupplementaryFigure1.pdf": 170.0,
               "SupplementaryFigure5.pdf": 170.0}

def copy_scaled(src, dst_path, target_mm):
    from pypdf import PdfReader, PdfWriter
    r = PdfReader(src); page = r.pages[0]
    factor = target_mm / (float(page.mediabox.width) / 72 * 25.4)
    page.scale_by(factor)
    w = PdfWriter(); w.add_page(page); w.write(dst_path)
    return factor

for src, dst in DOCS + FIGS:
    s = os.path.join(R, src)
    if not os.path.exists(s):
        raise SystemExit(f"MISSING: {src}")
    d = os.path.join(S, dst)
    if dst in SCALE_TO_MM:
        f = copy_scaled(s, d, SCALE_TO_MM[dst])
        print(f"  scaled {dst} to {SCALE_TO_MM[dst]:.0f} mm wide (x{f:.4f})")
    else:
        shutil.copyfile(s, d)

print(f"SUBMISSION/  ({len(DOCS) + len(FIGS)} files)\n")
for _, dst in DOCS:
    print(f"  {dst:<34}{os.path.getsize(os.path.join(S, dst))/1024:>8.0f} KB")
print()
for src, dst in FIGS:
    p = os.path.join(S, dst)
    info = subprocess.run(["pdfinfo", p], capture_output=True, text=True).stdout
    size = next((l.split(":", 1)[1].strip() for l in info.splitlines()
                 if l.startswith("Page size")), "?")
    pages = next((l.split(":", 1)[1].strip() for l in info.splitlines()
                  if l.startswith("Pages")), "?")
    mm = "?"
    try:
        w, h = [float(x) for x in size.split()[0:3:2]]
        mm = f"{w/72*25.4:.0f} x {h/72*25.4:.0f} mm"
    except Exception:
        pass
    print(f"  {dst:<34}{os.path.getsize(p)/1024:>8.0f} KB  {pages}p  {mm:<16} <- {os.path.basename(src)}")
