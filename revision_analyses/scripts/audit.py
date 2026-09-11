#!/usr/bin/env python3
"""Pre-submission audit. Run from build/:  python3 audit.py"""
import os, re, subprocess, sys
import docx
from docx.oxml.ns import qn

R = ".."
CLEAN = f"{R}/02_Revised_Manuscript_clean.docx"
TRACKED = f"{R}/03_Revised_Manuscript_tracked.docx"
SUPP = f"{R}/04_Supplementary_Methods.md"
SUB = f"{R}/SUBMISSION"
fails, warns = [], []

def fail(m): fails.append(m); print(f"  FAIL  {m}")
def warn(m): warns.append(m); print(f"  warn  {m}")
def ok(m):   print(f"  ok    {m}")

ps = [p.text.strip() for p in docx.Document(CLEAN).paragraphs]
iref = next(i for i, t in enumerate(ps) if re.match(r'^REFERENCES', t, re.I))
body = " ".join(ps[:iref])

print("\n== word counts ==")
marks = {}
for i, t in enumerate(ps):
    for k, pat in [("intro", r'^1\.\s*INTRODUCTION'), ("meth", r'^2\.\s*MATERIALS'),
                   ("res", r'^3\.\s*RESULTS'), ("disc", r'^4\.\s*DISCUSSION'),
                   ("end", r'^(ETHICS|CONFLICT OF INTEREST|ACKNOWLEDG|LIST OF ABBREV|CONSENT FOR)')]:
        if k not in marks and re.match(pat, t, re.I): marks[k] = i
wc = lambda a, b: sum(len(t.split()) for t in ps[a:b])
main = wc(marks['intro'], marks['end'])
(ok if 4000 <= main <= 6000 else fail)(f"main text {main} words (limit 4,000-6,000)")
abstract = next((len(c.text.split()) - 1 for tb in docx.Document(CLEAN).tables
                 for row in tb.rows for c in row.cells if c.text.strip().startswith("Abstract")), None)
(ok if abstract and abstract <= 250 else fail)(f"abstract {abstract} words (limit 250)")

print("\n== structure ==")
hs = [t for t in ps if re.match(r'^\d+(\.\d+)*\.\s+\S', t) and len(t.split()) <= 14]
nums = [h.split()[0].rstrip(".") for h in hs]
ok(f"headings: {' '.join(nums)}")
for i, t in enumerate(ps):
    if not (re.match(r'^\d+(\.\d+)*\.\s', t) and len(t.split()) <= 14): continue
    nxt = next((ps[j] for j in range(i + 1, len(ps)) if ps[j]), "")
    if re.match(r'^\d+(\.\d+)*\.\s', nxt) and len(nxt.split()) <= 14 \
       and nxt.split()[0].count(".") <= t.split()[0].count("."):
        fail(f"empty section: {t}")
else:
    ok("no section heading is left without text")

print("\n== run formatting ==")
_d = docx.Document(CLEAN)
_head = lambda t: bool(re.match(r'^\d+(\.\d+)*\.\s+\S', t.strip())) and len(t.split()) <= 14
_runs = lambda p: [(bool(r.bold) or bool(p.style.font.bold), (r.text or ""))
                   for r in p.runs if (r.text or "").strip()]
_bad_h, _bad_b, _bad_l = [], [], []
for i, p in enumerate(_d.paragraphs):
    t = p.text.strip()
    rs = _runs(p)
    if not t or not rs or i < 12 or t.isupper():
        continue
    if _head(t):
        if not all(b for b, _ in rs): _bad_h.append(t[:50])
    elif re.match(r'^((?:Supplementary )?(?:Figure|Table)\s+S?\d+)', t):
        if all(b for b, _ in rs): _bad_l.append(t[:50])
    elif all(b for b, _ in rs) and len(t.split()) > 8:
        _bad_b.append(f"[{i}] {t[:55]}")
(ok if not _bad_h else fail)(f"numbered headings bold: {len(_bad_h) or 'all'}"
                             + (f" -- {_bad_h}" if _bad_h else ""))
(ok if not _bad_b else fail)(f"body paragraphs wrongly bold: {len(_bad_b) or 'none'}"
                             + (f" -- {_bad_b}" if _bad_b else ""))
(ok if not _bad_l else fail)(f"legends wholly bold: {len(_bad_l) or 'none'}"
                             + (f" -- {_bad_l}" if _bad_l else ""))

print("\n== leftover objects ==")
_M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
_orphans = []
for i, p in enumerate(docx.Document(CLEAN).paragraphs):
    mt = "".join(t.text or "" for t in p._p.iter(f"{{{_M}}}t"))
    if mt and p.text.strip():
        _orphans.append(f"[{i}] {mt[:45]!r} in {p.text.strip()[:40]!r}")
(ok if not _orphans else fail)(
    f"equation objects stranded in rewritten paragraphs: {len(_orphans) or 'none'}")
for o in _orphans: print(f"        {o}")
import zipfile as _zf, re as _re
_z = _zf.ZipFile(CLEAN)
_heads = " ".join("".join(_re.findall(r'<w:t[^>]*>([^<]*)</w:t>', _z.read(n).decode("utf8", "ignore")))
                  for n in _z.namelist() if _re.match(r'word/header\d+\.xml', n))
(ok if "De novo design of CLDN4 binders for CAR-T" not in _heads else fail)("running head matches the revised title")
_body = _z.read("word/document.xml").decode("utf8", "ignore")
_stale = [v for v in ("1.67", "2.58", "1.33 ", "5.25", "1.82", "2.86") if v in _body]
(ok if not _stale else fail)(f"superseded p-values in the file: {_stale or 'none'}")

print("\n== references ==")
cited = set()
for m in re.finditer(r'\[([\d,\s]+)\]', body):
    cited |= {int(x) for x in m.group(1).split(",") if x.strip().isdigit()}
listed = {int(re.match(r'^\[?(\d+)', t).group(1)) for t in ps[iref:] if re.match(r'^\[?\d+[\].]', t)}
(ok if not listed - cited else fail)(f"{len(listed)} references; uncited: {sorted(listed - cited) or 'none'}")
(ok if not {c for c in cited - listed if c > 1} else fail)(
    f"citations with no entry: {sorted(c for c in cited - listed if c > 1) or 'none'}")
(ok if len(listed) >= 75 else warn)(f"{len(listed)} references (journal minimum 75)")

print("\n== figures and tables ==")
legends = {}
for t in ps:
    m = re.match(r'^((?:Supplementary )?(?:Figure|Table) S?\d+)[.\s]', t)
    if m: legends[m.group(1)] = t
calls = set()
for t in ps[:iref]:
    if re.match(r'^(Supplementary )?(Figure|Table) S?\d+[.\s]', t): continue
    calls |= {m.group(0) for m in re.finditer(r'(?:Supplementary )?(?:Figure|Table)\s+S?\d+', t)}
norm = lambda s: s.replace("Supplementary Table S", "Supplementary Table S")
for c in sorted(calls):
    if c not in legends and not c.startswith("Supplementary Table"):
        fail(f"cited with no legend: {c}")
for l in sorted(legends):
    if l not in calls: fail(f"legend never cited: {l}")
ok(f"{len(legends)} legends, all cited; {len(calls)} distinct callouts")

# panel letters: every letter in a legend must be called out, and vice versa
for name, leg in legends.items():
    if not name.startswith(("Figure", "Supplementary Figure")): continue
    inleg = set(re.findall(r'\(([A-F])(?:,\s*([A-F]))?\)', leg))
    inleg = {x for tup in inleg for x in tup if x}
    n = name.replace("Supplementary ", "S")
    called = set(re.findall(rf'{re.escape(n)}([A-F])', body)) | \
             {m for m in re.findall(r'\b\d([A-F])\b', body)} if inleg else set()
    called = set(re.findall(rf'{re.escape(name)}([A-F])', body))
    for extra in re.finditer(rf'{re.escape(name)}[A-F](?:,\s*\d?([A-F]))+', body):
        called |= {g for g in extra.groups() if g}
    missing = inleg - called
    if missing and inleg:
        warn(f"{name}: panels {sorted(missing)} in the legend but not called out in the text")

print("\n== supplementary ==")
smeth = {l.split(".")[0][3:] for l in open(SUPP, encoding="utf-8") if l.startswith("## S")}
scited = {m.group(1) for t in ps for m in re.finditer(r'Supplementary Methods\s+(S\d)', t)}
(ok if not scited - smeth else fail)(f"Supplementary Methods sections {sorted(smeth)}; dangling: {sorted(scited - smeth) or 'none'}")
(ok if not smeth - scited else warn)(f"uncited Supplementary Methods sections: {sorted(smeth - scited) or 'none'}")

print("\n== submission files ==")
want = ["01_Response_to_Reviewers.docx", "02_Revised_Manuscript_clean.docx",
        "03_Revised_Manuscript_tracked.docx", "04_Supplementary_Methods.docx",
        "05_Supplementary_Tables.xlsx", "Figure1.pdf", "Figure2.pdf", "Figure3.pdf"] + \
       [f"SupplementaryFigure{i}.pdf" for i in range(1, 6)]
for f in want:
    p = os.path.join(SUB, f)
    if not os.path.exists(p): fail(f"missing from SUBMISSION/: {f}"); continue
    if f.endswith(".pdf"):
        o = subprocess.run(["pdfinfo", p], capture_output=True, text=True).stdout
        sz = [l for l in o.splitlines() if l.startswith("Page size")][0].split(":")[1].strip()
        w, h = float(sz.split()[0]), float(sz.split()[2])
        mm = (w / 72 * 25.4, h / 72 * 25.4)
        tag = "ok   " if mm[0] <= 180 else "warn "
        if mm[0] > 180: warns.append(f"{f} is {mm[0]:.0f} mm wide (double column is 170 mm)")
        print(f"  {tag} {f:<30}{mm[0]:>6.0f} x {mm[1]:>4.0f} mm")
    else:
        ok(f"{f} ({os.path.getsize(p)/1024:.0f} KB)")

print("\n== tracked changes ==")
d = docx.Document(TRACKED)
ins = len(d.element.body.findall(".//" + qn("w:ins")))
dele = len(d.element.body.findall(".//" + qn("w:del")))
def accepted(p):
    out = []
    for n in p.iter():
        if n.tag == qn("w:t"):
            a, skip = n, False
            while a is not None:
                if a.tag == qn("w:del"): skip = True; break
                a = a.getparent()
            if not skip: out.append(n.text or "")
    return "".join(out).strip()
acc = [t for t in (accepted(p._p) for p in d.paragraphs) if t and not t.startswith("[Revision note:")]
cl = [t for t in ps if t]
diff = sum(1 for a, b in zip(acc, cl) if a.replace("\t", "") != b.replace("\t", ""))
(ok if len(acc) == len(cl) and diff == 0 else fail)(
    f"{ins} insertions, {dele} deletions; accept-all reproduces the clean file "
    f"({len(acc)} vs {len(cl)} paragraphs, {diff} mismatches)")

print("\n" + "=" * 62)
print(f"{len(fails)} failures, {len(warns)} warnings")
for w in warns: print(f"  warn  {w}")
sys.exit(1 if fails else 0)
