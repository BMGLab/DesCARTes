#!/usr/bin/env python3
"""Build the clean and tracked-changes revised manuscripts from the submitted DOCX."""
import copy, datetime, re, shutil, sys
from docx import Document
from docx.oxml.ns import qn, nsmap
from docx.oxml import OxmlElement

SRC = "/home/biolab/Projects/DesCARTes_wd/docs/Current_Gene_Therapy/CGT_STRIVE_final.docx"
CLEAN = "../02_Revised_Manuscript_clean.docx"
TRACKED = "../03_Revised_Manuscript_tracked.docx"
AUTHOR = "Author revision"
STAMP = datetime.datetime(2026, 9, 8, 12, 0, 0).strftime("%Y-%m-%dT%H:%M:%SZ")

ns = {}
for f in ["edits_part1.py", "edits_part2.py", "edits_part3.py", "edits_part4.py",
          "edits_part5.py", "edits_part6.py", "edits_part7.py", "edits_part8.py", "edits_part9.py", "edits_part10.py", "edits_part11.py", "edits_part12.py", "edits_part13.py", "edits_part14.py", "edits_part15.py", "edits_part16.py", "edits_part17.py", "edits_part18.py", "edits_part19.py", "edits_part20.py", "edits_part21.py", "edits_part22.py", "edits_part23.py"]:
    exec(open(f).read(), ns)
EDITS = {}
for name in ["EDITS", "EDITS2", "EDITS3", "EDITS4", "EDITS5", "EDITS6", "EDITS7",
             "EDITS8", "EDITS9", "EDITS10", "EDITS11", "EDITS12", "EDITS13", "EDITS14", "EDITS15", "EDITS16", "EDITS17", "EDITS18", "EDITS19", "EDITS20", "EDITS21", "EDITS22", "EDITS23"]:
    EDITS.update(ns[name])
ABSTRACT = ns["ABSTRACT"]
print(f"{len(EDITS)} paragraph edits loaded")

_uid = [1000]
def uid():
    _uid[0] += 1
    return str(_uid[0])

def para_text(p):
    return "".join(n.text or "" for n in p.findall(".//" + qn("w:t")))

OMML = "http://schemas.openxmlformats.org/officeDocument/2006/math"

def clear_runs(p):
    """Strip a paragraph's content before writing the replacement.

    This must include the OMML equation objects (<m:oMath>, <m:oMathPara>), not only the
    <w:r> runs. They are siblings of the runs, so removing runs alone left the original
    inline equations behind and the new text was appended after them - which is how
    superseded p-values (p = 1.67 x 10-31), an obsolete thermostat constant and the old
    R_specificity nomenclature survived into paragraphs that had been rewritten. The
    stand-alone display equations in the Methods live in their own paragraphs, which carry
    no text and are never edit targets, so they are untouched.
    """
    for child in list(p):
        if child.tag in (qn("w:r"), qn("w:ins"), qn("w:del"), qn("w:hyperlink"),
                         f"{{{OMML}}}oMath", f"{{{OMML}}}oMathPara"):
            p.remove(child)

def rpr_of(p):
    """Run properties of the run carrying the most text, not of the first run.

    The first run is often a short bold lead-in - "[12] " in a reference, "Figure 3." in a
    legend - so copying it turned whole paragraphs bold. The dominant run is the body
    formatting; bold is then set per segment by segments() below.
    """
    best, best_len = None, -1
    for r in p.findall(qn("w:r")):
        n = len("".join(t.text or "" for t in r.findall(qn("w:t"))).strip())
        if n > best_len:
            rpr = r.find(qn("w:rPr"))
            best, best_len = (copy.deepcopy(rpr) if rpr is not None else None), n
    return best

def with_bold(rpr, bold):
    """Copy of rpr with <w:b>/<w:bCs> forced on or off."""
    out = copy.deepcopy(rpr) if rpr is not None else OxmlElement("w:rPr")
    for tag in ("w:b", "w:bCs"):
        for el in out.findall(qn(tag)):
            out.remove(el)
        if bold:
            el = OxmlElement(tag); el.set(qn("w:val"), "1"); out.insert(0, el)
    return out

HEADING_RE = re.compile(r"^\d+(\.\d+)*\.\s+\S")
REF_RE = re.compile(r"^(\[\d+\]\s)")
LEGEND_RE = re.compile(r"^((?:Supplementary )?(?:Figure|Table)\s+S?\d+\.?\s)")

def segments(text):
    """Split replacement text into (chunk, bold) runs matching the journal's own style:
    numbered headings bold throughout; a reference's "[n] " marker bold; a legend's
    "Figure N." plus its title sentence bold, and each "(A)" panel marker bold; everything
    else regular."""
    t = text.strip()
    if HEADING_RE.match(t) and len(t.split()) <= 14:
        return [(text, True)]
    m = REF_RE.match(t)
    if m:
        return [(m.group(1), True), (text[len(m.group(1)):], False)]
    m = LEGEND_RE.match(t)
    if m:
        # only the "Figure N." label is bold; the rest of the legend, title sentence and
        # panel markers included, is regular
        head = m.group(1)
        return [(head, True), (text[len(head):], False)]
    return [(text, False)]

def make_run(text, rpr=None, deleted=False):
    r = OxmlElement("w:r")
    if rpr is not None:
        r.append(copy.deepcopy(rpr))
    t = OxmlElement("w:delText" if deleted else "w:t")
    t.set(qn("xml:space"), "preserve")
    t.text = text
    r.append(t)
    return r

def wrap(tag, children):
    el = OxmlElement(tag)
    el.set(qn("w:id"), uid())
    el.set(qn("w:author"), AUTHOR)
    el.set(qn("w:date"), STAMP)
    for c in children:
        el.append(c)
    return el

def set_clean(p, text):
    rpr = rpr_of(p)
    clear_runs(p)
    for chunk, bold in segments(text):
        p.append(make_run(chunk, with_bold(rpr, bold)))

def set_tracked(p, old_text, new_text):
    """Mark the old text deleted and the new text inserted, in place."""
    rpr = rpr_of(p)
    clear_runs(p)
    if old_text.strip():
        p.append(wrap("w:del", [make_run(old_text, rpr, deleted=True)]))
    if new_text.strip():
        p.append(wrap("w:ins", [make_run(chunk, with_bold(rpr, bold))
                               for chunk, bold in segments(new_text)]))

def insert_para_after(p, text, tracked):
    new_p = OxmlElement("w:p")
    ppr = p.find(qn("w:pPr"))
    if ppr is not None:
        new_p.append(copy.deepcopy(ppr))
    rpr = rpr_of(p)
    # the anchor may be a heading; the inserted paragraph takes its own formatting
    runs = [make_run(chunk, with_bold(rpr, bold)) for chunk, bold in segments(text)]
    if tracked:
        new_p.append(wrap("w:ins", runs))
        # mark the paragraph mark itself as inserted; in CT_ParaRPr the <w:ins>
        # element must be the FIRST child of <w:rPr>, and <w:rPr> the last child
        # of <w:pPr>, or Word/LibreOffice reject the document.
        ppr2 = new_p.find(qn("w:pPr"))
        if ppr2 is None:
            ppr2 = OxmlElement("w:pPr")
            new_p.insert(0, ppr2)
        rpr_mark = ppr2.find(qn("w:rPr"))
        if rpr_mark is None:
            rpr_mark = OxmlElement("w:rPr")
            ppr2.append(rpr_mark)
        for stale in rpr_mark.findall(qn("w:ins")):
            rpr_mark.remove(stale)
        rpr_mark.insert(0, wrap("w:ins", []))
    else:
        for r in runs:
            new_p.append(r)
        ppr2 = new_p.find(qn("w:pPr"))
        if ppr2 is not None:
            rpr_mark = ppr2.find(qn("w:rPr"))
            if rpr_mark is not None:
                for stale in rpr_mark.findall(qn("w:ins")):
                    rpr_mark.remove(stale)
    p.addnext(new_p)
    return new_p

def build(dest, tracked):
    shutil.copy(SRC, dest)
    doc = Document(dest)
    applied, missing = set(), []

    # ---- abstract lives inside the header table -------------------------------
    for tbl in doc.tables:
        for row in tbl.rows:
            for cell in row.cells:
                for p in cell.paragraphs:
                    t = para_text(p._p)
                    if t.strip().startswith("Abstract Background:") or \
                       (t.strip().startswith("Background: Chimeric antigen receptor")):
                        new = ("Abstract " + ABSTRACT) if t.strip().startswith("Abstract") else ABSTRACT
                        if tracked:
                            set_tracked(p._p, t, new)
                        else:
                            set_clean(p._p, new)
                        applied.add("__abstract__")

    # ---- body paragraphs -------------------------------------------------------
    for p in list(doc.paragraphs):
        t = para_text(p._p)
        ts = t.strip()
        if not ts:
            continue
        for key, val in EDITS.items():
            if key in applied:
                continue
            if ts.startswith(key[:60]) or key[:60] in ts[:200]:
                if val is None:                       # delete the paragraph
                    if tracked:
                        set_tracked(p._p, t, "")
                        ppr = p._p.find(qn("w:pPr"))
                        if ppr is None:
                            ppr = OxmlElement("w:pPr"); p._p.insert(0, ppr)
                        rpr_mark = ppr.find(qn("w:rPr"))
                        if rpr_mark is None:
                            rpr_mark = OxmlElement("w:rPr"); ppr.append(rpr_mark)
                        for stale in rpr_mark.findall(qn("w:del")):
                            rpr_mark.remove(stale)
                        rpr_mark.insert(0, wrap("w:del", []))
                    else:
                        p._p.getparent().remove(p._p)
                elif isinstance(val, list):
                    if tracked:
                        set_tracked(p._p, t, val[0])
                    else:
                        set_clean(p._p, val[0])
                    anchor = p._p
                    for extra in val[1:]:
                        anchor = insert_para_after(anchor, extra, tracked)
                else:
                    if tracked:
                        set_tracked(p._p, t, val)
                    else:
                        set_clean(p._p, val)
                applied.add(key)
                break

    # Reviewer 2 asked for the blank final page to be removed: drop trailing empty
    # paragraphs so the document ends on the last line of content.
    body = doc.element.body
    removed = 0
    for child in reversed(list(body.iterchildren())):
        if child.tag == qn("w:sectPr"):
            continue
        if child.tag == qn("w:p") and not "".join(
                n.text or "" for n in child.findall(".//" + qn("w:t"))).strip():
            body.remove(child); removed += 1
        else:
            break
    if removed:
        print(f"  removed {removed} trailing empty paragraph(s)")

    missing = [k for k in EDITS if k not in applied]
    doc.save(dest)
    return applied, missing

for dest, tracked in [(CLEAN, False), (TRACKED, True)]:
    applied, missing = build(dest, tracked)
    print(f"\n{dest}: applied {len(applied)}/{len(EDITS)+1}")
    if missing:
        print("  NOT MATCHED:")
        for m in missing:
            print("   -", m[:80])
