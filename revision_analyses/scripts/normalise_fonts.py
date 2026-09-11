#!/usr/bin/env python3
"""Clear the 'erewhon' font left on Table 2's runs so they inherit the document font.

erewhon is a LaTeX text font; it will not be installed on the typesetter's machine and
Word would substitute something arbitrary. The runs carry no other formatting intent, so
dropping the explicit face lets them take the document default. Text is untouched - the
script asserts it.

Run after revise_affiliations.py.
"""
import docx
from docx.oxml.ns import qn

STRIP = {"erewhon"}

for path in ("../02_Revised_Manuscript_clean.docx", "../03_Revised_Manuscript_tracked.docx"):
    doc = docx.Document(path)
    before = [p.text for p in doc.paragraphs]
    n = 0
    for p in doc.paragraphs:
        for r in p.runs:
            if r.font.name in STRIP:
                rpr = r._element.find(qn("w:rPr"))
                if rpr is not None:
                    for f in rpr.findall(qn("w:rFonts")):
                        rpr.remove(f); n += 1
    doc.save(path)
    after = [p.text for p in docx.Document(path).paragraphs]
    assert before == after, f"text changed in {path}"
    print(f"{path}: cleared {n} explicit font references, text unchanged")
