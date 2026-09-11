#!/usr/bin/env python3
"""Update the running head, which still carried the old title.

header2.xml reads "De novo design of CLDN4 binders for CAR-T therapy" - the framing the
revision moved away from, since the paper is now an in silico prioritization study. Run
after revise_affiliations.py.
"""
import re, shutil, zipfile

OLD = "De novo design of CLDN4 binders for CAR-T therapy"
NEW = "STRIVE: in silico design and off-target triage of CLDN4 scFv binders"

for path in ("../02_Revised_Manuscript_clean.docx", "../03_Revised_Manuscript_tracked.docx"):
    tmp = path + ".tmp"
    hits = 0
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if re.match(r"word/header\d+\.xml", item.filename):
                text = data.decode("utf-8")
                if OLD in text:
                    text = text.replace(OLD, NEW)
                    hits += 1
                    data = text.encode("utf-8")
            zout.writestr(item, data)
    shutil.move(tmp, path)
    print(f"{path}: running head updated in {hits} header part(s)")
