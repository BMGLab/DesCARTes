#!/usr/bin/env python3
"""
Revision analysis 6 -- Pan-normal-tissue CLDN4 protein expression at cell-type
resolution.  Addresses Reviewer 1 major point 1 and Reviewer 2 major point 3.

Source: Human Protein Atlas v22.0 immunohistochemistry, as bundled with
HPAanalyze 1.24.0 (hpa_histology_data$normal_tissue), antibody reliability
"Enhanced".  Extracted with scripts/R6_extract_hpa.R.

Outputs: output/TableS_CLDN4_normal_tissue.csv
         output/TableS_CLDN4_pan_cancer.csv
         output/hpa_normal_tissue_summary.txt
"""
import numpy as np
import pandas as pd

ORDER = ["Not detected", "Low", "Medium", "High"]
NUM = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}

n = pd.read_csv("output/cldn4_hpa_normal_all.csv")
n["level"] = pd.Categorical(n.level, categories=ORDER, ordered=True)
n["score"] = n.level.map(NUM)

# organ-system grouping for the safety discussion
SYSTEM = {
    "appendix": "Gastrointestinal", "colon": "Gastrointestinal", "duodenum": "Gastrointestinal",
    "rectum": "Gastrointestinal", "small intestine": "Gastrointestinal",
    "stomach 1": "Gastrointestinal", "stomach 2": "Gastrointestinal",
    "esophagus": "Gastrointestinal", "oral mucosa": "Gastrointestinal",
    "salivary gland": "Gastrointestinal", "liver": "Hepatobiliary/pancreas",
    "gallbladder": "Hepatobiliary/pancreas", "pancreas": "Hepatobiliary/pancreas",
    "kidney": "Renal/urinary", "urinary bladder": "Renal/urinary",
    "lung": "Respiratory", "bronchus": "Respiratory",
    "breast": "Reproductive/breast", "cervix": "Reproductive/breast",
    "endometrium 1": "Reproductive/breast", "endometrium 2": "Reproductive/breast",
    "fallopian tube": "Reproductive/breast", "ovary": "Reproductive/breast",
    "placenta": "Reproductive/breast", "prostate": "Reproductive/breast",
    "seminal vesicle": "Reproductive/breast", "epididymis": "Reproductive/breast",
    "testis": "Reproductive/breast", "vagina": "Reproductive/breast",
    "skin 1": "Skin", "skin 2": "Skin",
    "caudate": "CNS", "cerebellum": "CNS", "cerebral cortex": "CNS", "hippocampus": "CNS",
    "heart muscle": "Cardiac/muscle", "skeletal muscle": "Cardiac/muscle",
    "smooth muscle": "Cardiac/muscle",
    "bone marrow": "Haematopoietic/lymphoid", "lymph node": "Haematopoietic/lymphoid",
    "spleen": "Haematopoietic/lymphoid", "tonsil": "Haematopoietic/lymphoid",
    "adrenal gland": "Endocrine", "parathyroid gland": "Endocrine",
    "thyroid gland": "Endocrine",
    "adipose tissue": "Connective/other", "soft tissue 1": "Connective/other",
    "soft tissue 2": "Connective/other",
}
n["organ_system"] = n.tissue.map(SYSTEM).fillna("Other")
n = n.sort_values(["organ_system", "tissue", "score"], ascending=[True, True, False])
n[["organ_system", "tissue", "cell_type", "level", "reliability"]].to_csv(
    "output/TableS_CLDN4_normal_tissue.csv", index=False)

lines = []
def emit(s=""):
    print(s); lines.append(s)

emit("=" * 86)
emit("CLDN4 PROTEIN EXPRESSION ACROSS NORMAL HUMAN TISSUES (HPA v22.0, IHC, Enhanced)")
emit("=" * 86)
emit(f"{n.tissue.nunique()} tissues, {len(n)} tissue x cell-type entries")
emit("")
emit("Staining level distribution:")
for lv in ORDER[::-1]:
    c = int((n.level == lv).sum())
    emit(f"  {lv:14s} {c:4d} ({c/len(n)*100:4.1f}%)")
emit("")

pos = n[n.score >= 2]
emit(f"Cell types with MEDIUM or HIGH CLDN4 staining: {len(pos)} "
     f"across {pos.tissue.nunique()} tissues")
emit("")
emit("--- HIGH staining ---")
for _, r in n[n.level == "High"].iterrows():
    emit(f"  {r.organ_system:26s} {r.tissue:18s} {r.cell_type}")
emit("")
emit("--- MEDIUM staining ---")
for _, r in n[n.level == "Medium"].iterrows():
    emit(f"  {r.organ_system:26s} {r.tissue:18s} {r.cell_type}")
emit("")

emit("--- Per organ system: maximum CLDN4 level observed ---")
summ = (n.groupby("organ_system")
          .agg(n_entries=("level", "size"),
               max_level=("score", "max"),
               n_medium_or_high=("score", lambda s: int((s >= 2).sum())))
          .sort_values("max_level", ascending=False))
summ["max_level"] = summ.max_level.map({v: k for k, v in NUM.items()})
emit(summ.to_string())
emit("")
emit("Tissues where CLDN4 is NOT detected in any assessed cell type:")
neg = n.groupby("tissue").score.max()
emit("  " + ", ".join(sorted(neg[neg == 0].index)))
emit("")
emit("NOTE: HPA reports total cellular immunoreactivity and does not distinguish")
emit("membrane-accessible antigen from junction-sequestered or intracellular protein.")
emit("These data therefore bound the on-target/off-tumour risk from above; they do")
emit("not establish that CLDN4 is accessible to a CAR on these cells.")

# --- pan-cancer -----------------------------------------------------------------
c = pd.read_csv("output/cldn4_hpa_pathology_all.csv")
c.to_csv("output/TableS_CLDN4_pan_cancer.csv", index=False)
emit("")
emit("=" * 86)
emit(f"CLDN4 IHC ACROSS {c.cancer.nunique()} CANCER TYPES (HPA v22.0 pathology)")
emit("=" * 86)
cols = [x for x in ["high", "medium", "low", "not_detected"] if x in c.columns]
if cols:
    cc = c.set_index("cancer")[cols].fillna(0)
    cc["n_patients"] = cc.sum(axis=1)
    cc["pct_med_high"] = (cc[[x for x in ["high", "medium"] if x in cc]].sum(axis=1)
                          / cc.n_patients * 100).round(0)
    emit(cc.sort_values("pct_med_high", ascending=False).to_string())

open("output/hpa_normal_tissue_summary.txt", "w").write("\n".join(lines) + "\n")
print("\nSaved -> output/TableS_CLDN4_normal_tissue.csv, TableS_CLDN4_pan_cancer.csv")
