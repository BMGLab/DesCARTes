#!/usr/bin/env python3
"""
Figure 1A, rebuilt from UCSC Xena (Toil recompute of TCGA/GTEx, dataset
TcgaTargetGtex_rsem_gene_tpm, probe ENSG00000189143.9). This panel is now used in the
manuscript, both standalone and as panel A of the composite built by
R11_figure1_composite.py.

It replaces the originally submitted Figure 1A, which it did not reproduce:

    group        manuscript              Xena (log2 TPM+1)     Xena median TPM
    normal lung  5.34   (n = 50)         6.90  (n = 109)       119
    LUAD        13.05   (n = 510)        8.14  (n = 513)       280
    LUSC        11.97   (n = 484)        7.13  (n = 498)       139

The tumour sample counts agreed closely, so the cohorts were probably the same; the
normal set and the expression scale differed. The originally reported fold-changes
(207x LUAD, 99x LUSC) followed arithmetically from the original medians, but Xena TPM
gives 2.4x and 1.2x, and the single-cell analysis in Figure 1B gives roughly 2x between
normal epithelium and stage IV. A LUAD median of 13.05 on a log2(TPM+1) scale would
imply ~8,500 TPM, which would make CLDN4 one of the most abundant transcripts in the
cell. The Xena values are used instead, and the Results have been updated accordingly:
LUSC is not significantly different from normal lung after correction.

Outputs: figures/Figure1A_TCGA_Xena.{pdf,png}
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from statsmodels.stats.multitest import multipletests

BLUE=["#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7","#3987e5",
      "#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b"]
CMAP=LinearSegmentedColormap.from_list("seq_blue",BLUE)
SURF,INK,INK2="#fcfcfb","#0b0b0b","#52514e"
plt.rcParams.update({"font.size":8,"axes.edgecolor":"#c9c8c4","axes.labelcolor":INK2,
                     "xtick.color":INK2,"ytick.color":INK2,"figure.facecolor":SURF,
                     "axes.facecolor":SURF,"savefig.facecolor":SURF})

m=pd.read_csv("output/xena_cldn4.csv")
m["expr"]=np.log2(np.maximum(2**m.CLDN4.astype(float)-0.001,0)+1)
sel={"Normal lung":  m[(m._study=="TCGA")&(m._sample_type=="Solid Tissue Normal")
                       &(m.detailed_category.astype(str).str.contains("Lung",na=False))],
     "LUAD":         m[(m._study=="TCGA")&(m.detailed_category=="Lung Adenocarcinoma")
                       &(m._sample_type=="Primary Tumor")],
     "LUSC":         m[(m._study=="TCGA")&(m.detailed_category=="Lung Squamous Cell Carcinoma")
                       &(m._sample_type=="Primary Tumor")]}
order=list(sel)
data=[sel[k].expr.values for k in order]
# widest comparison last so its bracket sits highest and labels do not crowd
tests=[("Normal lung","LUAD"),("LUAD","LUSC"),("Normal lung","LUSC")]
p=[stats.mannwhitneyu(sel[a].expr,sel[b].expr,alternative="two-sided").pvalue for a,b in tests]
padj=multipletests(p,method="fdr_bh")[1]

fig,ax=plt.subplots(figsize=(4.0,3.4))
bp=ax.boxplot(data,patch_artist=True,widths=0.55,showfliers=False,
              medianprops=dict(color=INK,lw=1.3),
              whiskerprops=dict(color="#9a9994"),capprops=dict(color="#9a9994"),
              boxprops=dict(lw=0))
for patch,v in zip(bp["boxes"],np.linspace(0.20,0.80,len(order))):
    patch.set_facecolor(CMAP(v)); patch.set_edgecolor(SURF); patch.set_linewidth(2)
rng=np.random.default_rng(0)
for i,d in enumerate(data,start=1):
    ax.scatter(rng.normal(i,0.07,len(d)),d,s=3,color=INK,alpha=0.18,linewidths=0,zorder=3)
ax.set_xticks(range(1,len(order)+1))
ax.set_xticklabels([f"{k}\nn={len(sel[k])}" for k in order],fontsize=7.5)
ax.set_ylabel("CLDN4  log$_2$(TPM+1)",color=INK2)
y=max(v.max() for v in data)
for i,((a,b),q) in enumerate(zip(tests,padj)):
    ia,ib=order.index(a)+1,order.index(b)+1
    yy=y+0.6+0.55*i
    ax.plot([ia,ia,ib,ib],[yy-0.12,yy,yy,yy-0.12],lw=0.9,color=INK2)
    s="***" if q<0.001 else "**" if q<0.01 else "*" if q<0.05 else "ns"
    ax.text((ia+ib)/2,yy+0.03,s,ha="center",fontsize=8,color=INK)
ax.set_ylim(top=y+0.6+0.55*len(tests)+0.5)
ax.set_title("CLDN4 in TCGA lung cohorts (UCSC Xena / Toil recompute)",
             fontsize=8.5,color=INK,loc="left")
for s_ in ("top","right"): ax.spines[s_].set_visible(False)
fig.savefig("../figures/Figure1A_TCGA_Xena.pdf",bbox_inches="tight")
fig.savefig("../figures/Figure1A_TCGA_Xena.png",dpi=600,bbox_inches="tight")
print("wrote Figure 1A")
for k in order:
    t=2**sel[k].expr.values-1
    print(f"  {k:12s} n={len(sel[k]):4d}  median log2(TPM+1)={np.median(sel[k].expr):5.2f}  median TPM={np.median(t):7.1f}")
print("  BH-adjusted p:", ", ".join(f"{a} vs {b}: {q:.2e}" for (a,b),q in zip(tests,padj)))
