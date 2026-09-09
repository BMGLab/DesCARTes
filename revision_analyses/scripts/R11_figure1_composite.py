#!/usr/bin/env python3
"""
Composite Figure 1 (panels A-D) as a single vector figure.

  A  CLDN4 in TCGA lung cohorts (UCSC Xena, Toil recompute)
  B  Donor-level CLDN4 across tumour stage in the LuCA single-cell atlas
  C  CLDN4 immunohistochemistry, Human Protein Atlas
  D  CLDN4 immunohistochemistry, matched 45-patient cohort
  E  Representative CLDN4 immunohistochemistry, adjacent non-tumour lung
  F  Representative CLDN4 immunohistochemistry, NSCLC tumour

Panels E and F are the original micrographs, extracted at native resolution
(3024 x 4032 RGB) from the embedded images of Gocmenetal_figures_highquality.pdf and
cropped to the circular field. The panel assignment was verified against a 200 dpi
render of the source page rather than assumed from PDF image order, which is inverted
relative to the visual layout: the strongly stained field is the tumour (F) and the
weakly stained field is the adjacent normal lung (E).

NOTE: these panels still require scale bars, objective magnification and specimen
identifiers. That metadata is not recoverable from the PDF and has deliberately not
been invented here.

Outputs: figures/Figure1_composite.{pdf,png}
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcheck import assert_no_text_overlap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy import stats
from statsmodels.stats.multitest import multipletests

O, F = "output/", "../figures/"
BLUE=["#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7","#3987e5",
      "#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b"]
CMAP=LinearSegmentedColormap.from_list("seq_blue",BLUE)
SURF,INK,INK2="#fcfcfb","#0b0b0b","#52514e"
plt.rcParams.update({"font.size":7.5,"axes.edgecolor":"#c9c8c4","axes.labelcolor":INK2,
                     "xtick.color":INK2,"ytick.color":INK2,"figure.facecolor":SURF,
                     "axes.facecolor":SURF,"savefig.facecolor":SURF})
def bare(ax):
    for s in ("top","right"): ax.spines[s].set_visible(False)
def panel(ax,letter):
    ax.text(-0.20,1.10,letter,transform=ax.transAxes,fontsize=11,fontweight="bold",
            color=INK,va="top",ha="left")

fig=plt.figure(figsize=(7.4,10.4))
gs=fig.add_gridspec(4,2,height_ratios=[1.0,1.0,1.30,0.04],hspace=0.62,wspace=0.30)

# ------------------------------------------------------------------ A ------
axA=fig.add_subplot(gs[0,0])
m=pd.read_csv(O+"xena_cldn4.csv")
m["expr"]=np.log2(np.maximum(2**m.CLDN4.astype(float)-0.001,0)+1)
sel={"Normal\nlung": m[(m._study=="TCGA")&(m._sample_type=="Solid Tissue Normal")
                       &(m.detailed_category.astype(str).str.contains("Lung",na=False))],
     "LUAD":        m[(m._study=="TCGA")&(m.detailed_category=="Lung Adenocarcinoma")
                       &(m._sample_type=="Primary Tumor")],
     "LUSC":        m[(m._study=="TCGA")&(m.detailed_category=="Lung Squamous Cell Carcinoma")
                       &(m._sample_type=="Primary Tumor")]}
oA=list(sel); dA=[sel[k].expr.values for k in oA]
tests=[("Normal\nlung","LUAD"),("LUAD","LUSC"),("Normal\nlung","LUSC")]
pv=[stats.mannwhitneyu(sel[a].expr,sel[b].expr,alternative="two-sided").pvalue for a,b in tests]
padj=multipletests(pv,method="fdr_bh")[1]
bp=axA.boxplot(dA,patch_artist=True,widths=0.5,showfliers=False,
               medianprops=dict(color=INK,lw=1.2),whiskerprops=dict(color="#9a9994"),
               capprops=dict(color="#9a9994"),boxprops=dict(lw=0))
for pt,v in zip(bp["boxes"],np.linspace(0.20,0.80,len(oA))):
    pt.set_facecolor(CMAP(v)); pt.set_edgecolor(SURF); pt.set_linewidth(1.5)
rng=np.random.default_rng(0)
for i,d in enumerate(dA,start=1):
    axA.scatter(rng.normal(i,0.06,len(d)),d,s=1.6,color=INK,alpha=0.14,linewidths=0,zorder=3)
y=max(v.max() for v in dA)
for i,((a,b),q) in enumerate(zip(tests,padj)):
    ia,ib=oA.index(a)+1,oA.index(b)+1; yy=y+0.55+0.62*i
    axA.plot([ia,ia,ib,ib],[yy-0.14,yy,yy,yy-0.14],lw=0.8,color=INK2)
    s="***" if q<0.001 else "**" if q<0.01 else "*" if q<0.05 else "ns"
    axA.text((ia+ib)/2,yy+0.02,s,ha="center",fontsize=7.5,color=INK)
axA.set_ylim(top=y+0.55+0.62*len(tests)+0.55)
axA.set_xticks(range(1,len(oA)+1))
axA.set_xticklabels([f"{k}\nn={len(sel[k])}" for k in oA],fontsize=7)
axA.set_ylabel("CLDN4  log$_2$(TPM+1)",color=INK2)
axA.set_title("TCGA lung cohorts",fontsize=8,color=INK,loc="left")
bare(axA); panel(axA,"A")

# ------------------------------------------------------------------ B ------
axB=fig.add_subplot(gs[0,1])
prof=pd.read_csv(O+"luca_cldn4_profiles.csv")
don=prof.groupby(["donor_id","group"],observed=True)[["g","tot"]].sum().reset_index()
don["expr"]=np.log1p(don.g/don.tot*1e6)
ORD=["Epithelial (normal lung)","I","II","III","III or IV","IV"]
LAB=["Normal\nepith.","I","II","III","III/IV","IV"]
st=pd.read_csv(O+"luca_stage_stats_FINAL.csv"); st=st[st.unit=="donor-level"].set_index("stage")
dB=[don.loc[don.group==g,"expr"].values for g in ORD]
bp=axB.boxplot(dB,patch_artist=True,widths=0.52,showfliers=False,
               medianprops=dict(color=INK,lw=1.2),whiskerprops=dict(color="#9a9994"),
               capprops=dict(color="#9a9994"),boxprops=dict(lw=0))
for pt,c in zip(bp["boxes"],[CMAP(v) for v in np.linspace(0.18,0.88,len(ORD))]):
    pt.set_facecolor(c); pt.set_edgecolor(SURF); pt.set_linewidth(1.5)
for i,d in enumerate(dB,start=1):
    axB.scatter(rng.normal(i,0.065,len(d)),d,s=1.8,color=INK,alpha=0.22,linewidths=0,zorder=3)
yB=max(v.max() for v in dB)
for i,g in enumerate(ORD[1:],start=2):
    if g in st.index:
        sig=st.loc[g,"signif_BH"]
        axB.text(i,yB+0.30,sig,ha="center",va="bottom",fontsize=7.5,
                 color=INK if sig!="ns" else INK2)
axB.set_ylim(top=yB+1.05)
axB.set_xticks(range(1,len(ORD)+1))
axB.set_xticklabels([f"{l}\nn={len(d)}" for l,d in zip(LAB,dB)],fontsize=6.6)
axB.set_ylabel("CLDN4  log(CPM+1)",color=INK2)
axB.set_title("LuCA atlas, donor level",fontsize=8,color=INK,loc="left")
bare(axB); panel(axB,"B")

# --------------------------------------------------------------- C and D ---
LV=["Not detected","Low","Medium","High"]; cols=[CMAP(v) for v in (0.10,0.35,0.62,0.88)]
ihc=pd.read_csv(O+"TableS_IHC_pairs.csv"); t=ihc.tumour_intensity.value_counts()
sets=[("C",gs[1,0],{"Normal lung\nn=4":[3,1,0,0],"NSCLC\nn=11":[0,2,9,0]},
       "Human Protein Atlas",""),
      ("D",gs[1,1],{"Adjacent\nnormal\nn=45":[45,0,0,0],
                    "NSCLC\ntumour\nn=45":[0,int(t.get(1,0)),int(t.get(2,0)),int(t.get(3,0))]},
       "Matched cohort","all 45 pairs concordant; exact Wilcoxon\nsigned-rank p = 5.7 × 10$^{-14}$")]
for letter,slot,dat,title,note in sets:
    ax=fig.add_subplot(slot); keys=list(dat); bot=np.zeros(len(keys))
    for li,lv in enumerate(LV):
        vals=np.array([dat[k][li] for k in keys],float)
        pct=vals/np.array([sum(dat[k]) for k in keys])*100
        ax.bar(keys,pct,bottom=bot,color=cols[li],edgecolor=SURF,linewidth=1.5,width=0.55)
        for xi,(p_,b_) in enumerate(zip(pct,bot)):
            if p_>=9:
                ax.text(xi,b_+p_/2,f"{p_:.0f}%",ha="center",va="center",fontsize=6.8,
                        color="#ffffff" if li>=2 else INK)
        bot+=pct
    ax.set_ylabel("% of cases",color=INK2); ax.set_ylim(0,100)
    ax.tick_params(axis="x",labelsize=6.8)
    # the note goes in the free strip between the axes and the title, so it can never
    # collide with whatever is laid out below this panel
    if note:
        ax.set_title(title,fontsize=8,color=INK,loc="left",pad=22)
        ax.text(0,1.015,note,transform=ax.transAxes,fontsize=6,color=INK2,va="bottom")
    else:
        ax.set_title(title,fontsize=8,color=INK,loc="left")
    bare(ax); panel(ax,letter)

# shared legend for C and D, in its own row so it cannot overlap either panel
axL=fig.add_subplot(gs[3,:]); axL.axis("off")
axL.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(cols,LV)],
           ncol=4,fontsize=7,frameon=False,loc="center",title="CLDN4 staining",
           title_fontsize=7)

# ------------------------------------------------------------- E and F ----
import matplotlib.image as mpimg
MIC=[("E",gs[2,0],"source_micrographs/E_CLDN4_normal.png",
      "Adjacent non-tumour lung"),
     ("F",gs[2,1],"source_micrographs/F_CLDN4_tumour.png",
      "NSCLC tumour")]
for letter,slot,path,cap in MIC:
    ax=fig.add_subplot(slot)
    ax.imshow(mpimg.imread(F+path))
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(cap,fontsize=8,color=INK,loc="left")
    panel(ax,letter)

assert_no_text_overlap(fig, "Figure1_composite")

fig.savefig(F+"Figure1_composite.pdf",bbox_inches="tight")
fig.savefig(F+"Figure1_composite.png",dpi=600,bbox_inches="tight")
print("wrote Figure1_composite")
print("  A  TCGA (Xena):", ", ".join(f"{k.replace(chr(10),' ')} n={len(sel[k])} median={np.median(sel[k].expr):.2f}" for k in oA))
print("     BH p:", ", ".join(f"{a.replace(chr(10),' ')} vs {b}: {q:.2e}" for (a,b),q in zip(tests,padj)))
