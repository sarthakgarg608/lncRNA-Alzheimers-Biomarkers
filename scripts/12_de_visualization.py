"""Stage 11 descriptive figures from valid, unlocked statistical results."""
import sys
from pathlib import Path
if __package__ in (None,""):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.utils.stage11 import ARRAYS,BRAIN,COHORTS,FIG,OUT,inputs,read_result

INK="#335c70"
ACCENT="#b5462b"
UNITS={ds:("log2 fold change" if ds==BRAIN else "processed intensity difference") for ds in COHORTS}

def save(fig,name):
    fig.savefig(FIG/f"{name}.png",dpi=160,bbox_inches="tight")
    plt.close(fig)

def volcano(ds,r):
    fig,ax=plt.subplots(figsize=(7,4.6),layout="constrained")
    r=r[r.fdr.notna()]
    y=-np.log10(r.fdr.clip(lower=np.finfo(float).tiny))
    for significant,color,label in [(False,INK,"FDR ≥ 0.05"),(True,ACCENT,"FDR < 0.05")]:
        m=r.fdr.lt(.05).eq(significant)
        ax.scatter(r.effect[m],y[m],s=12,alpha=.6,c=color,label=label)
    ax.axhline(-np.log10(.05),ls="--",lw=.8,c="0.5")
    ax.axvline(0,lw=.6,c="0.6")
    best=r.nsmallest(1,"p_value").iloc[0]
    ax.annotate(best.symbol,(best.effect,y[best.name]),xytext=(5,-13),textcoords="offset points",
                ha="right" if best.effect>0 else "left",fontsize=9)
    ax.set(title=ds+" · adjusted AD versus Control",xlabel="AD − Control ("+UNITS[ds]+")",ylabel="−log10(BH FDR)")
    ax.legend(frameon=False)
    save(fig,ds+"_volcano")

def main():
    FIG.mkdir(parents=True,exist_ok=True)
    results={ds:read_result(ds) for ds in COHORTS}
    for ds,r in results.items():
        volcano(ds,r)
    pairs=pd.read_csv(OUT/"blood_replication_summary.csv").dropna(subset=["directional_p"])
    fig,ax=plt.subplots(figsize=(5.5,4.7),layout="constrained")
    for picked,col,label in [(False,INK,"Other matched probes"),(True,ACCENT,"Discovery FDR < 0.05")]:
        sub=pairs[pairs.fdr_d.lt(.05).eq(picked)]
        ax.scatter(sub.standardized_d,sub.standardized_r,s=18,c=col,alpha=.65,label=label)
    ax.axhline(0,c="0.6",lw=.8)
    ax.axvline(0,c="0.6",lw=.8)
    ax.set(title="Identical sequence and gene · blood cohorts",xlabel="Discovery AD coefficient / residual SD",ylabel="Replication AD coefficient / residual SD")
    ax.legend(frameon=False,fontsize=8)
    save(fig,"blood_effect_replication")
    for kind in ["pvalue_distribution","qc_sensitivity_comparison","top_lncrna_expression_boxplots"]:
        fig,axes=plt.subplots(1,3,figsize=(12,4),layout="constrained")
        for ax,(ds,r) in zip(axes,results.items()):
            ax.set_title(ds)
            if kind=="pvalue_distribution":
                ax.hist(r.p_value.dropna(),bins=np.linspace(0,1,21),color=INK)
                ax.set(xlabel="Raw two-sided p-value",ylabel="Valid tests")
            elif kind=="qc_sensitivity_comparison":
                v=r.dropna(subset=["fdr","sensitivity_fdr"])
                ax.scatter(v.effect,v.sensitivity_effect,s=8,alpha=.5,c=INK)
                lo=min(v.effect.min(),v.sensitivity_effect.min())
                hi=max(v.effect.max(),v.sensitivity_effect.max())
                ax.plot([lo,hi],[lo,hi],lw=1,c=ACCENT,ls="--")
                ax.set(xlabel="Primary "+UNITS[ds],ylabel="Sensitivity "+UNITS[ds])
            else:
                best=r.dropna(subset=["fdr"]).sort_values(["fdr","p_value"]).iloc[0]
                expr,meta,_=inputs(ds)
                ledger=pd.read_csv(OUT/f"{ds}_samples.csv",index_col="sample_id")
                vals=expr.loc[ledger.primary,best.feature_id]
                if ds==BRAIN:
                    factors=pd.read_csv(OUT/f"{ds}_primary_size_factors.csv",index_col="sample_id").iloc[:,0]
                    vals=np.log2(1+vals/factors)
                groups=[vals[meta.loc[vals.index,"diagnosis"].eq(d)] for d in ["Control","AD"]]
                ax.boxplot(groups,tick_labels=["Control","AD"],patch_artist=True,boxprops={"facecolor":"#b8d9d8"})
                ax.set_title(f"{ds}\n{best.symbol} · {best.feature_id}",fontsize=10)
                ax.set_ylabel("log2(1 + size-normalized count)" if ds==BRAIN else "Processed intensity")
        fig.suptitle({"pvalue_distribution":"P-values across valid tests","qc_sensitivity_comparison":"AD coefficients before and after severe-QC exclusions",
                      "top_lncrna_expression_boxplots":"Lowest-FDR feature per cohort · unadjusted expression by group"}[kind])
        save(fig,kind)
    cross=pd.read_csv(OUT/"brain_blood_comparison.csv")
    fig,axes=plt.subplots(1,2,figsize=(9,4),layout="constrained")
    for ax,ds in zip(axes,ARRAYS):
        r=cross[cross.dataset_blood.eq(ds)&cross.direction_comparison.ne("not comparable")]
        ax.scatter(r.effect_blood,r.effect_brain,s=13,alpha=.6,c=INK)
        ax.axhline(0,c="0.6",lw=.8)
        ax.axvline(0,c="0.6",lw=.8)
        ax.set(title=ds+" and brain",xlabel="Blood processed intensity difference",ylabel="Brain log2 fold change")
    fig.suptitle("Exploratory cross-tissue associations · one point per linked blood probe")
    save(fig,"brain_blood_comparison")

if __name__=="__main__":
    main()

