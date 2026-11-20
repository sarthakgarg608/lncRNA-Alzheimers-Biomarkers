"""PyDESeq2 adjusted AD versus Control Wald contrasts for brain raw counts."""
import sys
from pathlib import Path
if __package__ in (None,""):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from threadpoolctl import threadpool_limits
from scripts.utils.stage11 import (BRAIN,OUT,design,finalize,preflight,save_coefficients,validate_counts)

def eligible_genes(counts):
    validate_counts(counts)
    return counts.ge(10).sum().ge(10)

def fit(counts,meta,annotation,eligible,label):
    validate_counts(counts)
    x=design(meta,BRAIN)
    if not counts.index.equals(meta.index):
        raise ValueError("Count/design alignment failure")
    use=eligible.index[eligible]
    dds=DeseqDataSet(counts=counts.loc[:,use],metadata=meta,design=x,
                     size_factors_fit_type="poscounts",refit_cooks=False,n_cpus=1)
    dds.deseq2()
    contrast=x.columns.to_numpy()=="AD"
    stats=DeseqStats(dds,contrast=contrast.astype(float),independent_filter=False,
                     cooks_filter=True,quiet=True,n_cpus=1)
    stats.summary()
    model=stats.results_df.rename(columns={"log2FoldChange":"effect","lfcSE":"se","stat":"statistic","pvalue":"p_value"})
    result=annotation.join(model[["baseMean","effect","se","statistic","p_value"]])
    diagnostics=dds.var[["dispersions","_genewise_converged","_MAP_converged","_LFC_converged","_pvalue_cooks_outlier"]].copy()
    cooks=dds.layers["cooks"]
    diagnostics["max_cooks"]=np.nanmax(cooks,axis=0)
    diagnostics["max_cooks_sample"]=counts.index[np.nanargmax(cooks,axis=0)]
    diagnostics["dispersion_trend"]=dds.uns["disp_function_type"]
    result=result.join(diagnostics)
    converged=diagnostics[["_genewise_converged","_MAP_converged","_LFC_converged"]].fillna(False).all(axis=1)
    result["status"]="low counts"
    result.loc[use,"status"]="valid"
    bad=~np.isfinite(model[["effect","se","statistic","p_value"]]).all(axis=1) | model.se.le(0) | ~np.isfinite(diagnostics.dispersions) | diagnostics.dispersions.le(0)
    result.loc[bad.index[bad],"status"]="invalid estimate"
    result.loc[diagnostics.index[diagnostics._pvalue_cooks_outlier],"status"]="Cook outlier"
    result.loc[converged.index[~converged],"status"]="nonconverged"
    result["n_samples"]=len(meta)
    result["effect_unit"]="log2 fold change"
    result["covariates"]="age_centered + sex"
    result["dataset"]=BRAIN
    factors=dds.obs["size_factors"]
    if not np.isfinite(factors).all() or factors.min()<=0:
        raise ValueError("Invalid size factors")
    factors.to_csv(OUT/f"{BRAIN}_{label}_size_factors.csv",index_label="sample_id")
    result=finalize(result)
    result.drop(columns=annotation.columns).to_csv(OUT/f"{BRAIN}_{label}_diagnostics.csv")
    return result

def main():
    counts,meta,a,ledger=preflight(BRAIN)
    eligible=eligible_genes(counts.loc[ledger.primary])
    with threadpool_limits(limits=1):
        fits=[fit(counts.loc[ledger[key]],meta.loc[ledger[key],["age","sex","diagnosis"]],a,eligible,key)
              for key in ["primary","sensitivity"]]
    save_coefficients(BRAIN,*fits)

if __name__=="__main__":
    main()

