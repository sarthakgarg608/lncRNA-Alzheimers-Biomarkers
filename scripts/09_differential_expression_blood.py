"""OLS baseline for processed array probes, with chip-clustered inference."""
import sys
from pathlib import Path
if __package__ in (None,""):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scripts.utils.stage11 import ARRAYS,design,finalize,preflight,save_coefficients

def fit(ds,expr,meta,annotation):
    if ds not in ARRAYS or not expr.index.equals(meta.index):
        raise ValueError("Expected an aligned, unlocked blood microarray")
    x=design(meta,ds)
    chips=meta.title.str.extract(r"^(\d+)_[A-Za-z0-9]+$",expand=False)
    if chips.isna().any() or chips.nunique()<3:
        raise ValueError("Insufficient chip identity for clustered inference")
    records=[]
    for probe in expr:
        model=sm.OLS(expr[probe].astype(float),x).fit(cov_type="cluster",use_t=True,
                     cov_kwds={"groups":chips.to_numpy(),"use_correction":True,"df_correction":True})
        effect,se,t,p=model.params["AD"],model.bse["AD"],model.tvalues["AD"],model.pvalues["AD"]
        sigma=np.sqrt(model.mse_resid)
        valid=np.isfinite([effect,se,t,p,sigma]).all() and se>0 and sigma>0 and expr[probe].nunique()>1
        records.append((probe,effect,se,t,p,sigma,effect/sigma if valid else np.nan,"valid" if valid else "invalid estimate"))
    result=pd.DataFrame(records,columns=["feature_id","effect","se","statistic","p_value","residual_sd","standardized","status"]).set_index("feature_id")
    result=annotation.join(result)
    result["n_samples"]=len(meta)
    result["t_df"]=chips.nunique()-1
    result["design_rank"]=x.shape[1]
    result["effect_unit"]="processed intensity difference"
    result["covariates"]="age_centered + sex + center; chip-clustered SE"
    result["dataset"]=ds
    return finalize(result)

def main():
    for ds in ARRAYS:
        expr,meta,a,ledger=preflight(ds)
        primary=fit(ds,expr.loc[ledger.primary],meta.loc[ledger.primary],a)
        if ledger.primary.equals(ledger.sensitivity):
            sensitivity=primary.copy()
        else:
            sensitivity=fit(ds,expr.loc[ledger.sensitivity],meta.loc[ledger.sensitivity],a)
        save_coefficients(ds,primary,sensitivity)

if __name__=="__main__":
    main()

