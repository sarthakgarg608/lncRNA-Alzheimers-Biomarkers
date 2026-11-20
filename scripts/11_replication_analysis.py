"""Probe-level replication, QC contrasts and exploratory cross-tissue summaries."""
import sys
from pathlib import Path
if __package__ in (None,""):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from scripts.utils.stage11 import ARRAYS,BRAIN,COHORTS,OUT,bh,read_result

SUPPORT="Replicated with supporting evidence"
DIRECTIONAL="Directionally consistent but not statistically replicated"
INCONSISTENT="Inconsistent"
INSUFFICIENT="Insufficient technical evidence"

def probe_disagreement(result):
    signs=result[result.status.eq("valid")].groupby("gene_id").effect.agg(lambda s:s.min()<0<s.max())
    return result.gene_id.map(signs)

def replication(discovery,replication):
    fields=["feature_id","gene_id","sequence","symbol","effect","standardized","p_value","fdr","sensitivity_effect","sensitivity_status"]
    left=discovery[fields].assign(opposed_probes=probe_disagreement(discovery))
    right=replication[fields].assign(opposed_probes=probe_disagreement(replication))
    pairs=left.merge(right,on=["gene_id","sequence"],how="left",suffixes=("_d","_r"))
    valid=pairs.p_value_d.notna() & pairs.p_value_r.notna() & pairs.effect_d.ne(0)
    aligned=(pairs.effect_d*pairs.effect_r)>0
    directional=np.where(aligned,pairs.p_value_r/2,1-pairs.p_value_r/2)
    pairs["directional_p"]=pd.Series(directional,index=pairs.index).where(valid)
    pairs["directional_fdr"]=bh(pairs.directional_p)
    stable=(pairs.effect_d*pairs.sensitivity_effect_d>0)&(pairs.effect_r*pairs.sensitivity_effect_r>0)
    stable &= pairs.sensitivity_status_d.eq("valid") & pairs.sensitivity_status_r.eq("valid")
    support=aligned & pairs.fdr_d.lt(.05) & pairs.directional_fdr.lt(.05) & stable
    pairs["replication"]=np.select([~valid,support,aligned],[INSUFFICIENT,SUPPORT,DIRECTIONAL],default=INCONSISTENT)
    return pairs

def cross_tissue(blood,brain):
    keys=["feature_id","gene_id","symbol","dataset","effect","p_value","fdr","status","direction"]
    linked=brain[brain.gene_id.notna() & brain.gene_id.ne("")]
    cross=blood[keys].merge(linked[keys],on="gene_id",how="left",validate="many_to_one",suffixes=("_blood","_brain"))
    valid=cross.fdr_brain.notna() & cross.fdr_blood.notna()
    cross["direction_comparison"]=np.select([~valid,cross.effect_blood*cross.effect_brain>0,
                                             cross.effect_blood*cross.effect_brain<0],
                                            ["not comparable","consistent","opposite"],default="zero")
    cross["significance"]=np.select([~valid,cross.fdr_blood.lt(.05)&cross.fdr_brain.lt(.05),
                                      cross.fdr_blood.lt(.05),cross.fdr_brain.lt(.05)],
                                     ["not comparable","both tissues","blood only","brain only"],default="neither")
    cross["scope"]="exploratory cross-tissue biological context"
    return cross

def main():
    results={ds:read_result(ds) for ds in COHORTS}
    pairs=replication(*(results[d] for d in ARRAYS))
    pairs.to_csv(OUT/"blood_replication_summary.csv",index=False)
    blood=pd.concat([results[d] for d in ARRAYS],ignore_index=True)
    cross=cross_tissue(blood,results[BRAIN])
    cross.to_csv(OUT/"brain_blood_comparison.csv",index=False)
    candidates=[]
    audits=[]
    for ds,r in results.items():
        both=r.fdr.notna() & r.sensitivity_fdr.notna()
        primary=r.fdr.lt(.05)
        sensitive=r.sensitivity_fdr.lt(.05)
        ledger=pd.read_csv(OUT/f"{ds}_samples.csv")
        audits.append(dict(dataset=ds,primary_n=int(ledger.primary.sum()),sensitivity_n=int(ledger.sensitivity.sum()),
                           eligible_flagged=int((ledger.primary & ledger.qc_reasons.notna()).sum()),
                           valid_primary=r.fdr.notna().sum(),valid_sensitivity=r.sensitivity_fdr.notna().sum(),
                           paired_valid=both.sum(),same_sign=((r.effect*r.sensitivity_effect>0)&both).sum(),
                           primary_fdr05=primary.sum(),sensitivity_fdr05=sensitive.sum(),
                           retained=(primary & sensitive).sum(),lost=(primary & ~sensitive).sum(),
                           gained=(~primary & sensitive).sum()))
        cand=r[primary].copy()
        cand["probe_id"]=cand.feature_id if ds in ARRAYS else ""
        cand["stable_gene_id"]=cand.gene_id if ds in ARRAYS else "NCBIGene:"+cand.feature_id
        cand["qc_sensitivity"]=np.where(cand.sensitivity_fdr.lt(.05)&(cand.effect*cand.sensitivity_effect>0),
                                       "same sign and FDR < 0.05","primary finding not retained")
        if ds in ARRAYS:
            key="feature_id_d" if ds==ARRAYS[0] else "feature_id_r"
            rep=pairs.groupby(key).replication.agg(lambda x:" | ".join(sorted(set(x))))
            cand["replication_evidence"]=cand.feature_id.map(rep).fillna(INSUFFICIENT)
            context=cross[cross.dataset_blood.eq(ds)].set_index("feature_id_blood")
            cand["brain_evidence"]=cand.feature_id.map(context.direction_comparison+" / "+context.significance)
            cand["caveat"]="probe-level clustered OLS; related blood cohorts"
        else:
            rep=pairs.groupby("gene_id").replication.agg(lambda x:" | ".join(sorted(set(x))))
            cand["replication_evidence"]=cand.gene_id.map(rep).fillna("no discovery match")
            context=cross[cross.feature_id_brain.notna()].groupby("feature_id_brain").direction_comparison.agg(lambda x:" | ".join(sorted(set(x))))
            cand["brain_evidence"]=cand.feature_id.map(context).fillna("no linked blood probe")
            cand["caveat"]="postmortem tissue association; bulk cell composition unresolved"
        candidates.append(cand)
    pd.DataFrame(audits).to_csv(OUT/"qc_sensitivity_summary.csv",index=False)
    pd.concat(candidates,ignore_index=True).sort_values(["dataset","fdr","feature_id"]).to_csv(OUT/"candidate_lncrna_summary.csv",index=False)

if __name__=="__main__":
    main()

