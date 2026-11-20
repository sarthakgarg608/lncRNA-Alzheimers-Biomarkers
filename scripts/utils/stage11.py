"""Stage 11 input gates, sample eligibility and shared inference conventions."""
import csv
import gzip
import logging
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from .stage10 import ANN, BRAIN, CFG, PROC, ROOT, STATES, digest, matrix_path

ARRAYS=("GSE63061","GSE63060")
COHORTS=(*ARRAYS,BRAIN)
OUT=ROOT/"outputs/tables/stage11"
FIG=ROOT/"outputs/figures/stage11"
OUT.mkdir(parents=True,exist_ok=True)
logging.basicConfig(format="%(asctime)s %(message)s",level=logging.INFO)
log=logging.getLogger("stage11")

def unlocked(ds):
    if ds not in COHORTS or CFG["datasets"][ds]["locked"]:
        raise ValueError(f"Dataset outside unlocked Stage 11 protocol: {ds}")

def bh(p):
    p=pd.Series(p,copy=True,dtype=float)
    valid=p.notna()
    if not p[valid].between(0,1).all():
        raise ValueError("Invalid probability")
    q=pd.Series(np.nan,index=p.index)
    if valid.any():
        q.loc[valid]=multipletests(p[valid],method="fdr_bh")[1]
    return q

def design(meta,ds):
    unlocked(ds)
    if set(meta.diagnosis)!={"AD","Control"} or meta.age.isna().any():
        raise ValueError("Contrast requires complete age and both eligible diagnoses")
    cov=["sex"] if ds==BRAIN else ["sex","batch"]
    if meta[cov].isna().any().any() or meta.sex.nunique()!=2:
        raise ValueError("Incomplete or nonestimable covariates")
    x=pd.DataFrame({"intercept":1.0,"AD":meta.diagnosis.eq("AD").astype(float),
                    "age_centered":meta.age-meta.age.mean()},index=meta.index)
    x=x.join(pd.get_dummies(meta[cov],drop_first=True,dtype=float))
    if np.linalg.matrix_rank(x)!=x.shape[1] or len(x)<=x.shape[1]:
        raise ValueError(f"Rank deficient design for {ds}")
    return x

def validate_counts(frame):
    if not all(pd.api.types.is_integer_dtype(t) for t in frame.dtypes) or (frame.to_numpy()<0).any():
        raise ValueError("Brain requires nonnegative raw integer counts")

def inputs(ds):
    unlocked(ds)
    meta=pd.read_csv(matrix_path(ds).parent/"metadata.csv").set_index("sample_id")
    if ds==BRAIN:
        path=PROC/f"{ds}_lncrna_raw_counts.parquet"
        mapping=ANN/f"{ds}_gene_biotypes.csv"
        a=pd.read_csv(mapping,dtype=str,keep_default_na=False).query("category=='lncRNA'").copy()
        if not (a.GeneType.eq("ncRNA") & a.Status.eq("active")).all():
            raise ValueError("Ineligible brain annotation")
        a=a.rename(columns={"GeneID":"feature_id","EnsemblGeneID":"gene_id","Symbol":"symbol","evidence":"confidence"})
    else:
        path=PROC/f"{ds}_validated_lncrna.parquet"
        mapping=ANN/f"{ds}_probe_mapping.csv"
        a=pd.read_csv(mapping,keep_default_na=False)
        a=a[a.mapping_status.eq(STATES[0])].copy()
        if not (a.exact_genes.eq(a.validated_gene)&a.within_two_substitutions.eq(a.validated_gene)
                & a.sense_exact & ~a.conflicting_sequence_or_annotation & a.sequence.str.fullmatch("[ACGT]{50}")).all():
            raise ValueError("Invalid blood sequence evidence")
        ref=pd.read_parquet(ANN/"lncrna_reference.parquet").set_index("stable_id")
        a["symbol"]=ref.loc[a.validated_gene,"symbol"].to_numpy()
        a=a.rename(columns={"probe_id":"feature_id","validated_gene":"gene_id","mapping_status":"confidence"})
    a=a.set_index("feature_id")
    if not meta.index.is_unique or not a.index.is_unique:
        raise ValueError(f"Duplicate sample or feature identifier in {ds}")
    mat=pd.read_parquet(path)
    if not mat.index.equals(meta.index) or not mat.columns.equals(a.index) or not np.isfinite(mat.to_numpy()).all():
        raise ValueError(f"Invalid sample/feature alignment or nonfinite data in {ds}")
    original=pd.read_parquet(matrix_path(ds),columns=[*a.index,"sample_id"])
    original.columns=original.columns.astype(str)
    pd.testing.assert_frame_equal(mat,original)
    if ds==BRAIN:
        validate_counts(mat)
    fields=["gene_id","symbol","confidence"]+([] if ds==BRAIN else ["sequence"])
    log.info("%s verified %s, matrix SHA256 %s, annotation SHA256 %s",ds,mat.shape,digest(path),digest(mapping))
    return mat,meta,a[fields]

def eligibility(ds,meta):
    unlocked(ds)
    ledger=meta[["diagnosis","age","sex","batch"]].copy()
    binary=ledger.diagnosis.isin(["AD","Control"])
    cov=["age","sex"]+([] if ds==BRAIN else ["batch"])
    ledger["primary"]=binary & ledger[cov].notna().all(axis=1)
    qcdir=ROOT/"outputs/tables/stage09"
    flags=pd.read_csv(qcdir/"sample_qc_flags.csv").query("dataset_id==@ds")
    reasons=flags.groupby("sample_id").flag_reason.agg("|".join)
    ledger["qc_reasons"]=ledger.index.map(reasons).fillna("")
    metrics=pd.read_csv(qcdir/"sample_qc_summary.csv",usecols=["dataset_id","sample_id","library_size","pct_zero_genes"])
    metrics=metrics[metrics.dataset_id.eq(ds)].set_index("sample_id").loc[ledger.index]
    if ds==BRAIN:
        severe=(ledger.qc_reasons.str.contains("Extreme library size") & ledger.qc_reasons.str.contains("High zero fraction")
                & metrics.library_size.lt(metrics.library_size.median()) & metrics.pct_zero_genes.gt(metrics.pct_zero_genes.median()))
        ledger=ledger.join(metrics[["library_size","pct_zero_genes"]])
    else:
        severe=ledger.qc_reasons.str.contains("Extreme IQR") & ledger.qc_reasons.str.contains("Extreme median")
    ledger["sensitivity"]=ledger.primary & ~severe
    ledger["exclusion_reason"]=np.select([~binary,~ledger.primary,severe],
                                          ["outside AD/Control","missing covariate","sensitivity only: concurrent technical extremes"],default="")
    return ledger

def preflight(ds):
    mat,meta,a=inputs(ds)
    ledger=eligibility(ds,meta)
    with gzip.open(ROOT/CFG["datasets"][ds]["series_matrix_path"],"rt") as stream:
        for line in stream:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                ids=next(csv.reader([line],delimiter="\t"))[1:]
            if line.startswith('!Sample_characteristics_ch1\t"age:') and ds==BRAIN:
                ages=[v.split(":",1)[1].strip() for v in next(csv.reader([line],delimiter="\t"))[1:]]
                pd.testing.assert_series_equal(pd.to_numeric(pd.Series(ages,index=ids),errors="coerce"),meta.age,check_names=False)
            if ((ds==BRAIN and line.startswith('!Sample_characteristics_ch1\t"case id:')) or
                (ds in ARRAYS and line.startswith("!Sample_description"))):
                subjects=next(csv.reader([line],delimiter="\t"))[1:]
    if ids!=meta.index.tolist() or len(set(subjects))!=len(ids):
        raise ValueError("Raw GSM/subject alignment failure")
    ledger["reported_sample_code"]=subjects
    for fit in ["primary","sensitivity"]:
        x=design(meta.loc[ledger[fit]],ds)
        log.info("%s %s full rank %d, n=%d, condition number %.2f",ds,fit,x.shape[1],len(x),np.linalg.cond(x))
    ledger.to_csv(OUT/f"{ds}_samples.csv")
    return mat,meta,a,ledger

def finalize(result):
    valid=result.status.eq("valid")
    result.loc[~valid,"p_value"]=np.nan
    result["fdr"]=bh(result.p_value)
    result["direction"]=np.select([valid & result.effect.gt(0),valid & result.effect.lt(0),valid],
                                   ["higher in AD","lower in AD","zero"],default="unavailable")
    return result

def save_coefficients(ds,primary,sensitivity):
    unlocked(ds)
    combined=primary.join(sensitivity[["effect","p_value","fdr","status","n_samples"]].add_prefix("sensitivity_"))
    combined.to_csv(OUT/f"{ds}_differential_expression.csv")
    return combined

def read_result(ds):
    unlocked(ds)
    return pd.read_csv(OUT/f"{ds}_differential_expression.csv",dtype={"feature_id":str,"gene_id":str})

