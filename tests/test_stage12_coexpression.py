"""Synthetic algorithm checks and label-independent brain-only output reconciliation."""
from importlib import import_module
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from scipy import stats
from statsmodels.stats.multitest import multipletests

stage=import_module("scripts.13_brain_coexpression")

def test_sample_order_missing_values_and_counts():
    raw=pd.DataFrame([[1,2],[3,4]],index=["GSM1","GSM2"],columns=["10","20"])
    stage.check_axes(raw,raw.copy(),raw.columns)
    with pytest.raises(ValueError,match="sample mismatch"):
        stage.check_axes(raw,raw.iloc[::-1],raw.columns)
    with pytest.raises(ValueError,match="Missing count"):
        stage.check_axes(raw,raw.astype(float).mask(raw.eq(1)),raw.columns)
    with pytest.raises(ValueError,match="original counts"):
        stage.check_axes(raw,raw+1,raw.columns)
    with pytest.raises(ValueError,match="Missing or nonfinite"):
        stage.unit_residuals(np.array([[1.],[np.nan]]),np.ones((2,1)))

def test_global_bh_with_blocks_ties_and_undefined_entries():
    p=np.array([[.04,.02,np.nan],[.02,1.,.0005]])
    result=np.empty_like(p)
    stage.global_bh(p,result,chunk=2)
    valid=np.isfinite(p)
    np.testing.assert_allclose(result[valid],multipletests(p[valid],method="fdr_bh")[1])
    assert np.isnan(result[~valid]).all()

def test_partial_correlation_against_independent_regression():
    rng=np.random.default_rng(42)
    n=40
    dx=np.repeat([0,1],n//2)
    design=np.column_stack([np.ones(n),dx,rng.normal(size=n)])
    a=rng.normal(size=(n,4))+20*dx[:,None]
    b=rng.normal(size=(n,3))+20*dx[:,None]
    b[:,1]+=0.6*a[:,1]
    left,valid=stage.unit_residuals(np.column_stack([a,np.ones(n),dx]),design)
    right,_=stage.unit_residuals(b,design)
    assert valid.tolist()==[True]*4+[False,False]
    blocks=list(stage.correlation_blocks(left,right,n-design.shape[1]-1,block=1))
    p=np.vstack([values for _,_,values in blocks])
    r=np.vstack([values for _,values,_ in blocks])
    exp=[]
    for i in range(4):
        x=np.column_stack([design,a[:,i]])
        df=n-x.shape[1]
        for j in range(3):
            beta=np.linalg.lstsq(x,b[:,j],rcond=None)[0]
            err=b[:,j]-x@beta
            se=np.sqrt((err@err/df)*np.linalg.inv(x.T@x)[-1,-1])
            t=beta[-1]/se
            assert r[i,j]==pytest.approx(t/np.sqrt(df+t*t),abs=1e-10)
            exp.append(2*stats.t.sf(abs(t),df))
    np.testing.assert_allclose(p[:4].ravel(),exp)
    q=np.empty_like(p)
    stage.global_bh(p,q,chunk=3)
    np.testing.assert_allclose(q[:4].ravel(),multipletests(exp,method="fdr_bh")[1])
    assert np.isnan(q[4:]).all()

def test_poscounts_reproduces_installed_stage11_normalization():
    counts=pd.read_parquet(stage.PROC/f"{stage.BRAIN}_lncrna_raw_counts.parquet")
    ledger=pd.read_csv(stage.PRIOR/f"{stage.BRAIN}_samples.csv",index_col="sample_id")
    primary=counts.loc[ledger.primary]
    use=primary.ge(10).sum().ge(10)
    for fit in stage.FITS:
        saved=pd.read_csv(stage.PRIOR/f"{stage.BRAIN}_{fit}_size_factors.csv",index_col="sample_id")
        estimate=stage.poscounts(counts.loc[ledger[fit],use])
        assert estimate.index.equals(saved.index)
        np.testing.assert_allclose(estimate,saved.size_factors,rtol=1e-10)

def test_real_pair_table_and_separate_family_sizes():
    tab=pq.ParquetFile(stage.OUT/"all_pairs.parquet")
    audit=pd.read_csv(stage.OUT/"feature_audit.csv")
    summary=pd.read_csv(stage.OUT/"edges_per_lncrna.csv")
    coding=audit.category.eq("protein_coding")
    assert tab.metadata.num_rows==audit.loc[coding,"eligible"].sum()*audit.loc[~coding,"eligible"].sum()
    first=tab.read_row_group(0).slice(0,1).to_pandas().iloc[0]
    factors=pd.read_csv(stage.OUT/"size_factors.csv",index_col="sample_id")
    meta=pd.read_csv(stage.OUT/"sample_audit.csv",index_col="sample_id")
    for fit in stage.FITS:
        mask=audit[fit+"_variable_residual"].eq(True)
        assert summary[fit+"_valid"].sum()==(mask&coding).sum()*(mask&~coding).sum()
        s=factors[factors.fit.eq(fit)]
        x=stage.design(meta.loc[s.index],stage.BRAIN).to_numpy()
        cols=[]
        for field,label in [("lnc_entrez","lncrna"),("coding_entrez","protein_coding")]:
            gid=str(int(first[field]))
            mat=pd.read_parquet(stage.PROC/f"{stage.BRAIN}_{label}_raw_counts.parquet",columns=[gid])
            v=np.log2(1+mat.loc[s.index,gid].to_numpy()/s.whole.to_numpy())
            cols.append(v-x@np.linalg.lstsq(x,v,rcond=None)[0])
        corr=stats.pearsonr(*cols).statistic
        assert first[fit+"_r"]==pytest.approx(corr)
        df=len(s)-x.shape[1]-1
        assert first[fit+"_p"]==pytest.approx(2*stats.t.sf(abs(corr)*np.sqrt(df/(1-corr**2)),df))
        assert first[fit+"_q"]>=first[fit+"_p"]-1e-14


def test_local_go_background_and_global_correction():
    from scipy.stats import hypergeom
    background=pd.read_csv(stage.OUT/"enrichment_background.csv")
    all_tests=pd.read_parquet(stage.OUT/"GO_BP_all_tests.parquet")
    if all_tests.empty:
        assert pd.read_csv(stage.OUT/"GO_BP_significant.csv.gz").empty
        return
    q=multipletests(all_tests.p,method="fdr_bh")[1]
    np.testing.assert_allclose(all_tests.q,q)
    top=all_tests.loc[all_tests.p.idxmin()]
    foreground=pd.read_csv(stage.OUT/"enrichment_foreground.csv.gz",dtype={"lnc_entrez":str,"coding_entrez":str})
    members=foreground.loc[foreground.lnc_entrez.eq(top.lnc_entrez)&foreground.direction.eq(top.direction),"coding_entrez"]
    ref=pd.read_csv(
        stage.ROOT/stage.CFG["datasets"][stage.BRAIN]["annotation_path"],sep="\t",dtype=str,keep_default_na=False)
    target=set(ref.loc[ref.GOProcessID.str.split("///").map(lambda terms:top.go_id in terms),"GeneID"])
    universe=set(background.coding_entrez.astype(str))
    assert set(members)<=universe
    overlap=len(target.intersection(members))
    size=len(target&universe)
    assert (overlap,size,len(members))==(top.overlap,top.term_size,top.foreground_size)
    assert top.p==pytest.approx(hypergeom.sf(overlap-1,len(universe),size,len(members)),abs=0)

