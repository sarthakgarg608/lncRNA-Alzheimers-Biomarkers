from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm,t
from scripts.utils.stage11 import ARRAYS,BRAIN,COHORTS,OUT,ROOT,bh,inputs,preflight,read_result

@pytest.mark.parametrize("reader",[inputs,preflight,read_result])
def test_external_closed_before_io(reader,monkeypatch):
    def unexpected(*args,**kwargs):
        pytest.fail("locked dataset reached I/O")
    monkeypatch.setattr(pd,"read_csv",unexpected)
    monkeypatch.setattr(pd,"read_parquet",unexpected)
    with pytest.raises(ValueError,match="outside unlocked"):
        reader("GSE140829")

@pytest.mark.parametrize("ds",COHORTS)
def test_full_universe_and_adjusted_statistics(ds):
    result=read_result(ds)
    mat,meta,a=inputs(ds)
    assert result.feature_id.tolist()==mat.columns.tolist()
    valid=result[result.status.eq("valid")]
    assert valid.se.gt(0).all()
    np.testing.assert_allclose(valid.effect/valid.se,valid.statistic,rtol=1e-8)
    probability=2*(norm.sf(abs(valid.statistic)) if ds==BRAIN else t.sf(abs(valid.statistic),valid.t_df))
    np.testing.assert_allclose(probability,valid.p_value,rtol=1e-8,atol=1e-15)
    np.testing.assert_allclose(bh(result.p_value),result.fdr,equal_nan=True,atol=1e-13)
    if ds==BRAIN:
        assert valid[["_LFC_converged","_MAP_converged","_genewise_converged"]].all().all()
    ledger=pd.read_csv(OUT/f"{ds}_samples.csv",index_col="sample_id")
    expected=meta.diagnosis.isin(["AD","Control"]) & meta.age.notna() & meta.sex.notna()
    np.testing.assert_array_equal(ledger.primary,expected)
    assert result.n_samples.eq(expected.sum()).all()
    assert ledger.loc[ledger.primary,"qc_reasons"].notna().any()
    assert ledger.index.equals(mat.index)

def test_candidates_are_associations_only():
    cand=pd.read_csv(OUT/"candidate_lncrna_summary.csv",dtype={"feature_id":str})
    assert cand.fdr.lt(.05).all() and cand.status.eq("valid").all()
    expected={(s,f) for s in COHORTS for f in read_result(s).query("fdr < .05").feature_id}
    assert set(zip(cand.dataset,cand.feature_id))==expected
    assert not any("panel" in c.lower() or "training" in c.lower() for c in cand)


@pytest.mark.parametrize("ds",ARRAYS)
def test_independent_real_ols_and_chip_sandwich(ds):
    from scripts.utils.stage11 import design
    result=read_result(ds).set_index("feature_id")
    expr,meta,_=inputs(ds)
    ledger=pd.read_csv(OUT/f"{ds}_samples.csv",index_col="sample_id")
    meta=meta.loc[ledger.primary]
    x=design(meta,ds)
    y=expr.loc[meta.index,result.index[0]].to_numpy()
    beta=np.linalg.lstsq(x,y,rcond=None)[0]
    residual=y-x.to_numpy()@beta
    chip=meta.title.str.split("_").str[0]
    scores=x.mul(residual,axis=0).groupby(chip).sum().to_numpy()
    bread=np.linalg.inv(x.T@x)
    n,k=x.shape
    g=len(scores)
    variance=bread@(scores.T@scores)@bread*g/(g-1)*(n-1)/(n-k)
    pos=x.columns.get_loc("AD")
    assert result.iloc[0].effect==pytest.approx(beta[pos])
    assert result.iloc[0].se==pytest.approx(np.sqrt(variance[pos,pos]))
