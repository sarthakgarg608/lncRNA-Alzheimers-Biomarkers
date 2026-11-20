from importlib import import_module
import numpy as np
import pandas as pd
import pytest
from scripts.utils.stage11 import ARRAYS,design

fit=import_module("scripts.09_differential_expression_blood").fit

def test_adjusted_coefficient_and_units():
    rng=np.random.default_rng(42)
    n=80
    meta=pd.DataFrame({"diagnosis":["AD","Control"]*(n//2),"age":rng.integers(60,91,n),
                       "sex":rng.choice(["Female","Male"],n),"batch":rng.choice(["C1","C2"],n),
                       "title":[f"{100+i//4}_{i%4}" for i in range(n)]})
    x=design(meta,ARRAYS[0])
    np.testing.assert_array_equal(x.AD,meta.diagnosis.eq("AD"))
    noise=rng.normal(size=n)
    residual=noise-x@np.linalg.lstsq(x,noise,rcond=None)[0]
    expr=pd.DataFrame({"p1":5+2*x.AD+.05*x.age_centered+residual})
    a=pd.DataFrame({"gene_id":["SYNTHETIC"]},index=pd.Index(["p1"],name="feature_id"))
    out=fit(ARRAYS[0],expr,meta,a).iloc[0]
    assert out.effect==pytest.approx(2)
    assert out.t_df==19 and out.direction=="higher in AD"
    assert out.effect_unit=="processed intensity difference"
    with pytest.raises(ValueError,match="aligned"):
        fit(ARRAYS[0],expr.iloc[::-1],meta,a)

