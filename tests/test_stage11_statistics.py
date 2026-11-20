import numpy as np
import pandas as pd
import pytest
from scripts.utils.stage11 import bh,design,ARRAYS

def test_bh_preserves_missing_tests_and_order():
    q=bh(pd.Series([.04,.001,np.nan,.03]))
    np.testing.assert_allclose(q,[.04,.003,np.nan,.04],equal_nan=True)
    with pytest.raises(ValueError):
        bh([1.1])

def test_design_rejects_confounded_sex():
    m=pd.DataFrame({"diagnosis":["AD","Control"]*5,"age":range(10),"sex":["M","F"]*5,"batch":"A"})
    with pytest.raises(ValueError,match="Rank deficient"):
        design(m,ARRAYS[0])

def test_design_rejects_third_label():
    m=pd.DataFrame({"diagnosis":["AD","MCI","Control"],"age":[70,71,72]})
    with pytest.raises(ValueError):
        design(m,ARRAYS[0])

