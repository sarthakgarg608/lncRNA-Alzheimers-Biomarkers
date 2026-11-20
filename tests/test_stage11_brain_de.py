from importlib import import_module
import numpy as np
import pandas as pd
import pytest
from scripts.utils.stage11 import validate_counts

eligible=import_module("scripts.10_differential_expression_brain").eligible_genes

def test_fixed_low_count_boundary():
    counts=pd.DataFrame({"yes":[10]*10,"few":[10]*9+[9],"weak":[9]*10})
    assert eligible(counts).to_dict()=={"yes":True,"few":False,"weak":False}

@pytest.mark.parametrize("values",[[6.4,7.0],[-1,0],[1.0,np.nan]])
def test_reject_processed_or_invalid_counts(values):
    with pytest.raises(ValueError,match="integer counts"):
        validate_counts(pd.DataFrame({"g":values}))

