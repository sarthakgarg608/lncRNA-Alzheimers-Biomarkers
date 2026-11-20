from importlib import import_module
import pandas as pd

de=import_module("scripts.11_replication_analysis")

def test_requires_sequence_identity_and_preserves_opposed_probes():
    a=pd.DataFrame({"feature_id":["A","B","C"],"sequence":["A"*50,"C"*50,"G"*50],"gene_id":["G"]*3,
                    "symbol":["G"]*3,"status":["valid"]*3,"effect":[1,-1,.2],"standardized":[1,-1,.2],
                    "p_value":[.001]*3,"fdr":[.002]*3,"sensitivity_effect":[.8,-.9,.1],
                    "sensitivity_status":["valid"]*3})
    b=a.copy()
    b["feature_id"]=["D","E","F"]
    b.loc[1,"effect"]=1
    b.loc[2,"sequence"]="T"*50
    out=de.replication(a,b)
    assert out.replication.tolist()==[de.SUPPORT,de.INCONSISTENT,de.INSUFFICIENT]
    assert out.opposed_probes_d.all()
    b.loc[0,"sensitivity_effect"]=-1
    assert de.replication(a,b).replication.iloc[0]==de.DIRECTIONAL

