import gzip
from importlib import import_module
import pytest
from scripts.utils.stage10 import align

classify=import_module("scripts.07_annotate_blood_probes").classify

@pytest.mark.parametrize("exact,near,historical,sense,conflict,status",[
    (set(),set(),{"L"},False,False,"PROVISIONAL"),
    (set(),set(),set(),False,False,"UNMAPPED"),
    ({"L"},{"L"},{"L"},True,False,"VALIDATED_UNIQUE_LNCRNA"),
    ({"L"},{"L","P"},set(),True,False,"AMBIGUOUS"),
    ({"L"},{"L"},{"P"},True,False,"AMBIGUOUS"),
    ({"L"},{"L"},set(),True,True,"AMBIGUOUS"),
    ({"L"},{"L"},set(),False,False,"PROVISIONAL"),
    ({"P"},{"P"},set(),True,False,"NON_LNCRNA")])
def test_evidence(exact,near,historical,sense,conflict,status):
    assert classify(exact,near,historical,{"L"},sense,{"P"},conflict)==status

def test_seed_scanner_both_strands_replication_and_mismatches(tmp_path):
    import random
    rng=random.Random(42)
    seq="".join(rng.choices("ACGT",k=50))
    fasta=tmp_path/"transcripts.fa.gz"
    with gzip.open(fasta,"wt") as f:
        for gene,n in enumerate((0,0,1,2,3)):
            bases=list(seq)
            for pos in [2,20,40][:n]:
                bases[pos]=next(b for b in "ACGT" if b!=bases[pos])
            target="".join(bases)
            if n==2:
                target=target.translate(str.maketrans("ACGT","TGCA"))[::-1]
            f.write(f">ENST{gene}|ENSG{max(1,gene)}|-|-|T|G|50|lncRNA|\n{target[:23]}\n{target[23:]}\n")
    out=tmp_path/"hits"
    align([seq],fasta,out)
    rows=[line.split("\t") for line in out.read_text().splitlines()]
    assert [(r[1],r[2],r[4]) for r in rows]==[("ENSG1","0","+"),("ENSG2","1","+"),("ENSG3","2","-")]


def test_duplicate_platform_rows_remain_visible(monkeypatch):
    import pandas as pd
    from scripts.utils import stage10
    platform=pd.DataFrame({"ID":["ILMN_01","ILMN_01"],"SEQUENCE":["A"*50,"C"*50]})
    monkeypatch.setattr(stage10,"table",lambda path:platform)
    monkeypatch.setattr(stage10,"features",lambda ds:["ILMN_01","ILMN_02"])
    data=stage10.measured_platform(stage10.BLOOD[0])
    assert data.ID.tolist()==["ILMN_01","ILMN_01","ILMN_02"]
    assert data.sequence.tolist()==["A"*50,"C"*50,""]
    assert classify({"L"},{"L"},set(),{"L"},set(),True,data.sequence.iloc[:2].nunique()>1)=="AMBIGUOUS"

