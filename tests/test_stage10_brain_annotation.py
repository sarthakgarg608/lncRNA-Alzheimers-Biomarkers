from importlib import import_module
import pytest

classify=import_module("scripts.08_annotate_brain_genes").brain_class

@pytest.mark.parametrize("kind,name,active,duplicate,lnc,expected",[
    ("ncRNA","unknown",True,False,False,"unresolved"),
    ("ncRNA","long intergenic non-protein coding RNA TEST",True,False,False,"lncRNA"),
    ("ncRNA","test",True,False,True,"lncRNA"),
    ("ncRNA","microRNA TEST",True,False,True,"unresolved"),
    ("pseudo","test",True,False,True,"unresolved"),
    ("protein-coding","test",True,False,True,"unresolved"),
    ("protein-coding","test",True,False,False,"protein_coding"),
    ("ncRNA","test",True,True,True,"unresolved"),
    ("ncRNA","test",False,False,True,"unresolved"),
    ("ncRNA","small nucleolar RNA host gene TEST",True,False,True,"lncRNA"),
    ("snoRNA","test",True,False,False,"other"),
])
def test_gene_evidence(kind,name,active,duplicate,lnc,expected):
    assert classify(kind,name,active,duplicate,lnc,set())[0]==expected

