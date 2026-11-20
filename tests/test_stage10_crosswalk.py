import pytest
from scripts.utils.stage10 import decode_nuid, stable

def test_checksum_and_padding():
    # Synthetic ACG = 6, checksum 6, no padding.
    assert decode_nuid("SG")=="ACG"
    assert decode_nuid("TG")==""
    assert decode_nuid("AG")==""
    assert decode_nuid("ILMN_123")==""

def test_namespaces():
    assert stable("ENSG0001.4_PAR_Y")=="ENSG0001_PAR_Y"
    assert stable("ENSG0001.4")=="ENSG0001"
    assert stable("123")!="ENSG123"

