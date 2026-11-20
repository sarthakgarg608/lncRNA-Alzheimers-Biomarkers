import gzip
from importlib import import_module
import pytest

parse = import_module("scripts.06_build_lncrna_reference").parse_genes

def test_gene_records_and_observed_biotypes(tmp_path):
    path=tmp_path/"genes.gz"
    with gzip.open(path,"wt") as f:
        f.write("##description: test (GRCh38), version 50 (Ensembl 116)\n")
        for i,(kind,typ) in enumerate([("gene","lncRNA"),("exon","lncRNA"),("gene","TEC"),("gene","pseudogene")]):
            f.write(f'chr1\ttest\t{kind}\t1\t9\t.\t+\t.\tgene_id "ENSG{i}.2"; gene_name "TEST"; gene_type "{typ}";\n')
    genes,counts=parse(path)
    assert counts=={"gene":3,"exon":1}
    assert genes.loc[genes.biotype.eq("lncRNA"),"stable_id"].tolist()==["ENSG0"]
    assert genes.strand.tolist()==["+"]*3

