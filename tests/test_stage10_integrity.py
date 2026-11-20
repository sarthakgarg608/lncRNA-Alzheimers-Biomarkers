import json
import pandas as pd
import pytest
from scripts.utils.stage10 import ANN,BLOOD,BRAIN,CFG,OUT,PROC,ROOT,STATES,digest,features,matrix_path

@pytest.mark.parametrize("ds",BLOOD)
def test_real_blood_projection(ds):
    mapping=pd.read_csv(ANN/f"{ds}_probe_mapping.csv",keep_default_na=False)
    assert mapping.probe_id.tolist()==features(ds)
    assert mapping.probe_id.is_unique and mapping.platform_rows.eq(1).all()
    assert set(mapping.mapping_status)<=set(STATES)
    valid=mapping[mapping.mapping_status.eq(STATES[0])]
    assert valid.exact_genes.eq(valid.within_two_substitutions).all()
    assert valid.validated_gene.eq(valid.exact_genes).all() and valid.sense_exact.all()
    assert valid.validated_gene.str.fullmatch(r"ENSG\d+(?:_PAR_Y)?").all()
    ids=valid.probe_id.tolist()
    output=pd.read_parquet(PROC/f"{ds}_validated_lncrna.parquet")
    original=pd.read_parquet(matrix_path(ds),columns=[*ids,"sample_id"])
    pd.testing.assert_frame_equal(output,original)
    stage08=pd.read_parquet(matrix_path(ds).parent/"expression.parquet")
    pd.testing.assert_frame_equal(output,stage08.loc[ids].T,check_names=False)
    assert len(output)==CFG["datasets"][ds]["expected_sample_count"]

@pytest.mark.parametrize("category",["lncrna","protein_coding"])
def test_real_brain_counts(category):
    raw=pd.read_csv(ROOT/CFG["datasets"][BRAIN]["counts_path"],sep="\t",index_col=0).T
    output=pd.read_parquet(PROC/f"{BRAIN}_{category}_raw_counts.parquet")
    raw.columns=raw.columns.astype(str)
    assert output.index.tolist()==raw.index.tolist()
    assert all(pd.api.types.is_integer_dtype(t) for t in output.dtypes)
    assert (raw[output.columns].to_numpy()==output.to_numpy()).all()

def test_raw_digests_against_original_inventory():
    saved=pd.read_csv(ROOT/"outputs/tables/local_inventory.csv")
    for path in (ROOT/"data/raw").glob("*"):
        if path.is_file():
            assert digest(path)==saved.loc[saved.Filename.eq(path.name),"SHA256"].item()

def test_external_lock_and_crosswalk():
    lock=json.loads((matrix_path(BLOOD[2]).parent/"parsing_report.json").read_text())
    assert lock["validation"]["diagnosis_locked"] is True
    ext=pd.read_csv(OUT/"external_feature_compatibility.csv",keep_default_na=False)
    disc=pd.read_csv(ANN/f"{BLOOD[0]}_probe_mapping.csv",keep_default_na=False)
    pairs=set(disc.loc[disc.mapping_status.eq(STATES[0]),["sequence","validated_gene"]].itertuples(index=False,name=None))
    assert ext.identical_validated_measurement.tolist()==[p in pairs for p in zip(ext.sequence,ext.validated_gene)]
    assert not any("diagnos" in c.lower() for c in ext)


def test_brain_namespace_and_outputs_partition():
    brain=pd.read_csv(ANN/f"{BRAIN}_gene_biotypes.csv",dtype=str,keep_default_na=False)
    assert brain.GeneID.tolist()==features(BRAIN)
    assert brain.GeneID.is_unique
    assert brain.category.isin(["lncRNA","protein_coding","unresolved","other"]).all()
    for category in ["lncRNA","protein_coding"]:
        selected=brain.loc[brain.category.eq(category),"GeneID"].tolist()
        columns=pd.read_parquet(PROC/f"{BRAIN}_{category.lower()}_raw_counts.parquet").columns
        assert selected==columns.tolist()
    genes=brain[brain.category.eq("lncRNA")]
    assert genes.Status.eq("active").all() and genes.GeneType.eq("ncRNA").all()

def test_reference_and_annotation_source_integrity():
    reference=pd.read_parquet(ANN/"lncrna_reference.parquet")
    assert reference.biotype.eq("lncRNA").all() and reference.stable_id.is_unique
    sources=json.loads((OUT/"annotation_sources.json").read_text())
    assert len(sources)==6
    for item in sources:
        assert digest(ROOT/item["file"])==item["sha256"]
    for ds in BLOOD:
        m=pd.read_csv(ANN/f"{ds}_probe_mapping.csv",keep_default_na=False)
        validated=m[m.mapping_status.eq(STATES[0])]
        assert validated.validated_gene.isin(reference.stable_id).all()
        assert not validated.conflicting_sequence_or_annotation.any()
