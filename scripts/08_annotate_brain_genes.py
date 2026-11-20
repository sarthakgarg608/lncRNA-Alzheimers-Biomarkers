"""Assign brain gene biotypes using active NCBI records and their explicit crosswalk."""
import re
import numpy as np
import pandas as pd
# Support both direct PowerShell invocation and python -m scripts.<module>.
if __package__ in (None, ""):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.utils.stage10 import (ANN, BRAIN, CFG, OUT, ROOT, annotation, features,
                                 matrix_path, prepare, subset, table)

SHORT = re.compile(r"^(microRNA |small nucleolar RNA|small nuclear RNA|transfer RNA|ribosomal RNA)",re.I)
LONG = re.compile(r"\blong (?:intergenic non-protein coding|non[- ]coding) RNA\b",re.I)

def brain_class(gene_type, description, active, duplicate_link, gencode_lnc, tx_types):
    if not active or duplicate_link:
        return "unresolved","inactive or missing NCBI record" if not active else "shared Ensembl target"
    short = bool(SHORT.search(description)) and "host gene" not in description.lower()
    explicit_long = bool(LONG.search(description))
    coding = "protein_coding" in tx_types
    if gencode_lnc and (gene_type!="ncRNA" or short):
        return "unresolved","NCBI type or description conflicts with GENCODE lncRNA"
    if gene_type=="ncRNA":
        if coding or (explicit_long and tx_types and not gencode_lnc):
            return "unresolved","NCBI ncRNA and cross-release classification require review"
        if gencode_lnc or explicit_long:
            return "lncRNA","NCBI ncRNA + GENCODE lncRNA" if gencode_lnc else "explicit NCBI long RNA description"
        if short:
            return "other","NCBI short RNA description"
        return "unresolved","generic ncRNA without specific long RNA evidence"
    if gene_type=="protein-coding":
        if tx_types and all("pseudogene" in t for t in tx_types):
            return "unresolved","coding versus pseudogene conflict"
        return "protein_coding","active NCBI protein-coding"
    if gene_type in {"","unknown","other"} or (gene_type=="pseudo" and coding):
        return "unresolved","unspecific or conflicting classification"
    return "other","specific NCBI non-lncRNA category"

def main():
    prepare()
    src=ROOT/CFG["datasets"][BRAIN]["counts_path"]
    counts=pd.read_csv(src,sep="\t",index_col="GeneID").T
    if not all(pd.api.types.is_integer_dtype(t) for t in counts.dtypes) or counts.min().min()<0:
        raise ValueError("Raw counts must be nonnegative integers")
    ids=counts.columns.astype(str).tolist()
    if ids!=features(BRAIN) or len(set(ids))!=len(ids) or not all(i.isdecimal() for i in ids):
        raise ValueError("Raw GeneID order or Entrez namespace mismatch")
    saved=pd.read_parquet(matrix_path(BRAIN))
    if not counts.index.equals(saved.index) or not np.array_equal(counts.to_numpy(),saved.to_numpy()):
        raise ValueError("Stage 09 sample order or counts disagree with the original count file")
    ncbi=table(annotation(BRAIN))
    annotated=pd.DataFrame({"GeneID":ids}).merge(ncbi,how="left",validate="one_to_one").fillna("")
    ref=pd.read_parquet(ANN/"lncrna_reference.parquet").set_index("stable_id")
    catalog=pd.read_parquet(ANN/"transcript_gene_catalog.parquet")
    types=catalog.groupby("stable_id").transcript_types.agg("|".join).to_dict()
    dup=set(ncbi.loc[ncbi.EnsemblGeneID.ne("") & ncbi.EnsemblGeneID.duplicated(False),"EnsemblGeneID"])
    records=[]
    for row in annotated.itertuples():
        link=row.EnsemblGeneID
        if link and not re.fullmatch(r"ENSG\d+",link):
            raise ValueError(f"Unexpected NCBI crosswalk namespace {link}")
        lnc=link in ref.index
        tx=types.get(link,"")
        category,why=brain_class(row.GeneType,row.Description,row.Status=="active",link in dup,
                                 lnc,set(tx.split("|")) if tx else set())
        symbol=ref.at[link,"symbol"] if lnc else ""
        notes=[]
        if link and link not in types:
            notes.append("linked ENSG absent from transcript FASTA")
        if symbol and symbol!=row.Symbol:
            notes.append("symbol changed across releases")
        if category=="unresolved":
            notes.append(why)
        records.append((category,why,symbol,tx,"|".join(notes)))
    annotated[["category","evidence","gencode_lnc_symbol","gencode_transcript_types","review_note"]]=records
    annotated["source"]="NCBI Human.GRCh38.p13.annot.tsv.gz + GENCODE v50"
    annotated.to_csv(ANN/f"{BRAIN}_gene_biotypes.csv",index=False)
    annotated.loc[annotated.review_note.ne(""),["GeneID","Symbol","GeneType","Status","EnsemblGeneID",
                  "gencode_lnc_symbol","category","review_note"]].to_csv(OUT/"annotation_discrepancies.csv",index=False)
    annotated.category.value_counts().rename_axis("category").to_csv(OUT/f"{BRAIN}_lncrna_summary.csv",header=["measured_genes"])
    for category in ("lncRNA","protein_coding"):
        subset(BRAIN,annotated.loc[annotated.category.eq(category),"GeneID"].tolist(),category.lower()+"_raw_counts")

if __name__=="__main__":
    main()

