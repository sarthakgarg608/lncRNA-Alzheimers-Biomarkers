"""Build the gene-level GENCODE reference and audit the six local annotations."""
import gzip
import json
import re
from collections import Counter, defaultdict
import pandas as pd
# Support both direct PowerShell invocation and python -m scripts.<module>.
if __package__ in (None, ""):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.utils.stage10 import ANN, BLOOD, BRAIN, OUT, ROOT, annotation, decode_nuid, digest, prepare, source, stable, table

def parse_genes(path):
    rows = []
    records = Counter()
    with gzip.open(path, "rt") as stream:
        header = stream.readline()
        version = re.search(r"\((GRCh\d+)\), version (\d+) \(Ensembl (\d+)\)", header)
        if not version:
            raise ValueError("Unknown GENCODE assembly or release")
        for line in stream:
            if line.startswith("#"):
                continue
            chrom, provider, feature, start, end, _, strand, _, attrs = line.rstrip().split("\t")
            records[feature] += 1
            if feature != "gene":
                continue
            a = dict(re.findall(r'(\w+) "([^"]+)";', attrs))
            gid = a["gene_id"]
            if not re.fullmatch(r"ENSG\d+\.\d+(?:_PAR_Y)?", gid) or int(start)>int(end):
                raise ValueError(f"Unexpected GTF gene record {gid}")
            rows.append((gid,stable(gid),a["gene_name"],a["gene_type"],chrom,int(start),int(end),
                         strand,version[2],version[1],version[3]))
    genes = pd.DataFrame(rows, columns=["gene_id","stable_id","symbol","biotype","chromosome","start",
                                        "end","strand","release","assembly","ensembl_release"])
    if genes.stable_id.duplicated().any():
        raise ValueError("Duplicate GENCODE stable gene ID")
    return genes, records

def main():
    prepare()
    audits = []
    def record(path, **info):
        audits.append(dict(file=path.relative_to(ROOT).as_posix(), bytes=path.stat().st_size,
                           sha256=digest(path), **info))
    path = source("gencode_lncrna_gtf")
    genes, counts = parse_genes(path)
    genes.biotype.value_counts().rename_axis("gene_biotype").to_csv(OUT/"gencode_biotype_summary.csv",header=["genes"])
    genes.loc[genes.biotype.eq("lncRNA")].to_parquet(ANN/"lncrna_reference.parquet",index=False)
    record(path, source="GENCODE v50 / Ensembl 116", assembly="GRCh38", namespace="versioned ENSG",
           records=sum(counts.values()), feature_counts=dict(counts), sequence=False,
           columns="seqname source feature start end score strand frame attributes".split(),
           gene_missing_ids=int(genes.gene_id.eq("").sum()), duplicate_gene_ids=0,
           duplicated_symbol_rows=int(genes.symbol.duplicated(False).sum()))
    path = source("gencode_transcripts_fasta")
    ids, gene_ids, types = set(), set(), Counter()
    catalog = defaultdict(set)
    duplicate = length = expected = records = 0
    with gzip.open(path,"rt") as stream:
        for line in stream:
            if line.startswith(">"):
                if length != expected:
                    raise ValueError("FASTA header length disagrees with sequence")
                fields = line[1:].rstrip().split("|")
                if not fields[0].startswith("ENST") or not fields[1].startswith("ENSG"):
                    raise ValueError("Unrecognized transcript namespace")
                duplicate += fields[0] in ids
                ids.add(fields[0])
                gene_ids.add(fields[1])
                catalog[(stable(fields[1]),fields[5])].add(fields[7])
                types[fields[7]] += 1
                expected, length = int(fields[6]), 0
                records += 1
            else:
                length += len(line.strip())
    if length != expected or duplicate:
        raise ValueError("Invalid transcript FASTA")
    pd.DataFrame([(gene,symbol,"|".join(sorted(types))) for (gene,symbol),types in catalog.items()],
                 columns=["stable_id","symbol","transcript_types"]).to_parquet(ANN/"transcript_gene_catalog.parquet",index=False)
    record(path, source="GENCODE v50 transcript FASTA (registry)", assembly="GRCh38",
           columns="ENST ENSG HavanaTranscript HavanaGene transcript_name gene_name length transcript_type optional_regions".split(),
           namespace="versioned ENST/ENSG", sequence=True, records=records, missing_ids=0,
           duplicate_ids=duplicate, distinct_genes=len(gene_ids), transcript_types=dict(types))
    for ds in (*BLOOD,BRAIN):
        path = annotation(ds)
        frame = table(path)
        key = "GeneID" if ds==BRAIN else "ID"
        symbol = "GENE SYMBOL" if ds==BLOOD[2] else "Symbol"
        extra = {}
        if "nuID" in frame:
            paired = frame.loc[frame.nuID.ne("") & frame.SEQUENCE.ne("")]
            extra["paired_nuID_sequence"] = len(paired)
            extra["paired_exact_decodes"] = int(paired.nuID.map(decode_nuid).eq(paired.SEQUENCE.str.upper()).sum())
        if ds==BRAIN:
            linked=frame.EnsemblGeneID[frame.EnsemblGeneID.ne("")]
            extra["duplicated_ensembl_link_rows"]=int(linked.duplicated(False).sum())
        record(path, source=CFG_SOURCE[ds], assembly=("GRCh38.p13" if ds==BRAIN else
               "unspecified" if ds==BLOOD[2] else "NCBI build 36 version 3 (coordinate header)"),
               records=len(frame), columns=list(frame), missing_by_field=frame.eq("").sum().to_dict(),
               sequence="SEQUENCE" in frame, namespace=("Entrez / Ensembl" if ds==BRAIN else
               "nuID" if ds==BLOOD[2] else "ILMN"), duplicate_ids=int(frame[key].duplicated().sum()),
               duplicated_nonempty_symbol_rows=int(frame.loc[frame[symbol].ne(""),symbol].duplicated(False).sum()), **extra)
    (OUT/"annotation_sources.json").write_text(json.dumps(audits,indent=2),encoding="utf-8")

CFG_SOURCE=dict(zip((*BLOOD,BRAIN),("GEO GPL10558-50081","GEO GPL6947-13512","GEO GPL15988",
                                     "GEO NCBI Gene annotation, release date unspecified")))
if __name__=="__main__":
    main()


