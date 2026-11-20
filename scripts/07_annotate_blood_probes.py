"""Label-independent blood probe annotation using exhaustive local transcript alignment."""
from collections import defaultdict
import re
import pandas as pd
# Support both direct PowerShell invocation and python -m scripts.<module>.
if __package__ in (None, ""):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.utils.stage10 import (ANN, BLOOD, BRAIN, CFG, OUT, STATES, align, annotation,
                                 measured_platform, prepare, source, stable, subset, table)

def classify(exact, near, historical, lncs, sense, other, conflict=False):
    """Multiple isoforms of one gene represent a single gene-level target."""
    if conflict or len(near | historical)>1:
        return "AMBIGUOUS"
    if len(exact)==1 and exact==near and sense and exact<=lncs:
        return "VALIDATED_UNIQUE_LNCRNA"
    if len(exact)==1 and sense and exact<=other:
        return "NON_LNCRNA"
    return "PROVISIONAL" if near or historical else "UNMAPPED"

def main():
    prepare()
    ref = pd.read_parquet(ANN/"lncrna_reference.parquet")
    lncs = set(ref.stable_id)
    catalog = pd.read_parquet(ANN/"transcript_gene_catalog.parquet")
    symbols = catalog.groupby("symbol").stable_id.agg(set).to_dict()
    other = set(catalog.loc[catalog.transcript_types.str.contains(
        r"protein_coding|pseudogene|(?:^|\|)(?:miRNA|rRNA|snRNA|snoRNA|Mt_tRNA)(?:$|\|)"),"stable_id"])-lncs
    ncbi = table(annotation(BRAIN)).set_index("GeneID")
    if not ncbi.index.is_unique:
        raise ValueError("Duplicate NCBI GeneID in platform crosswalk")
    different_type = ncbi.Status.eq("active") & ncbi.GeneType.isin(
        ["protein-coding","pseudo","snRNA","snoRNA","rRNA","tRNA","scRNA"])
    disputed = set(ncbi.loc[different_type,"EnsemblGeneID"]) & lncs
    links = ncbi.EnsemblGeneID.to_dict()
    types = ncbi.GeneType.to_dict()
    platforms = {ds:measured_platform(ds) for ds in BLOOD}
    sequences = {seq for frame in platforms.values() for seq in frame.sequence.unique()
                 if re.fullmatch("[ACGT]{50}",seq)}
    align(sequences,source("gencode_transcripts_fasta"),OUT/"sequence_hits.tsv")
    hits = defaultdict(list)
    with (OUT/"sequence_hits.tsv").open() as stream:
        for line in stream:
            seq,gene,dist,tx,strand,biotype = line.rstrip().split("\t")
            hits[seq].append((stable(gene),int(dist),strand))
    results = {}
    for ds,frame in platforms.items():
        records = []
        for probe,rows in frame.groupby("ID",sort=False):
            # Preserve all source assignments even if a later platform release repeats an ID.
            def values(col):
                return {v for v in rows[col] if v} if col in rows else set()
            seqs = values("sequence")
            names = values("Symbol") | values("GENE SYMBOL")
            entrez = values("Entrez_Gene_ID")
            history = {g for n in names for g in symbols.get(n,())}
            history.update(links[e] for e in entrez if links.get(e))
            matches = [h for seq in seqs for h in hits[seq]]
            exact = {g for g,d,_ in matches if d==0}
            near = {g for g,_,_ in matches}
            sense = any(d==0 and strand=="+" for _,d,strand in matches)
            accessions = values("Accession") | values("GB_ACC") | values("RefSeq_ID") | values("Source_Reference_ID")
            candidates = (near|history)&lncs
            conflict = len(seqs)>1 or bool(candidates & disputed) or bool(candidates and (
                any(types.get(e) in {"pseudo","protein-coding"} for e in entrez) or
                any(re.match(r"[NX]M_",a) for a in accessions)))
            status = classify(exact,near,history,lncs,sense,other,conflict)
            join = lambda items: "|".join(sorted(items))
            records.append(dict(probe_id=probe,platform=CFG["datasets"][ds]["platform_id"],
                platform_rows=int(rows["_merge"].eq("both").sum()),symbol=join(names),entrez_id=join(entrez),
                accessions=join(accessions),sequence=join(seqs),historical_genes=join(history),
                exact_genes=join(exact),within_two_substitutions=join(near),sense_exact=sense,
                conflicting_sequence_or_annotation=conflict,disputed_ncbi_links=join(candidates & disputed),
                candidate_lncrna_genes=join(candidates),mapping_status=status,
                validated_gene=join(exact) if status==STATES[0] else ""))
        result = pd.DataFrame(records)
        result.to_csv(ANN/f"{ds}_probe_mapping.csv",index=False)
        results[ds]=result
        counts=result.mapping_status.value_counts().reindex(STATES,fill_value=0)
        counts.rename_axis("mapping_status").to_csv(OUT/f"{ds}_mapping_summary.csv",header=["measured_probes"])
        use=result.mapping_status.eq(STATES[0])
        if use.any():
            subset(ds,result.loc[use,"probe_id"].tolist(),"validated_lncrna")
        else:
            use=result.candidate_lncrna_genes.ne("") & result.mapping_status.eq("PROVISIONAL")
            subset(ds,result.loc[use,"probe_id"].tolist(),"provisional_lncrna_review")
    discovery,external = (results[ds] for ds in (BLOOD[0],BLOOD[2]))
    by_seq=discovery.groupby("sequence").probe_id.agg("|".join)
    validated=discovery[discovery.validated_gene.ne("")]
    compatible=set(zip(validated.sequence,validated.validated_gene))
    ext=external[["probe_id","sequence","validated_gene","mapping_status","candidate_lncrna_genes"]].copy()
    ext["discovery_identical_sequence_probes"]=ext.sequence.map(by_seq).fillna("")
    ext["shared_validated_gene"]=ext.validated_gene.isin(validated.validated_gene)
    ext["identical_validated_measurement"]=[pair in compatible for pair in zip(ext.sequence,ext.validated_gene)]
    ext.to_csv(OUT/"external_feature_compatibility.csv",index=False)

if __name__=="__main__":
    main()

