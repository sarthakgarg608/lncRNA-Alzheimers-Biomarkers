"""Summarize reference-supported targets and technical overlap, independent of outcomes."""
from itertools import combinations
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Support both direct PowerShell invocation and python -m scripts.<module>.
if __package__ in (None, ""):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.utils.stage10 import ANN,BLOOD,BRAIN,FIG,OUT,PROC,STATES,prepare

def chart(data,title,ylabel,name):
    fig,ax=plt.subplots(figsize=(9,4.7),layout="constrained")
    data.plot.bar(ax=ax,rot=0,color=(["#20716d","#df9960","#7173a3","#8495a1","#aec8bb"]
                                  if isinstance(data,pd.DataFrame) else "#20716d"))
    ax.set_title(title,loc="left")
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.spines[["right","top"]].set_visible(False)
    if isinstance(data,pd.DataFrame):
        ax.legend(fontsize=8,frameon=False)
    fig.savefig(FIG/f"{name}.png",dpi=150)
    plt.close(fig)

def main():
    prepare()
    inventory=[]
    coverage=[]
    concordance=[]
    pairs={}
    genes={}
    states={}
    for ds in BLOOD:
        m=pd.read_csv(ANN/f"{ds}_probe_mapping.csv",keep_default_na=False)
        states[ds]=m.mapping_status.value_counts().reindex(STATES,fill_value=0)
        v=m[m.mapping_status.eq(STATES[0])].copy()
        mat=pd.read_parquet(PROC/f"{ds}_validated_lncrna.parquet") if len(v) else pd.DataFrame()
        variation=mat.nunique()>1
        for gene,group in v.groupby("validated_gene"):
            for a,b in combinations(group.probe_id,2):
                concordance.append((ds,gene,a,b,mat[a].corr(mat[b])))
        for r in v.itertuples():
            inventory.append((ds,r.probe_id,"ILMN" if ds!=BLOOD[2] else "nuID",r.validated_gene,
                              r.sequence,r.mapping_status,variation[r.probe_id]))
        pairs[ds]=set(zip(v.validated_gene,v.sequence))
        genes[ds]=set(v.validated_gene)
        coverage.append(dict(dataset=ds,measured=len(m),platform_matches=int(m.platform_rows.gt(0).sum()),
                             candidate_probes=int(m.candidate_lncrna_genes.ne("").sum()),validated_probes=len(v),
                             validated_genes=len(genes[ds]),variable_probes=int(variation.sum()),
                             ambiguous=int(states[ds]["AMBIGUOUS"]),unmapped=int(states[ds]["UNMAPPED"])))
    brain=pd.read_csv(ANN/f"{BRAIN}_gene_biotypes.csv",dtype=str,keep_default_na=False)
    lnc=brain[brain.category.eq("lncRNA")]
    genes[BRAIN]=set(lnc.EnsemblGeneID)-{""}
    counts=pd.read_parquet(PROC/f"{BRAIN}_lncrna_raw_counts.parquet")
    variation=counts.nunique()>1
    for r in lnc.itertuples():
        inventory.append((BRAIN,r.GeneID,"Entrez",r.EnsemblGeneID,"",r.evidence,variation[r.GeneID]))
    pd.DataFrame(inventory,columns=["dataset","feature_id","namespace","ensembl_id","sequence","evidence",
                                   "nonconstant"]).to_csv(OUT/"validated_feature_inventory.csv",index=False)
    pd.DataFrame(concordance,columns=["dataset","ensembl_id","probe_a","probe_b","pearson_r"]).to_csv(
        OUT/"within_gene_probe_concordance.csv",index=False)
    cover=pd.DataFrame(coverage).set_index("dataset")
    cover.to_csv(OUT/"blood_coverage.csv")
    rows=[]
    for size in range(1,5):
        for group in combinations((*BLOOD,BRAIN),size):
            shared=set.intersection(*(genes[ds] for ds in group))
            same=set.intersection(*(pairs[ds] for ds in group)) if BRAIN not in group else None
            rows.append(dict(datasets=" + ".join(group),shared_linked_ensembl_genes=len(shared),
                             identical_sequence_gene_pairs=len(same) if same is not None else None,
                             genes_with_identical_sequence=len({g for g,_ in same}) if same is not None else None,
                             brain_entrez_lncrnas=len(lnc) if group==(BRAIN,) else None))
    overlap=pd.DataFrame(rows)
    overlap.to_csv(OUT/"lncrna_overlap_summary.csv",index=False)
    chart(pd.DataFrame(states).T.rename(columns={
        STATES[0]:"Unique lncRNA","PROVISIONAL":"Provisional","AMBIGUOUS":"Ambiguous",
        "NON_LNCRNA":"Other RNA","UNMAPPED":"Unmapped"}),"Blood probe assignments","Measured probes","mapping_confidence")
    chart(cover[["candidate_probes","validated_probes"]].rename(columns={"candidate_probes":"Any lncRNA evidence",
          "validated_probes":"Validated"}),"Blood probes with lncRNA evidence","Measured probes","candidate_vs_validated")
    chart(cover.ambiguous,"Probes with conflicting assignments or multiple targets","Measured probes","ambiguity")
    chart(brain.category.value_counts(),"Gene classification in fusiform gyrus RNA-seq","Measured Entrez genes","brain_biotypes")
    summary=pd.read_csv(OUT/"gencode_biotype_summary.csv").set_index("gene_biotype")
    chart(summary.genes,"Gene types in the GENCODE v50 long non-coding RNA GTF","Gene records","gencode_biotypes")
    view=overlap[overlap.datasets.isin([" + ".join(g) for n in (2,3) for g in combinations(BLOOD,n)])].copy()
    view.index=view.datasets.str.replace("GSE","").str.replace(" + ","\n",regex=False)
    chart(view[["shared_linked_ensembl_genes","genes_with_identical_sequence"]].rename(columns={
        "shared_linked_ensembl_genes":"Same gene","genes_with_identical_sequence":"Same gene and probe sequence"}),
          "Validated blood features across cohorts","Genes","shared_blood_genes")

if __name__=="__main__":
    main()

