"""Overrepresentation against the local NCBI Entrez/GO Biological Process snapshot."""
from collections import defaultdict
import re

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from .stage10 import BRAIN, annotation, digest, table
from .stage11 import bh

def enrichment(foreground,background,ann,out):
    path=annotation(BRAIN)
    ref=table(path).set_index("GeneID")
    columns=["GOProcessID","GOProcess"]
    pd.testing.assert_frame_equal(ref.loc[background,columns],ann.loc[background,columns],check_names=False)
    sets=defaultdict(set)
    names=defaultdict(set)
    annotated=set()
    for gene,(ids,labels) in zip(background,ref.loc[background,columns].itertuples(index=False,name=None)):
        if not ids:
            continue
        terms=ids.split("///")
        labels=labels.split("///")
        if len(terms)!=len(labels) or any(not re.fullmatch(r"GO:\d{7}",term) for term in terms):
            raise ValueError(f"Unresolved local GO ID/term mapping for {gene}")
        annotated.add(gene)
        for term,label in zip(terms,labels):
            sets[term].add(gene)
            names[term].add(label)
    pd.DataFrame({"coding_entrez":background,"symbol":ann.loc[background,"Symbol"].to_numpy(),
                  "has_GO_BP":background.isin(annotated)}).to_csv(out/"enrichment_background.csv",index=False)
    terms={term:genes for term,genes in sets.items() if 15<=len(genes)<=500}
    reverse=defaultdict(list)
    for index,genes in enumerate(terms.values()):
        for gene in genes:
            reverse[gene].append(index)
    universe=len(background)
    sizes=np.array([len(genes) for genes in terms.values()],dtype=int)
    ids=list(terms)
    common=pd.DataFrame({"go_id":ids,"term":[" | ".join(sorted(names[t])) for t in ids],"term_size":sizes})
    parts=[]
    lists=[]
    for (lnc,sign),group in foreground.groupby(["lnc_entrez","direction"]):
        genes=set(group.coding_entrez)
        if not genes.issubset(background):
            raise ValueError("Foreground outside measured coding universe")
        lists.append((lnc,sign,len(genes),len(genes & annotated),len(genes.intersection(reverse))))
        if len(genes)<10 or not ids:
            continue
        hits=np.bincount([idx for g in genes for idx in reverse.get(g,[])],minlength=len(ids))
        p=np.ones(len(ids))
        positive=hits>0
        p[positive]=hypergeom.sf(hits[positive]-1,universe,sizes[positive],len(genes))
        parts.append(common.assign(lnc_entrez=lnc,direction=sign,foreground_size=len(genes),
                                   overlap=hits,fold_enrichment=hits*universe/(sizes*len(genes)),p=p))
    pd.DataFrame(lists,columns=["lnc_entrez","direction","genes","mapped_any_BP","mapped_tested_BP"]).to_csv(out/"foreground_coverage.csv",index=False)
    cols=["go_id","term","term_size","lnc_entrez","direction","foreground_size","overlap","fold_enrichment","p"]
    tests=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=cols)
    tests["q"]=bh(tests.p)
    tests.to_parquet(out/"GO_BP_all_tests.parquet",compression="zstd",index=False)
    hits=tests[tests.q.lt(.05)].sort_values(["q","p"])
    hits.to_csv(out/"GO_BP_significant.csv.gz",index=False)
    return dict(source=path.name,sha256=digest(path),background=universe,annotated=len(annotated),
                distinct_terms=len(sets),eligible_terms=len(terms),eligible_lists=len(parts),
                tested=len(tests),significant=len(hits),distinct_significant_terms=hits.go_id.nunique())

