"""Stage 12 brain-only partial Pearson correlations, disk-backed global BH and local GO."""
import gc
import json
import logging
from pathlib import Path
import sys
import tempfile

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.special import stdtr
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

from scripts.utils.stage10 import ANN, BRAIN, CFG, PROC, ROOT, digest, matrix_path
from scripts.utils.stage11 import design, eligibility, validate_counts
from scripts.utils.stage12_go import enrichment

OUT = ROOT / "outputs/tables/stage12"
FIG = ROOT / "outputs/figures/stage12"
PRIOR = ROOT / "outputs/tables/stage11"
FITS = ("primary", "sensitivity")
LOG = logging.getLogger("stage12")

def check_axes(raw, frame, expected):
    if not raw.index.equals(frame.index) or not raw.index.is_unique:
        raise ValueError("GSM sample mismatch or duplicate")
    if not frame.columns.equals(expected) or not frame.columns.is_unique or not raw.columns.is_unique:
        raise ValueError("Feature order, uniqueness or annotation mismatch")
    if frame.isna().any().any():
        raise ValueError("Missing count")
    validate_counts(frame)
    if not np.array_equal(frame.to_numpy(), raw[frame.columns].to_numpy()):
        raise ValueError("Values differ from original counts")

def poscounts(counts):
    """Median log ratio over each sample's positive entries, zero log contribution otherwise."""
    x = counts.to_numpy()
    logx = np.zeros(x.shape, dtype=float)
    np.log(x, out=logx, where=x>0)
    logx -= logx.mean(axis=0)
    logsf = np.array([np.median(row[c>0]) for row,c in zip(logx,x)])
    if not np.isfinite(logsf).all():
        raise ValueError("Undefined whole-transcriptome size factor")
    return pd.Series(np.exp(logsf-logsf.mean()), index=counts.index)

def unit_residuals(values, x):
    if not np.isfinite(values).all() or not np.isfinite(x).all():
        raise ValueError("Missing or nonfinite input to residualization")
    q, _ = np.linalg.qr(x)
    residuals = values - q @ (q.T @ values)
    norms = np.linalg.norm(residuals,axis=0)
    valid = norms > 1e-12 * np.maximum(1,np.linalg.norm(values,axis=0))
    residuals[:,valid] /= norms[valid]
    residuals[:,~valid] = np.nan
    return residuals,valid

def correlation_blocks(left,right,df,block=32):
    if df<=0:
        raise ValueError("Insufficient residual degrees of freedom")
    for start in range(0,left.shape[1],block):
        r=np.clip(left[:,start:start+block].T @ right,-1,1)
        with np.errstate(divide="ignore",invalid="ignore"):
            t=np.abs(r)*np.sqrt(df/((1-r)*(1+r)))
        yield start,r,2*stdtr(df,-t)

def global_bh(p,q,chunk=500_000):
    """Sort once across the entire fit and carry the reverse minimum across chunks."""
    flat=p.reshape(-1)
    out=q.reshape(-1)
    order=np.argsort(flat,kind="quicksort")
    valid=int(np.isfinite(flat).sum())
    out[order[valid:]]=np.nan
    end=valid
    running=1.
    while end:
        start=max(0,end-chunk)
        idx=order[start:end]
        adjusted=flat[idx]*(valid/np.arange(start+1,end+1,dtype=float))
        adjusted=np.minimum.accumulate(np.minimum(adjusted,running)[::-1])[::-1]
        out[idx]=adjusted
        running=adjusted[0]
        end=start

def load_inputs():
    raw_path=ROOT/CFG["datasets"][BRAIN]["counts_path"]
    raw=pd.read_csv(raw_path,sep="\t",index_col="GeneID").T
    raw.columns=raw.columns.astype(str)
    if raw.isna().any().any() or not raw.columns.str.fullmatch(r"[1-9]\d*").all():
        raise ValueError("Invalid raw GeneID or count")
    validate_counts(raw)
    ann_path=ANN/f"{BRAIN}_gene_biotypes.csv"
    ann=pd.read_csv(ann_path,dtype=str,keep_default_na=False).set_index("GeneID")
    if not ann.index.equals(raw.columns):
        raise ValueError("Raw Entrez and Stage 10 annotation order differ")
    matrices=[]
    sources=[raw_path,ann_path]
    for category in ("lncRNA","protein_coding"):
        path=PROC/f"{BRAIN}_{category.lower()}_raw_counts.parquet"
        frame=pd.read_parquet(path)
        check_axes(raw,frame,ann.index[ann.category.eq(category)])
        matrices.append(frame)
        sources.append(path)
    if len(matrices[0].columns.intersection(matrices[1].columns)):
        raise ValueError("Overlapping feature types")
    ledger_path=PRIOR/f"{BRAIN}_samples.csv"
    meta_path=matrix_path(BRAIN).parent/"metadata.csv"
    ledger=pd.read_csv(ledger_path,index_col="sample_id")
    meta=pd.read_csv(meta_path,index_col="sample_id")
    if not all(raw.index.equals(i) for i in [meta.index,ledger.index]) or ledger.reported_sample_code.duplicated().any():
        raise ValueError("Sample or subject alignment failure")
    verified=eligibility(BRAIN,meta)
    for col in ["diagnosis","age","sex",*FITS]:
        pd.testing.assert_series_equal(ledger[col],verified[col])
    sources.extend([ledger_path,meta_path,ROOT/CFG["datasets"][BRAIN]["annotation_path"]])
    sources.extend(PRIOR/f"{BRAIN}_{fit}_size_factors.csv" for fit in FITS)
    pd.DataFrame([(p.relative_to(ROOT),digest(p)) for p in sources],columns=["input","sha256"]).to_csv(OUT/"sources.csv",index=False)
    return raw,ann,matrices,ledger

def fit_pairs(tmp,raw,frames,ledger,audit):
    outputs={}
    norm_tables=[]
    regressions=[]
    for fit in FITS:
        meta=ledger.loc[ledger[fit]]
        x=design(meta,BRAIN)
        sf=poscounts(raw.loc[meta.index])
        prior=pd.read_csv(PRIOR/f"{BRAIN}_{fit}_size_factors.csv",index_col="sample_id").size_factors
        if not sf.index.equals(prior.index) or not (prior>0).all():
            raise ValueError("Stage 11 size factor mismatch")
        norm=meta[["diagnosis","age","sex"]].assign(fit=fit,whole=sf,lnc_only=prior,
                flagged=meta.qc_reasons.notna(),log2_ratio=np.log2(sf/prior))
        norm_tables.append(norm)
        model=sm.OLS(norm.log2_ratio,x.assign(flagged=norm.flagged.astype(float))).fit(cov_type="HC3")
        ci=model.conf_int()
        for term in model.params.index:
            regressions.append((fit,term,model.params[term],*ci.loc[term],model.pvalues[term]))
        unit=[]
        for frame in frames:
            y=np.log2(1+frame.loc[meta.index].to_numpy()/sf.to_numpy()[:,None])
            u,valid=unit_residuals(y,x.to_numpy())
            unit.append(u)
            audit.loc[frame.columns,fit+"_variable_residual"]=valid
        shape=(unit[0].shape[1],unit[1].shape[1])
        r,p,q=[np.memmap(tmp/f"{fit}_{s}.bin",mode="w+",shape=shape,dtype="float64") for s in ("r","p","q")]
        df=len(meta)-x.shape[1]-1
        for start,block,prob in correlation_blocks(*unit,df):
            r[start:start+len(block)]=block
            p[start:start+len(block)]=prob
        LOG.info("%s: adjusting %d x %d pairs, df=%d",fit,*shape,df)
        global_bh(p,q)
        outputs[fit]=(r,p,q)
    pd.concat(norm_tables).to_csv(OUT/"size_factors.csv",index_label="sample_id")
    pd.DataFrame(regressions,columns=["fit","term","log2_ratio_beta","ci_low","ci_high","p_descriptive"]).to_csv(OUT/"factor_comparison.csv",index=False)
    return outputs

def summarize(out,lnc,coding,ann):
    """Write both fits together and derive summaries only after global correction."""
    rows=[]
    tops=[]
    members=[]
    coding_seen={s:np.zeros(len(coding),bool) for s in (*FITS,"stable")}
    ids=coding.to_numpy(dtype=np.int64)
    with pq.ParquetWriter(OUT/"all_pairs.parquet",pa.schema([
        ("lnc_entrez",pa.int64()),("coding_entrez",pa.int64()),
        *((f"{fit}_{k}",pa.float64()) for fit in FITS for k in ("r","p","q"))]),compression="zstd",use_dictionary=["lnc_entrez","coding_entrez"]) as writer:
        for start in range(0,len(lnc),24):
            stop=min(start+24,len(lnc))
            block={"lnc_entrez":np.repeat(lnc[start:stop].to_numpy(dtype=np.int64),len(ids)),
                   "coding_entrez":np.tile(ids,stop-start)}
            for fit in FITS:
                for k,m in zip(("r","p","q"),out[fit]):
                    block[f"{fit}_{k}"]=m[start:stop].ravel()
            writer.write_table(pa.table(block))
            for i in range(start,stop):
                gene=lnc[i]
                row={"lnc_entrez":gene,"symbol":ann.at[gene,"Symbol"]}
                a,_,qa=(v[i] for v in out["primary"])
                b,_,qb=(v[i] for v in out["sensitivity"])
                sig=qa<.05
                sens=qb<.05
                same=(a*b)>0
                stable=sig & sens & same
                for name,mask,r,q in [("primary",sig,a,qa),("sensitivity",sens,b,qb)]:
                    row[name+"_valid"]=int(np.isfinite(q).sum())
                    row[name+"_positive"]=int(np.sum(mask & (r>0)))
                    row[name+"_negative"]=int(np.sum(mask & (r<0)))
                    coding_seen[name] |= mask
                coding_seen["stable"] |= stable
                row.update(stable_positive=int(np.sum(stable & (a>0))),stable_negative=int(np.sum(stable & (a<0))),
                           lost=int(np.sum(sig & ~sens)),gained=int(np.sum(sens & ~sig)),
                           reversed_both_significant=int(np.sum(sig & sens & ~same)),
                           primary_significant_reversed=int(np.sum(sig & (a*b<0))))
                rows.append(row)
                for direction,sign in [("positive",1),("negative",-1)]:
                    ix=np.flatnonzero(stable & (sign*a>0))
                    strongest=ix[np.argsort(-sign*a[ix],kind="stable")[:3]]
                    for j in strongest:
                        tops.append((gene,ann.at[gene,"Symbol"],coding[j],ann.at[coding[j],"Symbol"],
                                     direction,a[j],qa[j],b[j],qb[j]))
                    strong=ix[(sign*a[ix]>=.75)&(sign*b[ix]>=.75)]
                    members.extend((gene,direction,str(coding[j])) for j in strong)
    counts=pd.DataFrame(rows)
    counts.to_csv(OUT/"edges_per_lncrna.csv",index=False)
    edges=pd.DataFrame(tops,columns=["lnc_entrez","lnc_symbol","coding_entrez","coding_symbol","direction","primary_r","primary_q","sensitivity_r","sensitivity_q"])
    edges.to_csv(OUT/"edge_summary.csv",index=False)
    foreground=pd.DataFrame(members,columns=["lnc_entrez","direction","coding_entrez"])
    foreground.to_csv(OUT/"enrichment_foreground.csv.gz",index=False)
    return edges,foreground,{k:int(v.sum()) for k,v in coding_seen.items()}

def plots(edges):
    norm=pd.read_csv(OUT/"size_factors.csv").query("fit=='primary'")
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    for group,c in [("Control","#226979"),("AD","#bc5932")]:
        s=norm[norm.diagnosis.eq(group)]
        ax[0].scatter(np.log2(s.lnc_only),np.log2(s.whole),s=12,alpha=.65,c=c,label=group)
    bounds=np.log2(norm[["whole","lnc_only"]].to_numpy())
    ax[0].plot([bounds.min(),bounds.max()],[bounds.min(),bounds.max()],"k--",lw=1)
    ax[0].set(xlabel="Stage 11 lncRNA-only factor (log2)",ylabel="Whole-transcriptome factor (log2)")
    ax[0].legend(frameon=False)
    labels=[]; groups=[]
    for group in ["Control","AD"]:
        for flag in [False,True]:
            groups.append(norm.loc[norm.diagnosis.eq(group)&norm.flagged.eq(flag),"log2_ratio"])
            labels.append(group+"\n"+("flagged" if flag else "unflagged"))
    ax[1].boxplot(groups,tick_labels=labels)
    ax[1].axhline(0,color="0.5",lw=.7)
    ax[1].set_ylabel("log2(whole / lncRNA-only size factor)")
    fig.suptitle("Normalization references in fusiform gyrus · primary samples")
    fig.savefig(FIG/"normalization_comparison.png",dpi=160); plt.close(fig)
    if edges.empty:
        return
    view=edges.assign(strength=np.minimum(edges.primary_r.abs(),edges.sensitivity_r.abs())).sort_values("strength",ascending=False)
    genes=view.lnc_entrez.drop_duplicates().head(4)
    view=view[view.lnc_entrez.isin(genes)].groupby("lnc_entrez",sort=False).head(2)
    fig,ax=plt.subplots(figsize=(10,5),layout="constrained")
    left={g:i for i,g in enumerate(genes)}
    partners=view.coding_entrez.drop_duplicates()
    right=dict(zip(partners,np.linspace(0,max(1,len(genes)-1),len(partners))))
    for r in view.itertuples():
        ax.plot([0,1],[left[r.lnc_entrez],right[r.coding_entrez]],lw=2,color="#247677" if r.primary_r>0 else "#b45137")
    for x,positions,idcol,symcol in [(0,left,"lnc_entrez","lnc_symbol"),(1,right,"coding_entrez","coding_symbol")]:
        for g,y in positions.items():
            ax.scatter(x,y,color="#283d4a",s=32,zorder=2)
            sym=view.loc[view[idcol].eq(g),symcol].iloc[0]
            ax.text(x+(.025 if x else -.025),y,sym,va="center",ha="left" if x else "right")
    ax.set(xlim=(-.65,1.65),ylim=(-.55,max(1,len(genes)-1)+.55),title="Selected reproducible associations · lncRNA (left) and coding gene (right)")
    ax.text(.5,-.43,"Partial correlations: teal positive · rust negative",ha="center",fontsize=9)
    ax.axis("off")
    fig.savefig(FIG/"coexpression_network.png",dpi=160); plt.close(fig)
    view.to_csv(OUT/"network_edges.csv",index=False)

def run(tmp):
    raw,ann,frames,ledger=load_inputs()
    ledger.to_csv(OUT/"sample_audit.csv")
    audits=[]
    for frame in frames:
        primary=frame.loc[ledger.primary]
        count=primary.ge(10).sum()
        constant=primary.nunique().eq(1)
        audits.append(pd.DataFrame({"n_ge10":count,"constant_raw":constant,
                                   "eligible":count.ge(10)&~constant,"missing":frame.isna().sum()}))
    audit=pd.concat(audits).join(ann[["Symbol","category","evidence"]])
    frames=[f.loc[:,audit.loc[f.columns,"eligible"].to_numpy()] for f in frames]
    out=fit_pairs(tmp,raw,frames,ledger,audit)
    audit.to_csv(OUT/"feature_audit.csv",index_label="entrez")
    edges,front,seen=summarize(out,frames[0].columns,frames[1].columns,ann)
    background=audit.index[audit.category.eq("protein_coding") & audit.primary_variable_residual.eq(True)]
    go=enrichment(front,background,ann,OUT)
    plots(edges)
    (OUT/"run_summary.json").write_text(json.dumps({"distinct_coding_with_edges":seen,"GO":go},indent=2))
    LOG.info("Completed both fits and local GO analysis")

def main():
    OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=OUT) as tmp, threadpool_limits(limits=2):
        run(Path(tmp))
        gc.collect()

if __name__=="__main__":
    main()

