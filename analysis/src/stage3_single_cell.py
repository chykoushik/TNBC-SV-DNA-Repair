import itertools
import json
import tarfile
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from .config import PROJECT_ROOT,CANDIDATE_GENES
from .statistics import bh
from .utils import sha256

CLASS_MAP={'Cancer Epithelial':'malignant_epithelial','Normal Epithelial':'normal_epithelial',
           'T-cells':'immune','B-cells':'immune','Plasmablasts':'immune','Myeloid':'immune',
           'CAFs':'stromal','PVL':'stromal','Endothelial':'stromal'}

def accumulate_triplets(triples,cell_group,gene_lookup,totals,sums):
    rows=triples[:,0]-1;cells=triples[:,1]-1;values=triples[:,2]
    if (rows<0).any() or (rows>=len(gene_lookup)).any() or (cells<0).any() or (cells>=len(cell_group)).any() or (values<0).any():raise ValueError('Invalid sparse coordinate/count')
    groups=cell_group[cells]
    totals+=np.bincount(groups,weights=values,minlength=len(totals))
    selected=gene_lookup[rows];keep=selected>=0
    packed=selected[keep]*len(totals)+groups[keep]
    sums+=np.bincount(packed,weights=values[keep],minlength=sums.size).reshape(sums.shape)

def pseudobulk(root=PROJECT_ROOT):
    original=root.parent.parent
    source=original/'dataset/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'
    loose=original/'results/sc_extracted/Wu_etal_2021_BRCA_scRNASeq'
    dest=root/'outputs/v2/stage3';dest.mkdir(parents=True,exist_ok=True)
    cache=dest/'single_cell_pseudobulk_counts.csv';manifest=dest/'single_cell_source.json'
    stamp={str(p):{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in [source,loose/'metadata.csv',loose/'count_matrix_barcodes.tsv',loose/'count_matrix_genes.tsv']}
    if cache.exists() and manifest.exists() and json.loads(manifest.read_text())['source_stat']==stamp:
        return pd.read_csv(cache)
    meta=pd.read_csv(loose/'metadata.csv').rename(columns={'Unnamed: 0':'barcode'})
    barcodes=pd.read_csv(loose/'count_matrix_barcodes.tsv',header=None,sep='\t')[0].tolist()
    genes=pd.read_csv(loose/'count_matrix_genes.tsv',header=None,sep='\t')[0].tolist()
    if len(set(barcodes))!=len(barcodes) or meta.barcode.duplicated().any():raise ValueError('Duplicate cell IDs')
    if set(barcodes)!=set(meta.barcode):raise ValueError('Cell metadata/barcode mismatch')
    meta=meta.set_index('barcode').loc[barcodes].reset_index()
    if not all(str(b).startswith(str(d)+'_') for b,d in zip(meta.barcode,meta['orig.ident'])):raise ValueError('Barcode/donor prefix mismatch')
    meta['cell_class']=meta.celltype_major.map(CLASS_MAP)
    if meta.cell_class.isna().any():raise ValueError('Unmapped cell class')
    if (meta.groupby('orig.ident').subtype.nunique()>1).any():raise ValueError('Donor subtype conflict')
    keys=list(zip(meta['orig.ident'],meta.cell_class));unique=sorted(set(keys));lookup={k:i for i,k in enumerate(unique)}
    cell_group=np.array([lookup[k] for k in keys],dtype=np.int64)
    gene_lookup=np.full(len(genes),-1,dtype=np.int64)
    for i,g in enumerate(CANDIDATE_GENES):
        positions=[j for j,s in enumerate(genes) if s==g]
        if len(positions)!=1:raise ValueError('Missing/ambiguous single-cell candidate: '+g)
        gene_lookup[positions[0]]=i
    totals=np.zeros(len(unique));sums=np.zeros((9,len(unique)));observed=0
    with tarfile.open(source,'r|gz') as archive:
        for member in archive:
            if not member.name.endswith('/count_matrix_sparse.mtx'):continue
            stream=archive.extractfile(member)
            if stream.readline().strip()!=b'%%MatrixMarket matrix coordinate integer general':raise ValueError('Unexpected Matrix Market encoding')
            line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            ng,nc,nnz=map(int,line.split())
            if (ng,nc)!=(len(genes),len(barcodes)):raise ValueError('Sparse dimensions do not match annotation')
            while True:
                chunk=b''.join(itertools.islice(stream,200000))
                if not chunk:break
                values=np.fromstring(chunk.decode('ascii'),sep=' ',dtype=np.int64)
                if values.size%3:raise ValueError('Malformed sparse triplets')
                triples=values.reshape(-1,3);accumulate_triplets(triples,cell_group,gene_lookup,totals,sums);observed+=len(triples)
                if observed%10000000==0:print(f'Single-cell pseudobulk: {observed:,}/{nnz:,} entries',flush=True)
            if observed!=nnz:raise ValueError('Truncated sparse matrix')
            break
        else:raise ValueError('Sparse matrix not found in source archive')
    rows=[]
    for i,(donor,cls) in enumerate(unique):
        subset=meta[(meta['orig.ident']==donor)&(meta.cell_class==cls)]
        r=dict(donor_sample_id=donor,cell_class=cls,subtype=subset.subtype.iloc[0],n_cells=len(subset),library_counts=totals[i])
        r.update(dict(zip(CANDIDATE_GENES,sums[:,i])));rows.append(r)
    result=pd.DataFrame(rows);result.to_csv(cache,index=False)
    manifest.write_text(json.dumps({'source_stat':stamp,'archive_sha256':sha256(source),'nnz':observed,'cells':len(meta),'donor_sample_units':len(meta['orig.ident'].unique()),
                                    'grouping':'Exact orig.ident units; suffixes retained; barcode prefix checked','archive_member':member.name},indent=2))
    return result

def paired_result(differences,label,gene='composite',population='TNBC',seed=20260908):
    d=np.asarray(differences,dtype=float);d=d[np.isfinite(d)]
    r=dict(component='single_cell',population=population,endpoint=label,gene=gene,N=len(d),effect=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,status='unavailable',reason='fewer than 5 paired donor/sample units')
    if len(d)<5:return r
    rng=np.random.default_rng(seed);boot=np.mean(d[rng.integers(0,len(d),size=(2000,len(d)))],axis=1)
    r.update(effect=float(d.mean()),ci_low=float(np.quantile(boot,.025)),ci_high=float(np.quantile(boot,.975)),
             p=float(wilcoxon(d,alternative='two-sided',method='auto').pvalue) if np.any(d!=0) else 1.,status='estimated',reason='')
    return r

def analyze(root=PROJECT_ROOT):
    raw=pseudobulk(root);dest=root/'outputs/v2/stage3'
    totals=[]
    for donor,sub in raw.groupby('donor_sample_id'):
        nm=sub[sub.cell_class!='malignant_epithelial'];r=nm[['n_cells','library_counts',*CANDIDATE_GENES]].sum().to_dict()
        r.update(donor_sample_id=donor,cell_class='nonmalignant_combined',subtype=sub.subtype.iloc[0]);totals.append(r)
    frame=pd.concat([raw,pd.DataFrame(totals)],ignore_index=True)
    # Library count consistency
    metadata=pd.read_csv(root.parent.parent/'results/sc_extracted/Wu_etal_2021_BRCA_scRNASeq/metadata.csv')
    metadata['cell_class']=metadata.celltype_major.map(CLASS_MAP)
    expected=metadata.groupby(['orig.ident','cell_class']).nCount_RNA.sum()
    actual=raw.set_index(['donor_sample_id','cell_class']).library_counts
    if set(expected.index)!=set(actual.index) or not np.allclose(expected.sort_index(),actual.sort_index(),rtol=0,atol=0):
        raise ValueError('Pseudobulk library totals disagree with source nCount_RNA')
    usable=(frame.n_cells>=20)&(frame.library_counts>0)
    log=np.log2(frame[list(CANDIDATE_GENES)].div(frame.library_counts,axis=0)*1e6+1)
    reference=usable&(frame.subtype=='TNBC')&(frame.cell_class!='nonmalignant_combined')
    mean=log.loc[reference].mean();sd=log.loc[reference].std(ddof=0)
    if (sd<=0).any():raise ValueError('Constant candidate in TNBC single-cell pseudobulk reference')
    z=log.sub(mean).div(sd);frame['candidate_mean_z']=z.mean(axis=1).where(usable)
    for gene in CANDIDATE_GENES:frame[gene+'_log2cpm1']=log[gene];frame[gene+'_z']=z[gene].where(usable)
    frame['eligible']=usable;frame.to_csv(dest/'single_cell_pseudobulk.csv',index=False)
    results=[]
    for pop,subset in [('TNBC',frame[frame.subtype=='TNBC']),('all_breast',frame)]:
        for cls in ['nonmalignant_combined','immune','stromal','normal_epithelial']:
            pivot=subset.pivot(index='donor_sample_id',columns='cell_class',values='candidate_mean_z')
            if cls not in pivot or 'malignant_epithelial' not in pivot:diff=[]
            else:diff=(pivot.malignant_epithelial-pivot[cls]).dropna()
            results.append(paired_result(diff,'malignant_vs_'+cls,population=pop))
        for gene in CANDIDATE_GENES:
            pivot=subset[subset.eligible].pivot(index='donor_sample_id',columns='cell_class',values=gene+'_log2cpm1')
            diff=(pivot.malignant_epithelial-pivot.nonmalignant_combined).dropna()
            results.append(paired_result(diff,'malignant_vs_nonmalignant_combined',gene=gene,population=pop))
    for pop in ['TNBC','all_breast']:
        for composite in [True,False]:
            subset=[r for r in results if r['population']==pop and (r['gene']=='composite')==composite]
            for r,q in zip(subset,bh([r['p'] for r in subset])):r['fdr']=q
    pd.DataFrame(results).to_csv(dest/'single_cell_results.csv',index=False)
    metadata={'all_donor_sample_n':raw.donor_sample_id.nunique(),'tnbc_donor_sample_n':raw[raw.subtype=='TNBC'].donor_sample_id.nunique(),
              'minimum_cells':20,'reference':'TNBC donor-by-disjoint-cell-class pseudobulks','gene_means':mean.to_dict(),'gene_sd':sd.to_dict(),
              'primary_result':results[0], 'library_count_validation':'All 85 donor/class totals equal source nCount_RNA sums exactly',
              'donor_identity_caveat':'orig.ident is the published donor/sample unit; no suffix stripping or independent cells-as-patients tests'}
    (dest/'single_cell_summary.json').write_text(json.dumps(metadata,indent=2,default=lambda x:int(x) if isinstance(x,np.integer) else x))
    return results,frame,metadata

if __name__=='__main__':pseudobulk()
