import csv
import gzip
import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
from statsmodels.duration.hazard_regression import PHReg
from .config import PROJECT_ROOT,CANDIDATE_GENES
from .covariates import PROLIFERATION_GENES
from .statistics import fit_model,bh

def gene_mean(frame):
    return frame.groupby('Hugo_Symbol',sort=True).mean(numeric_only=True)

def signature(expression,genes):
    missing=[g for g in genes if g not in expression.index]
    if missing:return pd.Series(dtype=float),{},missing
    x=expression.loc[list(genes)];sd=x.std(axis=1,ddof=0)
    if (sd<=0).any() or x.isna().any().any():return pd.Series(dtype=float),{},[g for g in genes if sd[g]<=0 or x.loc[g].isna().any()]
    z=x.sub(x.mean(axis=1),axis=0).div(sd,axis=0)
    return z.mean(axis=0),{'mean':x.mean(axis=1).to_dict(),'sd':sd.to_dict(),'genes':list(genes)},[]

def cox(data,time,event,predictor='candidate_perSD',label='OS',covariates=()):
    d=data[[time,event,predictor,*covariates]].replace([np.inf,-np.inf],np.nan).dropna()
    d=d[(d[time]>0)&d[event].isin([0,1])]
    r=dict(component='patient',endpoint=label,gene='composite',N=len(d),events=int(d[event].sum()),status='unavailable',reason='insufficient events/cases',beta=np.nan,se=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,hazard_ratio=np.nan,hr_ci_low=np.nan,hr_ci_high=np.nan)
    if len(d)<20 or d[event].sum()<max(10,10*(1+len(covariates))):return r
    x=d[[predictor,*covariates]].astype(float)
    if np.linalg.matrix_rank(x)!=len(x.columns):r['reason']='rank deficiency';return r
    try:
        with warnings.catch_warnings(record=True) as ws:
            fit=PHReg(d[time],x,status=d[event],ties='efron').fit()
        if ws or not np.isfinite(fit.params).all():r['reason']='Cox convergence/finite estimate problem';return r
        lo,hi=fit.conf_int()[0]
        r.update(status='estimated',reason='',beta=float(fit.params[0]),se=float(fit.bse[0]),ci_low=float(lo),ci_high=float(hi),p=float(fit.pvalues[0]),
                 hazard_ratio=float(np.exp(fit.params[0])),hr_ci_low=float(np.exp(lo)),hr_ci_high=float(np.exp(hi)))
        residual=np.asarray(fit.schoenfeld_residuals)[:,0];valid=(d[event].to_numpy()==1)&np.isfinite(residual)
        r['ph_diagnostic_log_time_spearman_p']=float(st.spearmanr(residual[valid],np.log(d.loc[valid,time])).pvalue) if valid.sum()>=5 else np.nan
    except (ValueError,np.linalg.LinAlgError) as exc:r['reason']=str(exc)
    return r

def metabric(root=PROJECT_ROOT):
    original=root.parent.parent;folder=original/'results/metabric_extracted/brca_metabric';dest=root/'outputs/v2/stage3'
    clinical=pd.read_csv(original/'dataset/brca_metabric_clinical_data.tsv',sep='\t')
    normalize=lambda s:s.astype(str).str.strip().str.lower().replace({'positve':'positive'})
    er=normalize(clinical['ER status measured by IHC']);reported=normalize(clinical['ER Status'])
    pr=normalize(clinical['PR Status']);her2=normalize(clinical['HER2 Status'])
    strict=(er=='negative')&(reported!='positive')&(pr=='negative')&(her2=='negative')&(clinical['Sample Type']=='Primary')
    clinical['clinical_TNBC']=strict
    clinical[['Patient ID','Sample ID','ER status measured by IHC','ER Status','PR Status','HER2 Status','clinical_TNBC']].to_csv(dest/'metabric_cohort_decisions.csv',index=False)
    cohort=clinical.loc[strict].copy().sort_values(['Patient ID','Sample ID'])
    if cohort['Patient ID'].duplicated().any():raise ValueError('METABRIC duplicate patients need source/QC adjudication')
    epath=folder/'data_mrna_illumina_microarray.txt'
    header=pd.read_csv(epath,sep='\t',nrows=0)
    samples=sorted(set(cohort['Sample ID'])&set(header.columns))
    e=pd.read_csv(epath,sep='\t',usecols=['Hugo_Symbol',*samples])
    gene_rows=e[['Hugo_Symbol',*samples]]
    duplicate_n=int((gene_rows.Hugo_Symbol.value_counts()>1).sum())
    expression=gene_mean(gene_rows);score,parameters,missing=signature(expression,CANDIDATE_GENES)
    cohort['candidate_mean_z']=cohort['Sample ID'].map(score);cohort['candidate_perSD']=cohort.candidate_mean_z/cohort.candidate_mean_z.std(ddof=0)
    prolif,_,pmissing=signature(expression,PROLIFERATION_GENES);cohort['proliferation']=cohort['Sample ID'].map(prolif)
    for g in CANDIDATE_GENES:
        if g in expression.index:cohort[g+'_z']=cohort['Sample ID'].map((expression.loc[g]-expression.loc[g].mean())/expression.loc[g].std(ddof=0))
    cpath=folder/'data_cna.txt';chead=pd.read_csv(cpath,sep='\t',nrows=0);cn_samples=sorted(set(samples)&set(chead.columns))
    cna=pd.read_csv(cpath,sep='\t',usecols=['Hugo_Symbol',*cn_samples])
    # Equal gene contributions
    dna=cna[['Hugo_Symbol',*cn_samples]].groupby('Hugo_Symbol',sort=True).median(numeric_only=True)
    valid=dna.notna();altered=(dna.abs()>=1)&valid
    fraction=altered.sum()/valid.sum();cohort['altered_gene_fraction']=cohort['Sample ID'].map(fraction)
    cohort['age']=pd.to_numeric(cohort['Age at Diagnosis'],errors='coerce')
    cohort['OS_months']=pd.to_numeric(cohort['Overall Survival (Months)'],errors='coerce')
    cohort['OS_event']=cohort['Overall Survival Status'].map({'0:LIVING':0,'1:DECEASED':1})
    cohort['RFS_months']=pd.to_numeric(cohort['Relapse Free Status (Months)'],errors='coerce')
    cohort['RFS_event']=cohort['Relapse Free Status'].map({'0:Not Recurred':0,'1:Recurred':1})
    results=[]
    r,_,_=fit_model(cohort,outcome='altered_gene_fraction',label='METABRIC independent CNA proxy')
    r.update(component='METABRIC',endpoint='altered_gene_fraction',gene='composite',validation_type='partial genomic replication; not HRD scar');results.append(r)
    for time,event,label in [('OS_months','OS_event','OS'),('RFS_months','RFS_event','RFS')]:
        r=cox(cohort,time,event,label=label);r['component']='METABRIC';r['validation_type']='secondary clinical outcome';results.append(r)
    for r,q in zip(results,bh([r['p'] for r in results])):r['fdr']=q
    r,_,_=fit_model(cohort,outcome='altered_gene_fraction',covariates=('age','proliferation'),label='CNA proxy adjusted for age/proliferation')
    r.update(component='METABRIC',endpoint='altered_gene_fraction_adjusted',gene='composite',fdr=np.nan,validation_type='prespecified sensitivity');results.append(r)
    r=dict(component='METABRIC',endpoint='HRD_scar',gene='composite',N=0,status='unavailable',reason='No independent HRD scar table local',beta=np.nan,p=np.nan,fdr=np.nan);results.append(r)
    for endpoint in ['altered_gene_fraction','OS','RFS']:
        exploratory=[]
        for gene in CANDIDATE_GENES:
            if endpoint=='altered_gene_fraction':
                r,_,_=fit_model(cohort,outcome=endpoint,predictor=gene+'_z',label='Exploratory individual gene CNA')
            else:
                r=cox(cohort,endpoint+'_months',endpoint+'_event',predictor=gene+'_z',label=endpoint)
            r.update(component='METABRIC',endpoint=endpoint,gene=gene,validation_type='exploratory gene; not weight optimization');exploratory.append(r)
        for r,q in zip(exploratory,bh([r['p'] for r in exploratory])):r['fdr']=q
        results.extend(exploratory)
    cohort.to_csv(dest/'metabric_analysis_cohort.csv',index=False)
    pd.DataFrame(results).to_csv(dest/'metabric_results.csv',index=False)
    summary=dict(clinical_tnbc_n=len(cohort),rna_n=int(cohort.candidate_mean_z.notna().sum()),candidate_genes_available=[g for g in CANDIDATE_GENES if g not in missing],missing_genes=missing,
                 source_reported_three_negative_n=int(((reported=='negative')&(pr=='negative')&(her2=='negative')).sum()),
                 reported_negative_but_IHC_positive_n=int(((reported=='negative')&(pr=='negative')&(her2=='negative')&(er=='positive')).sum()),
                 duplicate_gene_symbols_averaged=duplicate_n,STAG3_source_rows=int((e.Hugo_Symbol=='STAG3').sum()),
                 signature_parameters=parameters,proliferation_missing=pmissing,expression_units='source log2 Illumina HT12-v3 intensity; not logged again',
                 cna_proxy='fraction of unique assayed gene calls with absolute CNA>=1; not FGA or HRD',main_result=results[0])
    return results,cohort,summary

def geo_metadata(path):
    fields={};platforms=set()
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):break
            parts=next(csv.reader([line],delimiter='\t'))
            if not parts:continue
            if parts[0]=='!Sample_geo_accession':fields['GSM']=parts[1:]
            elif parts[0]=='!Sample_platform_id':platforms.update(parts[1:])
            elif parts[0]=='!Sample_characteristics_ch1':
                key=parts[1].split(': ',1)[0];fields[key]=[s.split(': ',1)[1] if ': ' in s else '' for s in parts[1:]]
    if platforms!={'GPL96'}:raise ValueError('GSE25066 platform is not exclusively GPL96')
    return pd.DataFrame(fields)

def gse25066(root=PROJECT_ROOT):
    source=root.parent.parent/'dataset/GSE25066_series_matrix.txt.gz';dest=root/'outputs/v2/stage3'
    metadata=geo_metadata(source)
    mask=(metadata.er_status_ihc=='N')&(metadata.pr_status_ihc=='N')&(metadata.her2_status=='N')
    cohort=metadata.loc[mask].copy();cohort['pCR']=cohort.pathologic_response_pcr_rd.map({'pCR':1,'RD':0})
    cohort['DRFS_years']=pd.to_numeric(cohort.drfs_even_time_years,errors='coerce')
    cohort['DRFS_event']=pd.to_numeric(cohort.drfs_1_event_0_censored,errors='coerce')
    paths=[root/'data/external/gpl96/GPL96.annot.gz',root.parent.parent/'dataset/GPL96.annot.gz']
    annotation=next((p for p in paths if p.exists()),None)
    summary=dict(clinical_tnbc_n=len(cohort),pcr_observed_n=int(cohort.pCR.notna().sum()),pcr_events=int(cohort.pCR.sum()),platform='GPL96',
                 genes_available=[],genes_missing=[],gene_availability_status='unassessable without authoritative GPL96 annotation',annotation_missing=annotation is None)
    results=[]
    if annotation is None:
        for endpoint in ['pCR','DRFS']:
            results.append(dict(component='GSE25066',endpoint=endpoint,gene='reduced_platform_composite',clinical_N=len(cohort),N=0,status='unavailable',reason='GPL96 annotation absent; probe-to-gene mapping not inferred',beta=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,fdr=np.nan))
    else:
        with gzip.open(annotation,'rt') as f:
            line=f.readline()
            while line and not line.startswith('ID\t'):line=f.readline()
            if not line:raise ValueError('GPL96 annotation header missing')
            mapping=pd.read_csv(f,sep='\t',names=line.rstrip().split('\t'),comment='!')
        symbol_col='Gene symbol'
        if symbol_col not in mapping:raise ValueError('GPL96 gene symbol column missing')
        # Unambiguous probe annotations
        mapping=mapping[mapping[symbol_col].isin(CANDIDATE_GENES)].copy()
        lookup=dict(zip(mapping.ID.astype(str),mapping[symbol_col]))
        with gzip.open(source,'rt') as f:
            for line in f:
                if line.startswith('!series_matrix_table_begin'):break
            reader=csv.reader(f,delimiter='\t');head=next(reader);selected=[]
            for row in reader:
                if row[0].startswith('!series_matrix_table_end'):break
                if row[0] in lookup:selected.append([lookup[row[0]],*map(float,row[1:])])
        expr=pd.DataFrame(selected,columns=['Hugo_Symbol',*head[1:]])
        expr=gene_mean(expr)[cohort.GSM.tolist()]
        available=[g for g in CANDIDATE_GENES if g in expr.index]
        if not available:raise ValueError('No candidate genes mapped by GPL96 annotation')
        score,params,missing=signature(expr,available)
        cohort['candidate_mean_z']=cohort.GSM.map(score);cohort['candidate_perSD']=cohort.candidate_mean_z/cohort.candidate_mean_z.std(ddof=0)
        r,_,_=fit_model(cohort,outcome='pCR',label='GSE25066 reduced-platform pCR',binary=True);r.update(component='GSE25066',endpoint='pCR',gene='reduced_platform_composite');results.append(r)
        r=cox(cohort,'DRFS_years','DRFS_event',label='DRFS');r.update(component='GSE25066',gene='reduced_platform_composite');results.append(r)
        for r,q in zip(results,bh([r['p'] for r in results])):r['fdr']=q
        summary.update(genes_available=available,genes_missing=[g for g in CANDIDATE_GENES if g not in available],gene_availability_status='verified GPL96 probe mapping',signature_parameters=params)
    cohort.to_csv(dest/'gse25066_clinical_cohort.csv',index=False);pd.DataFrame(results).to_csv(dest/'gse25066_results.csv',index=False)
    return results,cohort,summary
