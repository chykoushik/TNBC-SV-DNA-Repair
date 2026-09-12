import csv
import gzip
import json
import re
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from .config import CANDIDATE_GENES,COMPARISON_GENES
from .covariates import PROLIFERATION_GENES as PREVIOUS_PROLIFERATION
from .utils import sha256

STAGE4_PROLIFERATION=('MKI67','PCNA','CCNB1','CCNE1','MCM2','MCM4','MCM6','TOP2A','BIRC5','CDC20','CDC6')
PRIMARY=('HRDetect_probability','scarHRD_score')
SEED=20260909

def biological_key(identifier):
    s=re.findall(r'(?:^|\.)(S\d{6})(?=\.|$)',str(identifier))
    m=re.findall(r'(?:^|\.)(m\d*)(?=\.|$)',str(identifier))
    if len(s)!=1 or len(m)!=1:raise ValueError('Identifier must contain one S###### and one m component: '+str(identifier))
    return s[0],m[0]

def parse_soft(path):
    records=[];current=None;processing=set()
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            line=line.rstrip('\r\n')
            if line.startswith('^SAMPLE = '):
                if current:records.append(current)
                current={'sample':line.split(' = ',1)[1],'characteristics':[]}
            elif current and line.startswith('!Sample_title = '):current['title']=line.split(' = ',1)[1]
            elif current and line.startswith('!Sample_characteristics_ch1 = '):current['characteristics'].append(line.split(' = ',1)[1])
            elif line.startswith('!Sample_data_processing = '):processing.add(line.split(' = ',1)[1])
    if current:records.append(current)
    return records,sorted(processing)

def properties(record):
    out={}
    for c in record.get('characteristics',[]):
        if ': ' in c:
            key,value=c.split(': ',1)
            if key in out and out[key]!=value:raise ValueError('Conflicting GEO characteristic')
            out[key]=value
    return out

def map_samples(hrd,geo,expression_columns):
    pool=defaultdict(list);geo_titles=set();rows=[]
    if hrd.PD_ID.duplicated().any() or hrd.RBA.duplicated().any():raise ValueError('Duplicate WGS/RBA identity')
    seen_keys=set();hrd_s={biological_key(x)[0] for x in hrd.RBA}
    for r in geo:
        prop=properties(r);external=prop.get('scan-b external id','')
        if not external:continue
        title=r['title'].strip()
        if title in geo_titles:raise ValueError('Duplicate GEO title')
        geo_titles.add(title)
        try:key=biological_key(external)
        except ValueError:
            if set(re.findall(r'S\d{6}',external))&hrd_s:raise ValueError('Unresolvable GEO key for a cohort biological S identifier: '+external)
            continue  # Missing cohort identifiers
        pool[key].append((r,prop,external,title))
    for r in hrd.itertuples():
        key=biological_key(r.RBA)
        if key in seen_keys:raise ValueError('Multiple HRD patients share biological S/m key')
        seen_keys.add(key)
        candidates=pool[key]
        eligible=[x for x in candidates if not x[3].lower().endswith('repl')]
        ambiguous=len(eligible)>1
        if not candidates:
            rows.append(dict(PD_ID=r.PD_ID,HRD_RBA=r.RBA,S_identifier=key[0],m_component=key[1],GEO_accession='',GEO_external_ID='',F_expression_column='',mapping_status='unmatched',ambiguity_status=False,replicate_status=False,age=np.nan))
        for record,prop,external,title in candidates:
            repl=title.lower().endswith('repl')
            status='excluded_replicate' if repl else 'ambiguous' if ambiguous else 'mapped' if title in expression_columns else 'expression_column_missing'
            rows.append(dict(PD_ID=r.PD_ID,HRD_RBA=r.RBA,S_identifier=key[0],m_component=key[1],GEO_accession=record.get('geo',record['sample']),GEO_external_ID=external,F_expression_column=title,
                             mapping_status=status,ambiguity_status=ambiguous,replicate_status=repl,age=pd.to_numeric(prop.get('age at diagnosis'),errors='coerce')))
    result=pd.DataFrame(rows).sort_values(['PD_ID','replicate_status','F_expression_column']).reset_index(drop=True)
    mapped=result[result.mapping_status=='mapped']
    if mapped.F_expression_column.duplicated().any():raise ValueError('One expression column mapped to multiple WGS patients')
    return result

def require_unambiguous(mapping):
    if mapping.ambiguity_status.any():
        raise ValueError('Ambiguous mapping: aborting, see mapping audit')


def read_supplement(path):
    classes=pd.read_excel(path,sheet_name='HRD classifications');scores=pd.read_excel(path,sheet_name='HRD scores')
    if classes.PD_ID.duplicated().any() or scores.PD_ID.duplicated().any():raise ValueError('Duplicate PD_ID in supplement')
    if set(classes.PD_ID)!=set(scores.PD_ID):raise ValueError('Classification/score identity mismatch')
    data=classes.merge(scores,on='PD_ID',validate='one_to_one')
    for outcome in PRIMARY+('SBS3','SV3','delMH'):
        if outcome in data:data[outcome]=pd.to_numeric(data[outcome],errors='coerce')
    for field in ['HRDetect','scarHRD']:
        labels=data[field].dropna().unique()
        if not set(labels).issubset({'HRD','HRP'}):raise ValueError('Unexpected published class: '+field)
        data[field+'_binary']=data[field].map({'HRD':1.,'HRP':0.})
    if ((data.HRDetect_probability.dropna()<0)|(data.HRDetect_probability.dropna()>1)).any():raise ValueError('Probability out of bounds')
    return data

def stream_expression(path,columns,wanted=None,collect_moments=False):
    rows=[];moments=[];seen=set();kept=set(wanted or [])
    with gzip.open(path,'rt',encoding='utf-8') as f:header=next(csv.reader(f))
    first=header[0] or 'Unnamed: 0'
    for chunk in pd.read_csv(path,compression='gzip',usecols=[first,*columns],chunksize=256):
        names=chunk.iloc[:,0].astype(str)
        if names.duplicated().any() or set(names)&seen:raise ValueError('Duplicate gene symbol requires explicit annotation resolution')
        seen.update(names)
        values=chunk[columns].to_numpy(dtype=float)
        if not np.isfinite(values).all():raise ValueError('Nonfinite expression values')
        if collect_moments:
            means=values.mean(axis=1);variances=values.var(axis=1);detected=(values>=np.log2(.2)).mean(axis=1)
            moments.extend(zip(names,means,variances,detected))
        take=names.isin(kept).to_numpy()
        if take.any():rows.append(pd.DataFrame(values[take],index=names[take],columns=columns))
    expression=pd.concat(rows) if rows else pd.DataFrame(columns=columns)
    stats=pd.DataFrame(moments,columns=['gene','mean','variance','fraction_FPKM_ge_0_1']).set_index('gene') if collect_moments else None
    return expression,stats

def construct_predictors(expression,genes=CANDIDATE_GENES):
    if tuple(genes)!=tuple(CANDIDATE_GENES):raise ValueError('Candidate gene set is frozen')
    available=[g for g in genes if g in expression.index]
    missing=[g for g in genes if g not in expression.index]
    x=expression.loc[available];sd=x.std(axis=1,ddof=0);z=x.sub(x.mean(axis=1),axis=0).div(sd.replace(0,np.nan),axis=0)
    p=pd.DataFrame(index=expression.columns)
    p['candidate_mean_z']=z.mean(axis=0) if not missing and sd.gt(0).all() else np.nan
    p['candidate_median_z']=z.median(axis=0) if not missing and sd.gt(0).all() else np.nan
    for gene in available:p[gene+'_z']=z.loc[gene]
    for gene in ['HORMAD1','STAG3','REC8']:
        if gene+'_z' not in p:p[gene+'_z']=np.nan
    for label,genelist in [('proliferation_mean_z',STAGE4_PROLIFERATION),('previous_proliferation_mean_z',PREVIOUS_PROLIFERATION)]:
        present=[g for g in genelist if g in expression.index]
        e=expression.loc[present];esd=e.std(axis=1,ddof=0)
        p[label]=e.sub(e.mean(axis=1),axis=0).div(esd.replace(0,np.nan),axis=0).mean(axis=0) if present and esd.gt(0).all() else np.nan
    for g in COMPARISON_GENES:
        if g in expression.index:p[g+'_comparison_z']=(expression.loc[g]-expression.loc[g].mean())/expression.loc[g].std(ddof=0)
    params={'candidate_available':available,'candidate_missing':missing,'mean':x.mean(axis=1).to_dict(),'sd':sd.to_dict(),
            'proliferation_available':[g for g in STAGE4_PROLIFERATION if g in expression.index],
            'proliferation_missing':[g for g in STAGE4_PROLIFERATION if g not in expression.index]}
    return p,z,params

def null_neighborhoods(moments,genes=CANDIDATE_GENES):
    exclude=set(genes)|set(COMPARISON_GENES)|set(STAGE4_PROLIFERATION)|set(PREVIOUS_PROLIFERATION)|{'ESR1','PGR','ERBB2'}
    eligible=moments[(moments.variance>1e-8)&(moments.fraction_FPKM_ge_0_1>=.2)]
    pool=sorted(set(eligible.index)-exclude)
    if len(pool)<200:raise ValueError('Insufficient matched-gene background')
    names=pool+list(genes);coords=np.column_stack([rankdata(moments.loc[names,c])/len(names) for c in ['mean','variance']])
    neighborhoods={};distance={}
    for i,g in enumerate(genes):
        d=np.sqrt(((coords[:len(pool)]-coords[len(pool)+i])**2).sum(axis=1));indices=np.argsort(d,kind='stable')[:200]
        neighborhoods[g]=[pool[j] for j in indices];distance[g]={pool[j]:float(d[j]) for j in indices}
    return neighborhoods,distance,len(pool)

def sample_null_sets(neighbors,n=5000,seed=SEED):
    rng=np.random.default_rng(seed);sets=[]
    for _ in range(n):
        selected=[]
        for g in CANDIDATE_GENES:
            options=[x for x in neighbors[g] if x not in selected];selected.append(str(rng.choice(options)))
        sets.append(selected)
    return sets

def prepare(root):
    dest=root/'outputs/v2/stage4';dest.mkdir(parents=True,exist_ok=True)
    paths=[root/'data'/n for n in ['13058_2026_2325_MOESM2_ESM.xlsx','13058_2026_2325_MOESM3_ESM.xlsx','GSE96058_gene_expression.csv.gz','GSE96058_family.soft.gz']]
    stamp={str(p):{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in paths}
    cache=dest/'preparation_manifest.json'
    if cache.exists() and json.loads(cache.read_text())['input_stat']==stamp:
        print('Using validated cache',flush=True);return
    geo,processing=parse_soft(paths[3]);cached=root/'scanb_geo_samples.json'
    if cached.exists():
        previous=json.loads(cached.read_text())
        key=lambda r:(r['sample'],r['title'],properties(r).get('scan-b external id'))
        if sorted(map(key,previous))!=sorted(map(key,geo)):raise ValueError('Existing GEO JSON disagrees with authoritative local SOFT identity')
    else:cached.write_text(json.dumps(geo,indent=2))
    if not any('0.1 FPKM' in s and 'log2' in s for s in processing):raise ValueError('Expression transformation not confirmed in local GEO metadata')
    hrd=read_supplement(paths[0])
    with gzip.open(paths[2],'rt') as f:columns=next(csv.reader(f))[1:]
    mapping=map_samples(hrd,geo,columns);mapping.to_csv(dest/'scanb_sample_mapping.csv',index=False)
    matched=mapping[mapping.mapping_status=='mapped']
    audit={'supplement_n':len(hrd),'uniquely_mapped_n':len(matched),'unmatched_n':len(hrd)-matched.PD_ID.nunique(),
           'ambiguous_n':mapping.loc[mapping.ambiguity_status,'PD_ID'].nunique(),'replicate_rows_excluded':int(mapping.replicate_status.sum()),
           'matching_key':'biological S###### plus exact m/m2/m3/m4 component; k/k2 ignored',
           'source_processing':processing,'json_validated_against_soft':True}
    (dest/'scanb_mapping_audit.json').write_text(json.dumps(audit,indent=2))
    require_unambiguous(mapping)
    joined=hrd.merge(matched[['PD_ID','F_expression_column','GEO_accession','GEO_external_ID','age']],on='PD_ID',how='inner',validate='one_to_one').sort_values('PD_ID').reset_index(drop=True)
    cols=joined.F_expression_column.tolist();wanted=set(CANDIDATE_GENES)|set(COMPARISON_GENES)|set(STAGE4_PROLIFERATION)|set(PREVIOUS_PROLIFERATION)
    print(f'Streaming expression: {len(wanted)}/{len(cols)}',flush=True)
    expression,moments=stream_expression(paths[2],cols,wanted,True)
    predictors,z,params=construct_predictors(expression)
    joined=joined.merge(predictors.rename_axis('F_expression_column').reset_index(),on='F_expression_column',validate='one_to_one')
    joined.to_csv(dest/'scanb_analysis_cohort.csv',index=False)
    predictors.rename_axis('F_expression_column').reset_index().merge(joined[['PD_ID','RBA','F_expression_column']],on='F_expression_column',validate='one_to_one').to_csv(dest/'scanb_expression_predictors.csv',index=False)
    expression.to_csv(dest/'scanb_required_gene_expression.tsv.gz',sep='\t',compression='gzip',index_label='gene')
    moments.to_csv(dest/'scanb_gene_moments.csv')
    (dest/'scanb_standardization.json').write_text(json.dumps(params,indent=2))
    if params['candidate_missing']:raise ValueError('Missing frozen genes; complete primary signature cannot be constructed: '+','.join(params['candidate_missing']))
    neighbors,distances,pool_n=null_neighborhoods(moments)
    union=sorted(set(x for v in neighbors.values() for x in v))
    print(f'Streaming genes: {len(union)}',flush=True)
    null_expression,_=stream_expression(paths[2],cols,union,False)
    null_expression.to_csv(dest/'scanb_null_neighborhood_expression.tsv.gz',sep='\t',compression='gzip',index_label='gene')
    (dest/'scanb_null_neighborhoods.json').write_text(json.dumps({'neighbors':neighbors,'distances':distances,'eligible_pool_n':pool_n},indent=2))
    manifest={'input_stat':stamp,'input_sha256':{str(p):sha256(p) for p in paths},'gene_rows':len(moments),'mapped_expression_n':len(cols),
              'null_loaded_gene_n':len(union),'eligible_pool_n':pool_n,'expression_transform':'already log2(FPKM+0.1); identity used',
              'supplement3_use':'Inspected; gene lists/DEG results, no additional patient-level covariates; excluded from predictor selection'}
    cache.write_text(json.dumps(manifest,indent=2))
    print('Expression preparation complete',flush=True)

if __name__=='__main__':
    from .config import PROJECT_ROOT
    prepare(PROJECT_ROOT)
