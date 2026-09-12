from pathlib import Path
import numpy as np
import pandas as pd
from .utils import tcga_ids

SOURCES = {
 'tcga_hrd': ('TCGA.HRD_withSampleID.txt','sampleID',{'HRD_total':'HRD','HRD_LOH':'hrd-loh','LST':'lst1','TAI':'ai1'}),
 'tcga_aneuploidy': ('ABSOLUTE_scores.tsv',None,{'aneuploidy':'AS','LOH_fraction':'LOH_frac_altered'}),
 'tcga_fga': ('seg_based_scores.tsv','Sample',{'FGA':'frac_altered'}),
 'tcga_absolute': ('TCGA_mastercalls.abs_tables_JSedit.fixed.txt','sample',{'ploidy':'ploidy','WGD':'Genome doublings','purity':'purity'})}
OUTCOMES = ('HRD_total','HRD_LOH','LST','TAI','aneuploidy','FGA','LOH_fraction','ploidy','WGD')

def match_record(frame, clinical_id, rna_id):
    cid=tcga_ids(clinical_id)
    sub=frame[frame['_patient']==cid['patient_id']]
    if len(sub)==0:return None,'no_matching_outcome'
    exact=sub[sub['_source_id']==rna_id]
    if len(exact):sub=exact
    else:sub=sub[sub['_sample_key']==cid['sample_key']]
    if len(sub)==0:return None,'no_matching_primary_sample'
    cols=[c for c in frame.columns if not c.startswith('_')]
    if len(sub)>1 and any(sub[c].nunique(dropna=False)>1 for c in cols):return None,'conflicting_genomic_records'
    row=sub.sort_values('_source_id').iloc[0].copy()
    row['_source_id']=';'.join(sorted(sub['_source_id'].unique()))
    return row,'matched' if len(sub)==1 else 'identical_duplicate_values_collapsed'

def load_outcomes(root, cohort):
    out=pd.DataFrame(index=cohort.index)
    availability=[];conflicts=[]
    for source,(filename,idcol,mapping) in SOURCES.items():
        known=[root/'data/external'/source/filename,root.parent.parent/'dataset'/filename]
        paths=[p for p in known if p.is_file()]
        for name in mapping:out[name]=np.nan;out[name+'_sample_id']=''
        if len(paths)!=1:
            availability.append(dict(resource=source,filename=filename,status='missing' if not paths else 'ambiguous_source_paths',path=''))
            continue
        path=paths[0];frame=pd.read_csv(path,sep='\t',dtype=str)
        idcol=idcol or frame.columns[0]
        if not {idcol,*mapping.values()}.issubset(frame.columns):raise ValueError('Outcome schema mismatch: '+str(path))
        data=frame[list(mapping.values())].rename(columns={v:k for k,v in mapping.items()}).apply(pd.to_numeric,errors='coerce')
        data['_source_id']=frame[idcol].str.strip().str.upper(); patients=[];keys=[]
        for value in data['_source_id']:
            try:
                ids=tcga_ids(value);patients.append(ids['patient_id']);keys.append(ids['sample_key'] if ids['sample_type_code']=='01' else '')
            except (ValueError,AttributeError):patients.append('');keys.append('')
        data['_patient']=patients;data['_sample_key']=keys
        if 'WGD' in data:
            # Curated genome doublings
            data['WGD']=data['WGD'].map(lambda x:float(x>=1) if pd.notna(x) and x>=0 and x==int(x) else np.nan)
        for i,row in cohort.iterrows():
            selected,status=match_record(data,row['clinical_sample_id'],row['rna_sample_id'])
            if selected is not None:
                for name in mapping:out.loc[i,name]=selected[name];out.loc[i,name+'_sample_id']=selected['_source_id']
            if status!='matched':conflicts.append(dict(patient_id=row['patient_id'],resource=source,status=status))
        availability.append(dict(resource=source,filename=filename,status='loaded',path=str(path)))
    return out,availability,conflicts
