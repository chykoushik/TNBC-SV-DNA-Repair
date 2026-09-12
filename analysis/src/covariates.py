import numpy as np
import pandas as pd
from .config import CANDIDATE_GENES, COMPARISON_GENES
from .data_loading import load_clinical

# Nielsen 2010 proliferation
PROLIFERATION_GENES=('BIRC5','CCNB1','CDC20','NUF2','CEP55','NDC80','MKI67','PTTG1','RRM2','TYMS','UBE2C')
PROLIFERATION_SOURCE='https://unclineberger.org/peroulab/wp-content/uploads/sites/1008/2019/06/Nov-1-TAM-series-PAM50-Nielsen-et-al-CCR-2010.pdf'
assert not set(PROLIFERATION_GENES)&set(CANDIDATE_GENES)

def load_frozen(root):
    c=pd.read_csv(root/'outputs/v2/cohort/tcga_tnbc_patients.csv',keep_default_na=False)
    p=pd.read_csv(root/'outputs/v2/predictors/candidate_signature.csv')
    m=pd.read_csv(root/'outputs/v2/expression/expression_sample_manifest.csv',keep_default_na=False)
    m=m[m['selected'].astype(str).str.lower()=='true']
    if c.patient_id.duplicated().any() or p.patient_id.duplicated().any() or m.patient_id.duplicated().any():raise ValueError('Nonunique Stage 1 patient IDs')
    base=c[['patient_id','sample_id']].rename(columns={'sample_id':'clinical_sample_id'})
    base=base.merge(p,on='patient_id',how='left',validate='one_to_one')
    if not (base.clinical_sample_id==base.sample_id).all():raise ValueError('Frozen clinical/predictor sample conflict')
    base=base.rename(columns={'expression_sample_id':'rna_sample_id'})
    for r in m.to_dict('records'):
        match=base[base.patient_id==r['patient_id']]
        if len(match)!=1 or match.iloc[0].rna_sample_id!=r['expression_sample_id']:raise ValueError('RNA manifest conflict')
    x=pd.read_csv(root/'outputs/v2/expression/tcga_tnbc_expression.tsv.gz',sep='\t')
    samples=x.columns[3:].tolist()
    matched=base[base.rna_sample_id.notna()]
    if set(samples)!=set(matched.rna_sample_id):raise ValueError('Matrix/frozen predictor sample mismatch')
    # Explicit patient mapping
    sample_to_patient=dict(zip(matched.rna_sample_id,matched.patient_id))
    values=x[samples].astype(float);values.columns=[sample_to_patient[s] for s in samples]
    annotation=x.iloc[:,:3].copy()
    unique=~annotation.gene_symbol.duplicated(keep=False)
    expression=values.loc[unique].copy();expression.index=annotation.loc[unique,'gene_symbol']
    expression=expression.sort_index(axis=1)
    sd=expression.std(axis=1,ddof=0);z=expression.sub(expression.mean(axis=1),axis=0).div(sd.replace(0,np.nan),axis=0)
    candidate=z.loc[list(CANDIDATE_GENES)].T
    check=base.set_index('patient_id').loc[candidate.index,'candidate_mean_z']
    if not np.allclose(candidate.mean(axis=1),check,atol=1e-8,rtol=1e-8):raise ValueError('Reconstructed individual z-scores disagree with frozen composite')
    for g in CANDIDATE_GENES+COMPARISON_GENES:
        base[g+'_expr_z']=base.patient_id.map(z.loc[g] if g in z.index else {})
    missing=[g for g in PROLIFERATION_GENES if g not in z.index or z.loc[g].isna().any()]
    score=z.loc[list(PROLIFERATION_GENES)].mean(axis=0) if not missing else pd.Series(dtype=float)
    base['proliferation']=base.patient_id.map(score)
    return base,expression,z,annotation,missing

def add_clinical(base):
    clinical=pd.DataFrame(load_clinical())
    if clinical.sampleID.duplicated().any():raise ValueError('Clinical sample duplicates need adjudication')
    clinical=clinical.set_index('sampleID')
    for col,source in [('age','age_at_initial_pathologic_diagnosis'),('stage','pathologic_stage'),('PAM50','PAM50Call_RNAseq')]:
        base[col]=base.clinical_sample_id.map(clinical[source])
    base['age']=pd.to_numeric(base.age,errors='coerce')
    def stage(value):
        value=str(value).strip().upper().replace('STAGE ','')
        for s in ['IV','III','II','I']:
            if value in [s,s+'A',s+'B',s+'C']:return s
        return np.nan
    base['stage']=base.stage.map(stage)
    base['PAM50']=base.PAM50.replace({'':np.nan, 'NA':np.nan,'[Not Available]':np.nan})
    # Incomplete mutation ascertainment
    base['BRCA1_altered']=np.nan;base['BRCA2_altered']=np.nan
    return base
