import sys, json, os, subprocess, platform, io, unittest, csv
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.config import CANDIDATE_GENES, COMPARISON_GENES
from src.covariates import PROLIFERATION_GENES, load_frozen
from src.stage4_scanb import STAGE4_PROLIFERATION
from src.utils import sha256
from src.stage2 import clean
FINAL=ROOT/'outputs/v2/final'
FIGURES=['Figure1_study_design','Figure2_TCGA_primary_HRD','Figure3_SCANB_replication','Figure4_specificity_and_robustness','Figure5_cross_cohort_evidence']
TABLES=['Table1_cohorts','Table2_primary_results','Table3_robustness','Table4_external_validation','TableS1_gene_definitions','TableS2_all_models','TableS3_negative_results','TableS4_sample_mapping','TableS5_reproducibility']
DOCS=['V2_FINAL_METHODS.md','V2_LIMITATIONS.md','V2_FINAL_EVIDENCE_REVIEW.md']
LIMITATIONS=[
 'Proliferation confounding: SCAN-B unadjusted associations attenuate substantially and are no longer significant after proliferation adjustment; an independent biomarker is not established.',
 'No experimental perturbation, rescue or direct repair-function experiment establishes causal action of the candidate program.',
 'SCAN-B lacks patient-level BRCA1/2 alteration and purity covariates in the supplied tables; TCGA BRCA denominators are also unresolved.',
 'GSE25066 pCR and DRFS results are negative and test only four platform-covered genes, not the full nine-gene signature.',
 'DepMap has only 25 jointly measurable curated TNBC models and no genome-wide or DDR dependencies surviving genome-wide FDR; a negative screen does not establish absence of any dependency.',
 'Single-cell malignant-specificity test is negative in eight paired TNBC donor/sample units; donor identity and annotation limitations remain.',
 'METABRIC measures gene-weighted altered-gene fraction, not direct HRD scars or length-weighted fraction genome altered; its adjusted association is nonsignificant.',
 'SCAN-B matched-random scarHRD result is borderline and above 0.05 (P=0.05598880223955209). Matching distributions differ between TCGA and SCAN-B backgrounds.',
 'HRDetect and scarHRD share genomic information; SBS3, SV3 and delMH overlap HRDetect features and are not independent replications.',
 'Within-cohort expression standardization is cross-platform restandardization, not a transferable individual-patient assay or validated cutoff.',
 'SCAN-B nonlinear relationships and nonnormal residuals limit OLS interpretation; fractional-logit sensitivity also attenuates after proliferation adjustment.',
 'TCGA has 114 clinical patients, 113 RNA matches and 103 primary HRD cases; SCAN-B loses 10 of 235 to exact identity matching. Missingness may select patients.',
 'TCGA detailed aliquot barcodes are unavailable for many records; sample plus GDC file UUID is explicit, not inferred. Technical read-count ranking is not comprehensive pathology QC.',
 'GSE25066 DRFS Schoenfeld-residual/log-time diagnostic P=0.04965547910807285 flags possible proportional-hazards departure; this is a limited diagnostic.',
 'Observational multicohort association establishes neither clinical utility, treatment interaction, therapeutic selectivity nor real-time HR repair capacity.'
]

def readj(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def readc(p):return pd.read_csv(ROOT/p)
def savej(path,obj):path.write_text(json.dumps(clean(obj),indent=2,allow_nan=False),encoding='utf-8')
def rel(p):return str(Path(p).relative_to(ROOT)).replace('\\','/') if Path(p).is_relative_to(ROOT) else str(p)
def value(r,*keys,default=np.nan):
    for k in keys:
        x=r.get(k)
        if x is not None and not (isinstance(x,float) and np.isnan(x)) and x!='':return x
    return default
def fmt(x):
    if pd.isna(x):return 'NA'
    return f'{x:.3g}' if isinstance(x,(int,float,np.number)) else str(x)

def manuscript_files():
    found=[]
    for base,dirs,files in os.walk(ROOT.parent.parent):
        dirs[:]=[d for d in dirs if d not in {'dataset','data','results','outputs','.git','node_modules','__pycache__','.venv'}]
        for name in files:
            p=Path(base)/name
            if p.suffix.lower() in {'.tex','.bib','.docx','.doc','.odt'} or ('manuscript' in name.lower() and p.suffix.lower() in {'.md','.txt'}):found.append(p)
    return sorted(found)

def freeze_and_tests():
    protected=[]
    for folder in ['cohort','expression','predictors','stage2','stage3','stage4']:
        protected += [p for p in (ROOT/'outputs/v2'/folder).rglob('*') if p.is_file()]
    protected += list((ROOT/'src').glob('*.py'))+list((ROOT/'docs').glob('V2_METHODS_STAGE*.md'))+manuscript_files()
    baseline={str(p):sha256(p) for p in protected}
    baseline_path=FINAL/'freeze_baseline.json'
    if baseline_path.exists():
        old=readj(rel(baseline_path))
        changed=[p for p,h in old['protected'].items() if not Path(p).exists() or sha256(p)!=h]
        if changed:raise RuntimeError('Frozen baseline mismatch: '+str(changed))
    else:savej(baseline_path,{'created_utc':datetime.now(timezone.utc).isoformat(),'protected':baseline,'manuscripts':{str(p):sha256(p) for p in manuscript_files()}})
    historical=[]
    def check(path,expected,origin):
        p=Path(path);p=p if p.is_absolute() else ROOT/p
        actual=sha256(p) if p.exists() else None
        historical.append(dict(path=str(p),expected_sha256=expected,actual_sha256=actual,status='passed' if actual==expected else 'FAILED',reference=origin))
    for doc,key in [('stage2/stage2_summary.json','source_hashes'),('stage3/stage3_summary.json','stage1_stage2_source_hashes'),('stage4/stage4_summary.json','frozen_previous_input_hashes'),('stage4/preparation_manifest.json','input_sha256')]:
        for p,h in readj('outputs/v2/'+doc)[key].items():check(p,h,doc)
    for r in readj('outputs/v2/stage4/stage4_artifact_manifest.json'):check(r['path'],r['sha256'],'stage4_artifact_manifest.json')
    c=readj('outputs/v2/cohort/cohort_summary.json');check(c['clinical_file'],c['clinical_sha256'],'cohort_summary.json')
    for r in readc('outputs/v2/expression/expression_sample_manifest.csv').to_dict('records'):check(r['source_path'],r['sha256'],'expression_sample_manifest.csv')
    sc=readj('outputs/v2/stage3/single_cell_source.json');check(next(iter(sc['source_stat'])),sc['archive_sha256'],'single_cell_source.json')
    pd.DataFrame(historical).to_csv(FINAL/'historical_hash_validation.csv',index=False)
    if any(r['status']!='passed' for r in historical):raise RuntimeError('Historical frozen hash mismatch; see historical_hash_validation.csv')
    load_frozen(ROOT)  # Identity consistency
    log=io.StringIO();suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    result=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
    (FINAL/'tests.log').write_text(log.getvalue(),encoding='utf-8')
    failures=[dict(test=t.id(),traceback=s,classification='infrastructure/acquisition' if 'test_acquisition.' in t.id() else 'scientific/package') for t,s in result.failures+result.errors]
    report=dict(run=result.testsRun,passed=result.testsRun-len(failures)-len(result.skipped),failed=len(failures),skipped=len(result.skipped),failures=failures,scientific_tests_passed=all(f['classification']=='infrastructure/acquisition' for f in failures),existing_read_only_load_frozen='passed',historical_hash_checks=len(historical))
    savej(FINAL/'test_results.json',report)
    if not report['scientific_tests_passed']:raise RuntimeError('Scientific or package test failure')
    return report

def evidence():
    rows=[]
    def add(r,source,i,stage,cohort,kind='secondary/exploratory'):
        endpoint=value(r,'endpoint','outcome',default='HRD_total');model=value(r,'model',default='Cox PH' if endpoint in ['OS','RFS','DRFS'] else 'paired Wilcoxon' if cohort=='GSE176078' else 'summary')
        predictor=value(r,'predictor','gene',default='candidate_mean_z')
        status=value(r,'status',default='estimated');beta=value(r,'beta');effect=value(r,'effect','beta');lo=value(r,'ci_low');hi=value(r,'ci_high');measure='OLS beta per signature unit';ci='conventional 95% CI'
        if pd.notna(r.get('odds_ratio',np.nan)):effect=r['odds_ratio'];lo=r['or_ci_low'];hi=r['or_ci_high'];measure='odds ratio per predictor SD'
        elif pd.notna(r.get('hazard_ratio',np.nan)):effect=r['hazard_ratio'];lo=r['hr_ci_low'];hi=r['hr_ci_high'];measure='hazard ratio per predictor SD'
        elif cohort=='GSE176078':measure='paired mean difference';ci='paired-donor bootstrap 95% CI'
        elif 'Fractional' in str(model):measure='logit conditional mean coefficient';ci='HC0 sandwich 95% CI'
        if 'bootstrap' in source:measure='OLS beta; bootstrap uncertainty';effect=value(r,'observed_beta','median_beta');ci='percentile bootstrap 95% CI'
        if 'permutation' in source or 'null_signature_summary' in source or 'random_gene_null_summary' in source:measure='observed absolute t';effect=value(r,'observed_abs_t');beta=value(r,'observed_beta');lo=hi=np.nan
        p=value(r,'p','empirical_p');robust=value(r,'hc3_p');q=value(r,'fdr');significance=q if pd.notna(q) else robust if pd.notna(robust) else p
        supported=status=='estimated' and pd.notna(significance) and significance<.05
        limitation='Observational; endpoint-specific complete cases; no causal or clinical utility inference.'
        if cohort=='METABRIC':limitation+=' CNA is a gene-weighted DNA proxy, not direct HRD; survival endpoints are secondary.'
        if cohort=='GSE25066':limitation+=' Reduced four-gene platform score only; full nine-gene validation unavailable.';predictor='four-gene mean z (STAG3, REC8, SYCP2, MSH4)'
        if cohort=='GSE176078':limitation+=' Donor/sample pairs are units; no independent HRD annotation.'
        if cohort=='SCAN-B':limitation+=' Proliferation attenuation; BRCA/purity unavailable; correlated genomic endpoints.'
        if endpoint in ['proliferation','OS','RFS','DRFS','pCR','malignant_vs_nonmalignant_combined'] or cohort=='GSE176078':claim='CONTEXT_ONLY'
        elif kind=='primary' and supported:claim='QUALIFIED_SUPPORT'
        elif not supported:claim='NO_SUPPORT_AT_THRESHOLD' if status=='estimated' else 'UNAVAILABLE'
        else:claim='EXPLORATORY_SUPPORT' if pd.notna(beta) and beta>0 else 'CONTEXT_ONLY'
        if cohort in ['TCGA','SCAN-B'] and predictor not in ['candidate_mean_z','candidate_median_z','leave_one_out_score']:claim='CONTEXT_ONLY'
        interpretation='Association detected; scope limited to this model.' if supported else 'No statistical support at the recorded threshold; not proof of no effect.' if status=='estimated' else 'Not estimable: '+str(value(r,'reason',default='unavailable'))
        row=dict(stage=stage,cohort=cohort,population=value(r,'population',default='clinically defined TNBC'),N=value(r,'N',default=103 if cohort=='TCGA' else np.nan),endpoint=endpoint,predictor=predictor,analysis_type=model,beta=beta,effect_measure=measure,effect=effect,ci_low=lo,ci_high=hi,ci_type=ci,p=p,robust_p=robust,adjusted_p=value(r,'primary_holm_p'),fdr=q,r_squared=value(r,'r_squared'),spearman_rho=value(r,'spearman_rho'),events=value(r,'HRD_events','events'),status=status,interpretation=interpretation,primary_secondary=kind,supports_main_claim=claim,limitation=limitation,source_output=source,source_row=i,robust_ci_low=value(r,'hc3_ci_low'),robust_ci_high=value(r,'hc3_ci_high'),robust_adjusted_p=value(r,'primary_hc3_holm_p'),omitted_gene=value(r,'omitted_gene',default=''),AUC=value(r,'AUC_apparent'),seed=value(r,'seed'),iterations=value(r,'iterations'),percentile=value(r,'percentile'))
        if cohort=='GSE25066' and endpoint=='pCR':row['events']=readj('outputs/v2/stage3/stage3_summary.json')['gse25066']['pcr_events']
        rows.append(row)
    specs={
      'stage2':(['primary_hrd_model','adjusted_models','sensitivity_analysis','single_gene_results','secondary_outcomes','canonical_comparison','leave_one_gene_out','proliferation_analysis'],'TCGA'),
      'stage4':(['scanb_primary_results','scanb_adjusted_results','scanb_secondary_results','scanb_gene_sensitivity','scanb_leave_one_gene_out','scanb_canonical_results','scanb_permutation_results','scanb_bootstrap_results','scanb_null_signature_summary'],'SCAN-B')}
    for folder,(files,cohort) in specs.items():
        for name in files:
            source=f'outputs/v2/{folder}/{name}.csv'
            for i,r in enumerate(readc(source).to_dict('records')):add(r,source,i,int(folder[-1]),cohort,'primary' if name in ['primary_hrd_model','scanb_primary_results'] else 'sensitivity' if any(x in name for x in ['adjust','sensitivity','leave_one','permutation','bootstrap','null','proliferation']) else 'secondary/exploratory')
    for name,cohort in [('metabric_results','METABRIC'),('gse25066_results','GSE25066'),('single_cell_results','GSE176078')]:
        source=f'outputs/v2/stage3/{name}.csv'
        for i,r in enumerate(readc(source).to_dict('records')):add(r,source,i,3,cohort,'external/contextual')
    for name in ['permutation_summary','bootstrap_summary','random_gene_null_summary']:
        source='outputs/v2/stage2/'+name+'.json';add(readj(source),source,0,2,'TCGA','sensitivity')
    dep=readc('outputs/v2/stage3/depmap_results.csv');ddr=dep[dep.annotation.fillna('').str.contains('DDR',case=False)]
    for label,d in [('genome-wide',dep),('DDR annotation',ddr)]:
        r=dict(endpoint='CRISPR '+label,N=25,effect=int((d.fdr_adjusted<.05).sum()),status='estimated',model='screen summary; proliferation-adjusted genome-wide FDR',predictor='candidate_mean_z',fdr=d.fdr_adjusted.min())
        add(r,'outputs/v2/stage3/depmap_results.csv',-1,3,'DepMap','external/contextual')
        rows[-1].update(effect_measure='number of dependencies with genome-wide FDR <0.05',interpretation=f"{int(r['effect'])} discoveries among {len(d)} annotated/tested dependencies; fdr column is minimum genome-wide adjusted q, not a global test.",supports_main_claim='NO_SUPPORT_AT_THRESHOLD',limitation='25 curated TNBC models; low power; annotations do not define a new restricted FDR family.',population='curated TNBC cell-line models',tested_dependencies=len(d))
    out=pd.DataFrame(rows);out.insert(0,'analysis_id',['E'+str(i+1).zfill(4) for i in range(len(out))]);return out,dep

def claims(e):
    entries=[
      ('A','Candidate program is associated with HRD/genomic instability in TNBC','TCGA HRD and SCAN-B WGS endpoints; METABRIC CNA proxy','Proliferation attenuation; TCGA aneuploidy/ploidy/WGD not supported','SUPPORTED_WITH_QUALIFICATION','Retain as observational association with specific genomic endpoints.'),
      ('B','Unadjusted association independently replicates in SCAN-B','Both WGS endpoints positive with HC3 and Holm support in 225 patients','Independent cohort does not mean confounder-independent or statistically independent endpoints','SUPPORTED','Retain explicitly as independent-cohort unadjusted replication.'),
      ('C','Association is independent of proliferation','TCGA adjusted HRD retains support','SCAN-B both adjusted tests nonsignificant; METABRIC adjusted test nonsignificant','NOT_SUPPORTED','Remove general independence claim; show attenuation prominently.'),
      ('D','Candidate program predicts chemotherapy response','No supporting GSE25066 pCR evidence','pCR P=0.943483; only four platform genes, no treatment interaction','NOT_SUPPORTED','Remove response prediction claim; report reduced-platform negative evidence.'),
      ('E','Candidate program predicts survival','No composite survival support','GSE25066 DRFS and METABRIC OS/RFS nonsignificant','NOT_SUPPORTED','Remove prognostic claim; retain negative results and precision.'),
      ('F','Candidate program creates a CRISPR dependency phenotype','No genome-wide FDR-supported dependency','Zero adjusted genome-wide and DDR discoveries; observational expression-dependency regressions','NOT_SUPPORTED','Remove selective genetic vulnerability claim.'),
      ('G','Candidate program is malignant-cell specific','Program measurable in pseudobulks','Eight TNBC paired units: P=0.546875; mixed cellular expression','NOT_SUPPORTED','Remove malignant-specificity claim; retain contextual result.'),
      ('H','Candidate program is an independent clinical HRD biomarker','Unadjusted genomic associations','Proliferation confounding; no clinical calibration, utility or prospective validation','NOT_SUPPORTED','Remove clinical biomarker and therapeutic selection language.'),
      ('I','Program may represent a transcriptional correlate of HRD-associated tumor state','TCGA discovery and SCAN-B direct genomic replication','Observational association; shared proliferation state plausible','SUPPORTED','Use as central conservative framing.'),
      ('J','HORMAD1 alone explains the entire composite signal','HORMAD1 is individually strong','TCGA and SCAN-B leave-HORMAD1-out associations remain significant; residual correlation cannot partition causal contribution','NOT_SUPPORTED','Reject entire-signal assertion; report individual strength and leave-one-out stability.'),
      ('K','Meiotic/cohesin transcription mechanistically causes HRD','No perturbation evidence','Cross-sectional expression and genomic scars cannot establish direction or mechanism','NOT_SUPPORTED','Remove causal verbs; reserve as an experimentally testable hypothesis.')]
    ledger=pd.DataFrame(entries,columns=['claim_id','proposed_claim','evidence_supporting','evidence_against_or_limiting','allowed_strength','manuscript_action'])
    refs={'A':'TCGA;SCAN-B;METABRIC','B':'SCAN-B','C':'TCGA;SCAN-B;METABRIC','D':'GSE25066','E':'GSE25066;METABRIC','F':'DepMap','G':'GSE176078','H':'TCGA;SCAN-B;GSE25066','I':'TCGA;SCAN-B','J':'TCGA;SCAN-B','K':'TCGA;SCAN-B'}
    ledger['evidence_analysis_ids']=ledger.claim_id.map(lambda k:';'.join(e.loc[e.cohort.isin(refs[k].split(';')),'analysis_id']))
    return ledger

def cohorts():
    return pd.DataFrame([
      ['TCGA-BRCA',114,113,103,'patients','HRD_total','outputs/v2/stage2/stage2_summary.json','clinical_n;rna_n;HRD_overlap_n'],
      ['SCAN-B',235,225,225,'patients','HRDetect; scarHRD','outputs/v2/stage4/stage4_summary.json','mapping.supplement_n;mapping.uniquely_mapped_n'],
      ['METABRIC',258,258,258,'patients','altered gene fraction; OS; RFS','outputs/v2/stage3/stage3_summary.json','metabric.clinical_tnbc_n'],
      ['GSE25066',178,178,112,'patients','pCR (37 events); DRFS N=114 (37 events)','outputs/v2/stage3/gse25066_results.csv','rows 0/1; stage3_summary.gse25066.clinical_tnbc_n'],
      ['DepMap',34,25,25,'curated TNBC models','CRISPR gene effect','outputs/v2/stage3/depmap_summary.json','tnbc_model_n;analysis_model_n'],
      ['GSE176078',10,10,8,'TNBC donor/sample units','malignant vs pooled nonmalignant','outputs/v2/stage3/single_cell_summary.json','tnbc_donor_sample_n;primary_result.N']
    ],columns=['cohort','eligible_N','expression_or_context_N','primary_analysis_N','unit','endpoint','source_output','source_field'])

def table(name,df,latex=True):
    df.to_csv(FINAL/'tables'/f'{name}.csv',index=False)
    if latex:
        # Compact table views
        cols=[c for c in ['analysis_id','cohort','endpoint','N','effect','ci_low','ci_high','p','robust_p','fdr','status'] if c in df]
        view=df[cols] if len(cols)>2 else df
        (FINAL/'tables'/f'{name}.tex').write_text(view.to_latex(index=False,escape=True,longtable=True,float_format=lambda x:f'{x:.5g}'),encoding='utf-8')

def make_tables(e,dep,c,test):
    table('Table1_cohorts',c)
    table('Table2_primary_results',e[e.primary_secondary=='primary'])
    table('Table3_robustness',e[e.primary_secondary=='sensitivity'])
    external=e[(e.stage==3)&((e.predictor.isin(['candidate_mean_z','composite']))|e.cohort.isin(['GSE25066','DepMap']))]
    table('Table4_external_validation',external)
    genes=[]
    for role,items in [('candidate',CANDIDATE_GENES),('canonical comparison',COMPARISON_GENES),('TCGA/METABRIC proliferation',PROLIFERATION_GENES),('SCAN-B/DepMap proliferation',STAGE4_PROLIFERATION)]:
        for gene in items:genes.append(dict(role=role,gene=gene,source_output='src/config.py' if role in ['candidate','canonical comparison'] else 'src/covariates.py' if role.startswith('TCGA') else 'src/stage4_scanb.py'))
    table('TableS1_gene_definitions',pd.DataFrame(genes))
    allmodels=[e]
    for adjusted in [False,True]:
        suf='_adjusted' if adjusted else '';d=pd.DataFrame({'cohort':'DepMap','endpoint':dep.dependency,'N':dep['N'+suf],'effect':dep['beta'+suf],'ci_low':dep['ci_low'+suf],'ci_high':dep['ci_high'+suf],'p':dep['p'+suf],'fdr':dep['fdr'+suf],'analysis_type':'proliferation-adjusted OLS' if adjusted else 'unadjusted OLS','source_output':'outputs/v2/stage3/depmap_results.csv','source_row':np.arange(len(dep))})
        allmodels.append(d)
    table('TableS2_all_models',pd.concat(allmodels,ignore_index=True),latex=False)
    negative=e[(e.status=='estimated')&((e.supports_main_claim=='NO_SUPPORT_AT_THRESHOLD')|((e.cohort.isin(['GSE25066','GSE176078']))&(e.fdr>=.05))|((e.endpoint.isin(['OS','RFS']))&(e.p>=.05)))]
    table('TableS3_negative_results',negative)
    mapping=[]
    for cohort,p in [('TCGA','outputs/v2/expression/expression_sample_manifest.csv'),('SCAN-B','outputs/v2/stage4/scanb_sample_mapping.csv')]:
        d=readc(p);d.insert(0,'cohort',cohort);d['source_output']=p;mapping.append(d)
    table('TableS4_sample_mapping',pd.concat(mapping,ignore_index=True),latex=False)
    table('TableS5_reproducibility',pd.DataFrame([{'check':k,'result':json.dumps(v),'source_output':'outputs/v2/final/test_results.json'} for k,v in test.items() if k!='failures']))

def figures(e):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    from scipy.stats import t
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':400})
    blue='#286C8E';orange='#B9582A';gray='#59616B'
    def save(fig,name):
        for ext in ['pdf','png']:fig.savefig(FINAL/'figures'/f'{name}.{ext}',bbox_inches='tight')
        plt.close(fig)
    def scatter(ax,data,outcome,r,label):
        d=data[['candidate_mean_z',outcome]].dropna();x=d.candidate_mean_z.to_numpy();y=d[outcome].to_numpy();b=r['beta'];a=y.mean()-b*x.mean();X=np.column_stack([np.ones(len(x)),x]);inv=np.linalg.inv(X.T@X)
        residual=y-a-b*x;h=(X@inv*X).sum(axis=1);cov=inv@(X.T@((residual/(1-h))[:,None]**2*X))@inv
        grid=np.linspace(x.min(),x.max(),160);G=np.column_stack([np.ones(len(grid)),grid]);fit=a+b*grid;half=t.ppf(.975,len(x)-2)*np.sqrt((G@cov*G).sum(axis=1))
        ax.scatter(x,y,s=20,alpha=.65,color=blue,edgecolors='none');ax.plot(grid,fit,color=orange,lw=1.7);ax.fill_between(grid,fit-half,fit+half,color=orange,alpha=.15)
        ax.set(xlabel='Candidate mean z-score',ylabel=label,title=f"N={len(d)} | HC3 P={r['hc3_p']:.3g}")
        ax.margins(.07);ax.grid(alpha=.15)
    fig,ax=plt.subplots(figsize=(12,7));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    boxes=[(.03,.64,.40,.25,'TCGA discovery\n114 clinical / 113 RNA / 103 HRD\nFrozen nine-gene predictor; genomic scars'),(.56,.64,.40,.25,'SCAN-B direct WGS replication\n235 supplied / 225 mapped\nHRDetect and scarHRD'),(.03,.30,.28,.23,'METABRIC orthogonal DNA\n258 clinical TNBC\nGene-weighted CNA proxy'),(.36,.30,.28,.23,'GSE25066 clinical outcomes\n178 TNBC; four platform genes\npCR N=112 / DRFS N=114'),(.69,.30,.28,.23,'Functional and cellular context\nDepMap: 25 TNBC models\nSingle-cell: 8 TNBC pairs')]
    for x,y,w,h,text in boxes:ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.01',facecolor='#EDF3F5',edgecolor=blue));ax.text(x+w/2,y+h/2,text,ha='center',va='center',linespacing=1.7)
    ax.annotate('',(.55,.765),(.44,.765),arrowprops=dict(arrowstyle='->',color=gray,lw=1.5))
    for end in [.17,.50,.83]:ax.annotate('',(end,.54),(.75,.63),arrowprops=dict(arrowstyle='->',color=gray))
    ax.text(.5,.13,'Parallel evidence streams; arrows do not imply patient overlap.\nRNA defines the predictor; independent DNA defines genomic outcomes.\nPositive and negative evidence retained; no outcome-driven gene selection.',ha='center',va='center',linespacing=1.6)
    ax.set_title('Study design and frozen analysis populations',fontsize=16,pad=15);save(fig,FIGURES[0])
    data=readc('outputs/v2/stage2/analysis_cohort.csv');rs=readc('outputs/v2/stage2/primary_hrd_model.csv').to_dict('records')+readc('outputs/v2/stage2/secondary_outcomes.csv').to_dict('records')[:3]
    fig,axes=plt.subplots(2,2,figsize=(10,8),layout='constrained')
    for ax,r in zip(axes.flat,rs):scatter(ax,data,r['outcome'],r,('Primary: ' if r['outcome']=='HRD_total' else 'Secondary component: ')+r['outcome'])
    fig.suptitle('TCGA: frozen candidate activity and genomic HRD\nShaded bands: pointwise HC3 95% uncertainty in the OLS mean',fontsize=13);save(fig,FIGURES[1])
    data=readc('outputs/v2/stage4/scanb_analysis_cohort.csv');rs=readc('outputs/v2/stage4/scanb_primary_results.csv').to_dict('records');fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for ax,r in zip(axes,rs):scatter(ax,data,r['outcome'],r,r['outcome'])
    fig.suptitle('SCAN-B: independent-cohort unadjusted replication\nPointwise HC3 95% bands; OLS bands are not constrained to probability bounds',fontsize=12);save(fig,FIGURES[2])
    fig=plt.figure(figsize=(13,10),layout='constrained');gs=fig.add_gridspec(3,3,height_ratios=[1,1.9,.9]);s2=readj('outputs/v2/stage2/stage2_summary.json');s4=readj('outputs/v2/stage4/stage4_summary.json')
    groups=[('TCGA HRD',s2['primary'],s2['proliferation_adjusted']),('SCAN-B HRDetect',s4['primary'][0],s4['proliferation_adjusted'][0]),('SCAN-B scarHRD',s4['primary'][1],s4['proliferation_adjusted'][1])]
    for j,(title,a,b) in enumerate(groups):
        ax=fig.add_subplot(gs[0,j])
        for y,r,color in [(1,a,blue),(0,b,orange)]:ax.errorbar(r['beta'],y,xerr=[[r['beta']-r['hc3_ci_low']],[r['hc3_ci_high']-r['beta']]],fmt='o',color=color,capsize=3);ax.text(.02,.91 if y==1 else .04,f"HC3 P={r['hc3_p']:.3g}",transform=ax.transAxes,fontsize=9)
        ax.axvline(0,color=gray,lw=.8);ax.set(yticks=[0,1],yticklabels=['Proliferation adjusted','Unadjusted'],ylim=(-.6,1.6),title=title,xlabel='Beta (endpoint-specific units)')
        loo=readc('outputs/v2/stage2/leave_one_gene_out.csv') if j==0 else readc('outputs/v2/stage4/scanb_leave_one_gene_out.csv')
        if j:loo=loo[loo.outcome==a['outcome']]
        ax=fig.add_subplot(gs[1,j]);ax.errorbar(loo.beta,np.arange(len(loo)),xerr=[loo.beta-loo.hc3_ci_low,loo.hc3_ci_high-loo.beta],fmt='o',color=blue,capsize=2);ax.axvline(0,color=gray,lw=.8);ax.set(yticks=np.arange(len(loo)),yticklabels=loo.omitted_gene,xlabel='Leave-one-out beta; HC3 95% CI',title='Omitted gene');ax.invert_yaxis()
    ax=fig.add_subplot(gs[2,:]);ax.axis('off');control=[]
    for name,perm,rnd in [('TCGA',s2['permutation'],s2['random_gene_null']),('SCAN-B HRDetect',s4['permutation'][0],s4['matched_random'][0]),('SCAN-B scarHRD',s4['permutation'][1],s4['matched_random'][1])]:control.append([name,f"{perm['empirical_p']:.4g}",str(perm['iterations']),f"{rnd['empirical_p']:.4g}",f"{rnd['percentile']:.2f}%"])
    tab=ax.table(cellText=control,colLabels=['Endpoint','Permutation P','Permutations','Matched-random P','Null percentile'],loc='center',cellLoc='center');tab.auto_set_font_size(False);tab.set_fontsize(10);tab.scale(1,1.7)
    fig.suptitle('Specificity and robustness: attenuation is part of the result\nProliferation sets differ between TCGA and SCAN-B; each panel retains its own effect scale.',fontsize=13);save(fig,FIGURES[3])
    select=e[(e.primary_secondary=='primary')|((e.cohort=='SCAN-B')&(e.analysis_type=='proliferation'))|((e.cohort=='METABRIC')&(e.predictor=='candidate_mean_z'))|e.cohort.isin(['GSE25066','DepMap'])|((e.cohort=='GSE176078')&(e.predictor=='composite')&(e.population=='TNBC')&(e.endpoint=='malignant_vs_nonmalignant_combined'))].copy()
    fig,ax=plt.subplots(figsize=(15,max(6,.48*len(select)+2)));ax.axis('off');cells=[]
    for r in select.to_dict('records'):
        effect=f"{fmt(r['effect'])} [{fmt(r['ci_low'])}, {fmt(r['ci_high'])}]" if pd.notna(r['ci_low']) else fmt(r['effect'])
        cells.append([r['cohort'],str(r['endpoint']).replace('malignant_vs_nonmalignant_combined','malignant - nonmalignant'),int(r['N']),effect,'OR' if 'odds' in r['effect_measure'] else 'HR' if 'hazard' in r['effect_measure'] else 'count' if r['cohort']=='DepMap' else 'difference' if r['cohort']=='GSE176078' else 'beta',fmt(r['robust_p'] if pd.notna(r['robust_p']) else r['p']),fmt(r['fdr']),'Adjusted' if 'adjust' in r['analysis_type'].lower() or r['analysis_type']=='proliferation' else 'Unadjusted/context'])
    tab=ax.table(cellText=cells,colLabels=['Cohort','Endpoint','N','Effect [95% CI]','Scale','HC3 P / P','BH q','Model'],loc='center',cellLoc='left',colWidths=[.10,.23,.045,.23,.065,.085,.065,.16]);tab.auto_set_font_size(False);tab.set_fontsize(9);tab.scale(1,1.9)
    for i,r in enumerate(select.to_dict('records'),1):
        p=r['fdr'] if pd.notna(r['fdr']) else r['robust_p'] if pd.notna(r['robust_p']) else r['p'];color='#E4F0F1' if pd.notna(p) and p<.05 else '#F4E9E2'
        for j in range(8):tab[(i,j)].set_facecolor(color)
    ax.set_title('Cross-cohort evidence: positive and negative results\nSeparate effect measures; no pooled numerical axis or meta-analysis',fontsize=14,pad=20)
    fig.text(.5,.035,'CIs are conventional unless the scale is a paired bootstrap difference. HC3 P is used when available.\nDepMap q is the minimum genome-wide adjusted q; effect is the number of discoveries. GSE25066 uses four genes only.',ha='center',fontsize=9);save(fig,FIGURES[4])

def numbers(e,c):
    records=[]
    fields=['N','beta','effect','ci_low','ci_high','p','robust_p','adjusted_p','robust_adjusted_p','fdr','r_squared','spearman_rho','events','robust_ci_low','robust_ci_high','AUC','iterations','percentile','tested_dependencies']
    for r in e.to_dict('records'):
        for field in fields:
            if field in r and pd.notna(r[field]):
                source=r['source_output'];locator=str(r['source_row'])
                if r['cohort']=='GSE25066' and field=='events' and r['endpoint']=='pCR':source='outputs/v2/stage3/stage3_summary.json';locator='gse25066.pcr_events'
                if r['cohort']=='DepMap':locator=('all rows' if 'genome-wide' in r['endpoint'] else 'rows annotated DDR')+'; N=joint models; effect=count(fdr_adjusted<0.05); fdr=min(fdr_adjusted); tested_dependencies=row count'
                records.append(dict(section='Abstract' if r['primary_secondary']=='primary' else 'Discussion' if r['supports_main_claim']=='NO_SUPPORT_AT_THRESHOLD' else 'Results',analysis_id=r['analysis_id'],analysis_name=f"{r['cohort']}: {r['endpoint']}; {r['analysis_type']}; {r['predictor']}",field=field,exact_value=float(r[field]),publication_value=fmt(r[field]),source_output=source,source_row=locator,evidence_field=field))
    for i,r in c.iterrows():
        for field in ['eligible_N','expression_or_context_N','primary_analysis_N']:records.append(dict(section='Methods',analysis_id='COHORT_'+str(i),analysis_name=r['cohort']+' '+field,field=field,exact_value=float(r[field]),publication_value=str(r[field]),source_output=r['source_output'],source_row=r['source_field'],evidence_field=field))
    out=pd.DataFrame(records);out.to_csv(FINAL/'manuscript_numbers.csv',index=False)
    text='# Authoritative manuscript numbers\n\nGenerated from frozen outputs. No manuscript text has been edited. Exact values are IEEE double-precision representations read from saved CSV/JSON. Publication rounding is three significant digits. Conventional and robust intervals/P are distinct; beta units vary by endpoint. The machine-readable companion is manuscript_numbers.csv.\n'
    for section in ['Abstract','Methods','Results','Discussion']:
        text+='\n## '+section+'\n\n| Analysis | Field | Exact | Publication | Source output and row |\n|---|---|---:|---:|---|\n'
        for r in records:
            if r['section']==section:text+=f"| {r['analysis_id']}: {r['analysis_name']} | {r['field']} | {r['exact_value']!r} | {r['publication_value']} | `{r['source_output']}`; row/key {r['source_row']} |\n"
    (FINAL/'MANUSCRIPT_NUMBERS.md').write_text(text,encoding='utf-8')

def documents(summary,ledger):
    text='''# Final integrated computational methods

This document describes the executed frozen Stages 1-4. This freeze assembles existing results and redraws their relationships; it does not search for new associations, reselect cohorts, tune thresholds, download data or edit a manuscript. Source rows are traceable through master_evidence_table.csv. All figures use every complete case used by their corresponding existing model; HC3 pointwise bands are reconstructed for display without changing coefficients or inference.

## Clinical definition and expression

TCGA TNBC requires negative clinical ER, PR and HER2 IHC, with positive HER2 ISH conflicts excluded. Unknown/equivocal receptors are not negative. Code 01 and concordant sample annotation identify primary tumors. ESR1/PGR/ERBB2 expression was not used to define TNBC. GIPS was not used. The frozen cohort has 114 patients; 113 have usable RNA, and 103 have the primary HRD measurement. The missing RNA patient remains in the clinical cohort. These complete-case observations are the executed discovery analysis, with incomplete coverage documented; no reference refit is made during this freeze.

GDC augmented STAR counts use tpm_unstranded, transformed log2(TPM+1), GRCh38/GENCODE v36. Ensembl identifiers, symbols, sample IDs and available aliquots are retained. Exact clinical/sample matching precedes a highest-assigned-read-count technical tie-break, then lexical sample/aliquot/file ordering. One patient has two eligible files; no biological-quality guarantee follows from this technical rule. Missing full aliquot barcodes use an explicitly flagged sample/file UUID key. The original exon matrix is not used.

Nine fixed genes: HORMAD1, HORMAD2, STAG3, REC8, SMC1B, SYCP2, SYCP3, MSH4, MSH5. Each is z-standardized with population SD (ddof=0) within the specified cohort reference and averaged equally; the median and HORMAD1/STAG3/REC8 z-scores are sensitivities. All nine are required for the main TCGA, SCAN-B, METABRIC and DepMap score. Canonical comparisons are separate and never define outcomes. No gene or coefficient is chosen from validation associations. Clinical prediction of a new individual is not established by cohort restandardization.

## TCGA discovery and independent genomic measurements

Published TCGA HRD_total is the primary DNA scar outcome, with HRD_LOH, LST, TAI, aneuploidy, FGA, LOH_fraction, ploidy and WGD secondary. Publication tables are TCGA.HRD_withSampleID.txt, ABSOLUTE_scores.tsv, seg_based_scores.tsv and TCGA_mastercalls.abs_tables_JSedit.fixed.txt. WGD uses the curated genome-doubling field, not an optimized ploidy threshold. Exact primary sample keys and source identities are retained; discordant duplicate outcomes are excluded. Candidate RNA never defines DNA outcomes.

OLS reports beta, conventional t CI/P, R-squared, Pearson and Spearman correlations; HC3 sandwich t CI/P is separate. M1 adds age/stage, M2 purity, M3 proliferation, with same-complete-case unadjusted comparisons. The most adjusted model has N=97. Missing BRCA1/2 status is not wild type. Binary models use logistic regression and predictor scaling to one model-specific SD. Pseudo-R-squared is not variance explained by OLS.

## Independent SCAN-B replication

The 235 clinical-TNBC cases in Nacer et al. 2026 (doi:10.1186/s13058-026-02325-5) supply direct WGS HRDetect_probability and scarHRD_score. Supplement sheets merge on PD_ID. RBA maps to GEO through exact S###### and m/m2/m3/m4 components; k/k2 is ignored. No numeric RBA-to-F match or m-variant rescue occurs. Titles ending repl are excluded; ambiguous nonreplicate mapping aborts. There are 225 unique matches, 10 unmatched patients, zero ambiguities and nine excluded replicate records.

GSE96058 expression is already log2(FPKM+0.1), as confirmed by local SOFT metadata. It is not double-logged or mislabeled TPM. Standardization uses all 225 mapped cases independently of endpoint availability. Compressed data are streamed; only requested genes and matched-null neighborhoods are retained. The primary full-nine-gene score is unchanged. All nine genes are present.

HRDetect and scarHRD are two co-primary continuous endpoints; Holm adjusts this two-test family, separately for conventional and HC3 P. They are independent of candidate RNA, but not statistically independent of one another. SBS3, SV3 and delMH are secondary WGS measures sharing HRDetect inputs. Published HRD/HRP labels are used verbatim for secondary logistic regression; OR is per one SD and AUC is apparent, not held-out. Fractional-logit binomial quasi-likelihood with HC0 covariance accommodates continuous probabilities including boundaries; it is a sensitivity, not a replacement for OLS. Age, proliferation and PAM50 sensitivities are separate; PAM50 never selects the cohort. Patient-level BRCA status and purity are unavailable.

## Proliferation definitions

TCGA and METABRIC use BIRC5, CCNB1, CDC20, NUF2, CEP55, NDC80, MKI67, PTTG1, RRM2, TYMS, UBE2C. DepMap and primary SCAN-B adjustment use MKI67, PCNA, CCNB1, CCNE1, MCM2, MCM4, MCM6, TOP2A, BIRC5, CDC20, CDC6. Both are disjoint from candidates. SCAN-B also records the earlier TCGA/METABRIC definition as a named sensitivity. Stage 3 does not have one uniform proliferation set: its patient and functional modules differ. No attempt is made to silently harmonize or rerun them at freeze.

## METABRIC

The 258 clinical TNBC primary tumors require ER IHC negative, PR and HER2 negative, no conflicting positive reported ER. Source typo normalization precedes classification. Duplicate expression symbols are averaged in log2 Illumina HT12-v3 intensity, including both STAG3 rows. All nine genes are restandardized within 258 patients. The genomic endpoint is the fraction of unique assayed genes with absolute discrete CNA at least one, using median duplicate DNA calls. It is gene-weighted, not length-weighted FGA or direct HRD. Primary OLS/HC3 association is accompanied by age/proliferation adjustment. OS and RFS are continuous-score Cox models with Efron ties and HR per reference-cohort predictor SD.

## GSE25066

There are 178 clinically defined TNBC samples: er_status_ihc=N, pr_status_ihc=N, her2_status=N. Expression-based receptor rescue fields are not used. Completed local GPL96 annotation maps only STAG3, REC8, SYCP2 and MSH4. All unambiguous represented probes are averaged per gene, with within-cohort z-standardization. This previously executed four-gene score is explicitly reduced-platform; it is not the frozen full-nine-gene validation. Observed pCR logistic analysis includes 112 patients and 37 events. DRFS Cox analysis uses 114 patients and 37 events, source years and explicit event flags, Efron ties and score per cohort SD. Neither endpoint supports the clinical claim. The DRFS residual/log-time diagnostic flags possible PH departure. No full-nine-gene result is inferred.

## Single-cell context

GSE176078 uses published annotations, exact orig.ident donor/sample units and barcode checks. Sparse integer counts are summed per donor/class, with at least 20 cells and a positive library total; library sums agree with source nCount_RNA. Pseudobulk log2(CPM+1) is z-standardized over eligible TNBC donor-by-disjoint-class profiles. The primary malignant-minus-pooled-nonmalignant comparison is paired by eight eligible TNBC donor/sample units from ten TNBC units (26 all-breast units). Use paired two-sided Wilcoxon and a 2,000 paired-donor bootstrap mean-difference CI, seed 20260908. Cells are not independent replicates. No single-cell HRD state is inferred from the RNA score.

## Functional context

DepMap eligibility uses literal curated TNBC/triple-negative non-expression annotations. Of 34 annotated TNBC models, 25 have joint expression and CRISPR coverage (102 breast models in metadata). The frozen nine-gene score is standardized within these models. Every measurable genome-wide CRISPR gene effect is tested by OLS with and without the functional proliferation set. Negative beta means stronger dependency at higher signature. Conventional t inference and separate genome-wide BH corrections are retained; DDR annotations are applied after testing and do not create a restricted new FDR family. No dependency survives adjusted genome-wide FDR. No added drug screen or therapeutic efficacy analysis is part of this freeze.

## Multiplicity, empirical controls and diagnostics

TCGA secondary genomic endpoints, individual candidates and canonical analyses have separately documented BH families. METABRIC composite CNA/OS/RFS and individual-gene families are separate; GSE25066 pCR/DRFS form a two-test BH family. Single-cell composite contrasts and individual genes use separate families by population. SCAN-B uses BH for three secondary genomic endpoints, two binary labels, and separate nine/13-gene families per primary endpoint. Leave-one-out and median are declared sensitivities; nominal significance is not an additional independent confirmation.

TCGA: 5,000 matched nine-gene sets (seed 20260907), 5,000 outcome-label permutations (20260908), 2,000 paired-case bootstrap resamples (20260909). SCAN-B: 5,000 matched sets (20260909), 10,000 permutations per endpoint (20260919/20260920), 2,000 bootstrap resamples (20260929/20260930). Fixed standardization is retained through resampling. Permutations and random signatures use absolute OLS t and plus-one empirical P. Random controls select one gene per candidate from 200 nearest mean/variance percentile neighbors without replacement within a set. Eligibility requires expression at least 0.1 in at least 20% of the reference and nonzero variation. TCGA uses verified protein-coding genes; SCAN-B uses reported expressed symbols without a verified coding filter. These backgrounds are not interchangeable. Matching does not balance every biological property and unadjusted nulls do not control confounding.

Shapiro, Breusch-Pagan, quadratic RESET, residual plots, leverage and Cook distance characterize model limitations. Primary observations are not removed; exclusions above Cook distance 4/N are separately labeled sensitivities. All leave-one-gene-out effects are reported without choosing an omission. Figure 4 shows proliferation attenuation and the full leave-one-out series.

## Reconciliation and reproducibility

Earlier Stage 1 methods contain a provisional-score caution and an erroneous no-duplicates sentence despite one duplicate group. Later frozen analyses use the documented 113-sample reference; coverage and technical identity limitations remain explicit. Earlier Stage 3 narrative describes GPL96/DepMap as unavailable; the current saved result tables and executed adapters establish their later completion. Those original documents are preserved, and the final methods resolve the stale narrative rather than deleting negative results. The Stage 4 earlier-set terminology is refined here: it refers to the TCGA/METABRIC set, whereas Stage 3 DepMap already used the other set.

The package preserves pre-freeze scientific files, checks historical hashes where recorded, runs existing tests and read-only frozen-cohort validation, and records code/software/input/output SHA256 values. Newly established hashes prove integrity from this freeze onward; absent historical hashes cannot prove earlier raw-file history. TableS2 is the full machine model table; a giant DepMap LaTeX rendering is intentionally omitted. TableS4 CSV retains all mapping columns. Tables with compact TeX views retain full fields in CSV. Freeze timestamps are UTC; historical output modification times are provenance approximations, not inferred execution dates. No internet verification or new data acquisition was performed during this final assembly.
'''
    text+='\nSoftware at freeze: Python '+platform.python_version()+'; '+', '.join(k+' '+v for k,v in summary['software'].items())+'.\n'
    (ROOT/'docs/V2_FINAL_METHODS.md').write_text(text,encoding='utf-8')
    limitations='# Final limitations\n'
    for title,indices in [('Methodological limitations',[0,7,8,9,10,11,12,13]),('Biological limitations',[1,5]),('Validation limitations',[2,3,4,6]),('Clinical limitations',[14])]:
        limitations+='\n## '+title+'\n\n'+'\n\n'.join(LIMITATIONS[i] for i in indices)+'\n'
    limitations+='\nThe candidate signature replicated against WGS-derived HRD measures in SCAN-B in unadjusted analyses, but the association attenuated substantially and was no longer statistically significant after proliferation adjustment. Therefore the analysis supports a relationship with HRD-associated tumor state but does not establish a proliferation-independent HRD biomarker. Negative tests do not prove equivalence or biological absence.\n'
    (ROOT/'docs/V2_LIMITATIONS.md').write_text(limitations,encoding='utf-8')
    review='''# Final evidence review

## What is novel?

The contribution of this package is the frozen nine-gene, clinically defined TNBC analysis across independent genomic, clinical and functional evidence streams, including direct WGS replication and explicit negative controls. Literature-wide priority or novelty is not established by a computational freeze; no new literature search was conducted.

## What replicated?

TCGA HRD association replicates unadjusted in the independent SCAN-B cohort against HRDetect and scarHRD. METABRIC adds an orthogonal CNA proxy association. TCGA and SCAN-B leave-one-out analyses do not require HORMAD1 to retain nominal robust support. Independence of cohorts does not establish independence from proliferation or independence between genomic endpoints.

## What failed?

SCAN-B and METABRIC adjustment attenuates the candidate association. SCAN-B scarHRD matched-random P exceeds 0.05. GSE25066 reduced-platform pCR/DRFS, METABRIC composite survival, DepMap genome-wide/DDR and single-cell malignant specificity are unsupported. Missing drug evidence is unavailable, not a negative drug experiment. These outcomes remain visible in TableS3 and Figure 5.

## Defensible conclusions and removals

A fixed meiotic/cohesin transcriptional signature is an observational correlate of HRD-associated tumor state in clinically defined TNBC, with unadjusted independent-cohort WGS replication and important proliferation confounding. Remove proliferation-independent clinical biomarker, chemotherapy prediction, survival prediction, selective dependency, malignant-specificity and causal HRD claims. HORMAD1 is strong individually, but the assertion that it explains the entire composite is unsupported. Claim ledger A-K governs exact wording.

## What remains observational, and what needs experiments?

All expression-to-genomic and expression-to-dependency associations remain observational. A causal proposal would require controlled candidate perturbation and rescue, proliferation-matched controls, orthogonal HR repair assays, temporal genome-instability measurements and reproducible selective genetic/drug effects in appropriate TNBC models. Those experiments are not performed or promised by this package.

## Strongest manuscript framing

Frame the work as multicohort characterization of a meiotic/cohesin transcriptional correlate of HRD-associated tumor state, with proliferation and negative functional/clinical validation setting its limits. Lead with precise replicated effects, then show attenuation and null findings at equal visibility. Do not present retrospective association as a deployable diagnostic or mechanism.

## Realistic publication scope

Editorial judgment from the present evidence favors a specialist cancer genomics, translational bioinformatics or rigorous computational oncology venue accepting transparent observational replication and negative results. A high-selectivity mechanistic or clinical-validation positioning is not justified by these data alone. No named journal ranking, current impact-factor comparison or acceptance probability is asserted. Fit depends on reproducibility, incremental contribution, cohort independence and candid limitations; acceptance remains uncertain.
'''
    (ROOT/'docs/V2_FINAL_EVIDENCE_REVIEW.md').write_text(review,encoding='utf-8')

def manifest(summary,test):
    paths=set()
    for folder in ['cohort','expression','predictors','stage2','stage3','stage4']:
        paths.update(p for p in (ROOT/'outputs/v2'/folder).rglob('*') if p.is_file())
    paths.update((ROOT/'src').glob('*.py'));paths.update((ROOT/'tests').glob('*.py'));paths.update([ROOT/'requirements.txt',ROOT/'environment.yml',ROOT/'config/data_sources.yaml'])
    clinical=readj('outputs/v2/cohort/cohort_summary.json');paths.add(Path(clinical['clinical_file']))
    paths.update(Path(p) for p in readc('outputs/v2/expression/expression_sample_manifest.csv').source_path)
    paths.update(Path(p) for p in readj('outputs/v2/stage4/preparation_manifest.json')['input_sha256'])
    paths.update(Path(p) for p in readj('outputs/v2/stage3/single_cell_source.json')['source_stat'])
    original=ROOT.parent.parent
    paths.update(original/p for p in ['dataset/brca_metabric_clinical_data.tsv','dataset/GSE25066_series_matrix.txt.gz','results/metabric_extracted/brca_metabric/data_mrna_illumina_microarray.txt','results/metabric_extracted/brca_metabric/data_cna.txt','results/metabric_extracted/brca_metabric/meta_mrna_illumina_microarray.txt'])
    urls={}
    for r in readj('config/data_sources.yaml')['resources']:
        p=ROOT/'data/external'/r['id']/r['filename']
        if p.exists() and (r['id'].startswith('tcga_') and r['id']!='tcga_rna' or r['id'].startswith('depmap_') or r['id']=='gpl96'):
            if r['id'] in ['tcga_hrd','tcga_aneuploidy','tcga_fga','tcga_absolute','depmap_effect','depmap_model','depmap_expression','depmap_default','depmap_profiles','depmap_condition','depmap_readme','gpl96']:paths.add(p);urls[str(p)]=r['url']
    def rec(p):return dict(path=str(p),sha256=sha256(p),bytes=p.stat().st_size,modified_utc=datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat(),source_url=urls.get(str(p)),historical_integrity='See historical_hash_validation.csv; otherwise hash first established at freeze')
    inputs=[rec(p) for p in sorted(paths) if p.exists()]
    outputs=[rec(p) for p in sorted(FINAL.rglob('*')) if p.is_file() and p.name not in ['reproducibility_manifest.json','validation_report.json']]+[rec(ROOT/'docs'/n) for n in DOCS]+[rec(ROOT/'scripts'/n) for n in ['build_final_package.py','validate_final_package.py']]
    commit=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True)
    savej(FINAL/'reproducibility_manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),package_versions=summary['software'],candidate_genes=list(CANDIDATE_GENES),cohort_Ns=summary['cohort_Ns'],random_seeds=summary['random_seeds'],git_commit=commit.stdout.strip() if commit.returncode==0 else None,git_worktree_status=subprocess.run(['git','status','--porcelain'],cwd=ROOT,capture_output=True,text=True).stdout,known_test_failures=test['failures'],inputs=inputs,outputs=outputs,source_accessions=['TCGA-BRCA','GSE96058','GSE25066','GPL96','GSE176078','METABRIC','DepMap'],source_references=['https://doi.org/10.1186/s13058-026-02325-5','https://gdc.cancer.gov/about-data/publications/pancanatlas'],timestamp_note='Freeze execution timestamp recorded; input mtime is not asserted to be original analysis time.',exclusions_from_self_hash=['reproducibility_manifest.json','validation_report.json']))

def main():
    for p in [FINAL,FINAL/'figures',FINAL/'tables']:p.mkdir(parents=True,exist_ok=True)
    print('Validating package',flush=True);test=freeze_and_tests()
    e,dep=evidence();e.to_csv(FINAL/'master_evidence_table.csv',index=False);ledger=claims(e);ledger.to_csv(FINAL/'claim_ledger.csv',index=False);c=cohorts()
    s2=readj('outputs/v2/stage2/stage2_summary.json');s3=readj('outputs/v2/stage3/stage3_summary.json');s4=readj('outputs/v2/stage4/stage4_summary.json')
    summary=dict(status='assembled',cohort_Ns={'TCGA_clinical':114,'TCGA_RNA':113,'TCGA_HRD':103,'SCANB_supplement':235,'SCANB_matched':225,'METABRIC':258,'GSE25066_clinical':178,'GSE25066_pCR':112,'GSE25066_DRFS':114,'DepMap_curated_TNBC':34,'DepMap_analysis':25,'single_cell_TNBC_units':10,'single_cell_TNBC_pairs':8},candidate_genes=list(CANDIDATE_GENES),TCGA=s2,SCANB=s4,external_primary=e[(e.stage==3)&(e.predictor.isin(['candidate_mean_z','composite'])|e.cohort.isin(['GSE25066','DepMap']))].to_dict('records'),major_limitations=LIMITATIONS,software={n:version(n) for n in ['numpy','pandas','scipy','statsmodels','matplotlib','openpyxl']},random_seeds={'TCGA_matched':20260907,'TCGA_permutation':20260908,'TCGA_bootstrap':20260909,'single_cell_bootstrap':20260908,'SCANB_matched':20260909,'SCANB_permutation':[20260919,20260920],'SCANB_bootstrap':[20260929,20260930]},claim_ledger_counts=ledger.allowed_strength.value_counts().to_dict(),test_results=test,analysis_source='Frozen Stage 1-4 outputs; no new inferential tests')
    savej(FINAL/'final_statistical_summary.json',summary)
    print(f'Assembling evidence: {len(e)}',flush=True)
    make_tables(e,dep,c,test);figures(e);numbers(e,c);documents(summary,ledger)
    (FINAL/'figures/FIGURE_LEGENDS.md').write_text('''# Figure legends

Figure 1. Frozen study design and cohort-specific analysis populations. Arrows indicate evidence streams, not patient overlap. Single-cell N counts donor/sample pairs; DepMap N counts models.

Figure 2. TCGA candidate mean z-score against primary HRD_total and three secondary scar components. Every endpoint-complete case is shown. Orange lines use the saved OLS coefficient; shading is pointwise HC3 95% uncertainty in the fitted mean reconstructed from the identical saved cohort. No outlier removal or axis clipping occurs. Components share the primary scar biology; they are not independent primary tests.

Figure 3. Independent-cohort SCAN-B WGS outcomes against frozen candidate activity, N=225 for both. Shading is pointwise HC3 mean uncertainty, not a prediction interval. Primary OLS and its bands are not constrained to 0-1, exposing the limitation of a linear probability-mean relationship. Fractional-logit sensitivity is in the adjusted-results table. Reported HC3 P values are raw; two-endpoint Holm values are tabulated separately.

Figure 4. Top: unadjusted and proliferation-adjusted saved effects with HC3 95% CIs, each endpoint on its own scale. Middle: all nine leave-one-gene-out effects, with no chosen omission. Bottom: saved unadjusted permutation and matched-random comparisons. Gene-set nulls use 5,000 sets in each cohort. Differences in background filters and proliferation sets are described in methods. The loss of adjusted SCAN-B significance prevents a general independence claim.

Figure 5. Evidence plot using separate effect measures rather than an inappropriate pooled effect axis. Blue rows have support at the displayed BH q or HC3/conventional P threshold; tan rows do not. CIs are conventional except the donor-paired bootstrap difference. OR/HR are per predictor SD; continuous betas are per original signature unit and retain outcome units. DepMap count is discoveries at adjusted genome-wide FDR<0.05, and the q shown is the minimum genome-wide adjusted q, not a global hypothesis P. GSE25066 is four-gene reduced-platform evidence. Nonsignificance is not equivalence.
''',encoding='utf-8')
    manifest(summary,test)
    print('Package assembled',flush=True)

if __name__=='__main__':main()
