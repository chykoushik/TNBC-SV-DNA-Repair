import sys,json,csv
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from src.utils import sha256
from src.config import CANDIDATE_GENES
from build_final_package import FIGURES,TABLES,DOCS,manuscript_files
FINAL=ROOT/'outputs/v2/final'
FROZEN=('HORMAD1','HORMAD2','STAG3','REC8','SMC1B','SYCP2','SYCP3','MSH4','MSH5')

def forbidden_fields(obj):
    if isinstance(obj,dict):return any('gips' in str(k).lower() or forbidden_fields(v) for k,v in obj.items())
    if isinstance(obj,list):return any(forbidden_fields(v) for v in obj)
    return False

def validate():
    checks=[]
    def check(name,ok,detail=''):checks.append(dict(check=name,passed=bool(ok),detail=detail))
    required=['master_evidence_table.csv','claim_ledger.csv','final_statistical_summary.json','MANUSCRIPT_NUMBERS.md','manuscript_numbers.csv','reproducibility_manifest.json','freeze_baseline.json','test_results.json']
    expected=[FINAL/n for n in required]+[FINAL/'figures'/f'{n}.{ext}' for n in FIGURES for ext in ['pdf','png']]+[FINAL/'tables'/f'{n}.csv' for n in TABLES]+[ROOT/'docs'/n for n in DOCS]
    check('all_expected_files',all(p.is_file() and p.stat().st_size>0 for p in expected),';'.join(str(p) for p in expected if not p.exists()))
    if not checks[-1]['passed']:return finish(checks)
    readj=lambda p:json.loads(p.read_text(encoding='utf-8'))
    summary=readj(FINAL/'final_statistical_summary.json');ledger=pd.read_csv(FINAL/'claim_ledger.csv');e=pd.read_csv(FINAL/'master_evidence_table.csv');numbers=pd.read_csv(FINAL/'manuscript_numbers.csv');manifest=readj(FINAL/'reproducibility_manifest.json');baseline=readj(FINAL/'freeze_baseline.json')
    check('candidate_genes_exact',tuple(CANDIDATE_GENES)==FROZEN and tuple(summary['candidate_genes'])==FROZEN and tuple(manifest['candidate_genes'])==FROZEN)
    genes=pd.read_csv(FINAL/'tables/TableS1_gene_definitions.csv');check('candidate_table_exact',tuple(genes.loc[genes.role=='candidate','gene'])==FROZEN)
    clinical=pd.read_csv(ROOT/'outputs/v2/cohort/tcga_tnbc_patients.csv');s4=readj(ROOT/'outputs/v2/stage4/stage4_summary.json');s2=readj(ROOT/'outputs/v2/stage2/stage2_summary.json')
    check('TCGA_N_frozen',summary['cohort_Ns']['TCGA_clinical']==len(clinical)==s2['clinical_n'] and summary['cohort_Ns']['TCGA_HRD']==s2['primary']['N'])
    check('SCANB_N_stage4',summary['cohort_Ns']['SCANB_matched']==s4['mapping']['uniquely_mapped_n'] and set(e.loc[(e.cohort=='SCAN-B')&(e.primary_secondary=='primary'),'N'])=={s4['mapping']['uniquely_mapped_n']})
    bad=[]
    for p in FINAL.rglob('*.csv'):
        with p.open(encoding='utf-8',newline='') as f:header=next(csv.reader(f))
        if any('gips' in k.lower() for k in header):bad.append(str(p))
    check('no_retired_scientific_fields',not bad and not forbidden_fields(summary),str(bad))
    expected_labels={'A':'SUPPORTED_WITH_QUALIFICATION','B':'SUPPORTED','C':'NOT_SUPPORTED','D':'NOT_SUPPORTED','E':'NOT_SUPPORTED','F':'NOT_SUPPORTED','G':'NOT_SUPPORTED','H':'NOT_SUPPORTED','I':'SUPPORTED','J':'NOT_SUPPORTED','K':'NOT_SUPPORTED'}
    check('claim_ledger_negative_evidence',ledger.set_index('claim_id').allowed_strength.to_dict()==expected_labels and ledger.evidence_against_or_limiting.notna().all())
    negative=pd.read_csv(FINAL/'tables/TableS3_negative_results.csv')
    check('negative_results_visible',{'GSE25066','GSE176078','DepMap','SCAN-B'}.issubset(set(negative.cohort)) and len(negative[(negative.cohort=='SCAN-B')&(negative.analysis_type=='proliferation')])==2)
    changed=[p for p,h in baseline['protected'].items() if not Path(p).exists() or sha256(p)!=h]
    check('protected_inputs_outputs_unchanged',not changed,str(changed))
    current={str(p):sha256(p) for p in manuscript_files()};check('no_manuscript_modified',current==baseline['manuscripts'],f"{len(current)} manuscript/bibliography files checked; newly created files also checked")
    failed=[]
    for r in manifest['inputs']+manifest['outputs']:
        p=Path(r['path'])
        if not p.exists() or sha256(p)!=r['sha256']:failed.append(str(p))
    check('manifest_hashes_match',not failed,str(failed))
    check('all_numbers_have_source_references',numbers.source_output.notna().all() and all((ROOT/p).is_file() for p in numbers.source_output.unique()))
    failed=[];indexed=e.set_index('analysis_id')
    for r in numbers.to_dict('records'):
        if not r['analysis_id'].startswith('COHORT_'):
            saved=indexed.loc[r['analysis_id'],r['evidence_field']]
            if not np.isclose(float(saved),r['exact_value'],rtol=1e-12,atol=1e-14):failed.append(r['analysis_id']+':'+r['field'])
    check('manuscript_numbers_equal_evidence',not failed,str(failed[:20]))
    # Source comparison
    failed=[]
    for r in e.to_dict('records'):
        if r['source_row']<0 or not r['source_output'].endswith('.csv'):continue
        source=pd.read_csv(ROOT/r['source_output']).iloc[int(r['source_row'])]
        for dest,key in [('beta','beta'),('p','p'),('robust_p','hc3_p'),('N','N'),('fdr','fdr')]:
            if key in source and pd.notna(source[key]) and not np.isclose(float(r[dest]),float(source[key]),rtol=1e-12,atol=1e-14):failed.append(r['analysis_id']+':'+dest)
    check('evidence_equal_frozen_sources',not failed,str(failed))
    test=readj(FINAL/'test_results.json');check('scientific_tests_pass',test['scientific_tests_passed'],f"{test['passed']} passed; {test['failed']} infrastructure/acquisition failures")
    check('five_vector_and_high_resolution_figures',all((FINAL/'figures'/f'{n}.pdf').read_bytes().startswith(b'%PDF') and (FINAL/'figures'/f'{n}.png').read_bytes().startswith(b'\x89PNG') for n in FIGURES))
    return finish(checks)

def finish(checks):
    report=dict(status='PASSED' if all(r['passed'] for r in checks) else 'FAILED',checks=checks)
    (FINAL/'validation_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2));return report

if __name__=='__main__':sys.exit(0 if validate()['status']=='PASSED' else 1)
