import json
import platform
from importlib.metadata import version
from pathlib import Path
import numpy as np
import pandas as pd
from .config import PROJECT_ROOT,CANDIDATE_GENES,COMPARISON_GENES
from .utils import sha256
from .covariates import load_frozen,add_clinical,PROLIFERATION_GENES,PROLIFERATION_SOURCE
from .outcomes import load_outcomes,OUTCOMES
from .statistics import fit_model,bh,unavailable
from .null_models import matched_random,permutation,bootstrap,leave_one_out,SEED
from .figures_stage2 import render

def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value

def savejson(path,data):path.write_text(json.dumps(clean(data),indent=2,allow_nan=False)+'\n',encoding='utf-8')

def add_fdr(rows):
    q=bh([r.get('p',np.nan) for r in rows])
    for row,value in zip(rows,q):row['fdr']=value
    return rows

def run():
    root=PROJECT_ROOT;dest=root/'outputs/v2/stage2';dest.mkdir(parents=True,exist_ok=True)
    frozen=['outputs/v2/cohort/tcga_tnbc_patients.csv','outputs/v2/predictors/candidate_signature.csv',
            'outputs/v2/expression/expression_sample_manifest.csv','outputs/v2/expression/tcga_tnbc_expression.tsv.gz','outputs/v2/predictors/standardization_parameters.json']
    frozen_hashes={p:sha256(root/p) for p in frozen}
    protocol=dict(primary_outcome='HRD_total',primary_predictor='candidate_mean_z',secondary_outcomes=list(OUTCOMES[1:]),
                  random_sets=5000,permutations=5000,bootstrap=2000,seed=SEED,proliferation_genes=PROLIFERATION_GENES,
                  minimum_complete_n=10,adjusted_minimum_patients_per_parameter=10,vif_warning=5,
                  random_neighbors=200,random_matching_features=['reference mean log2(TPM+1)','reference variance log2(TPM+1)'],
                  outcome_dichotomization=False,signature_weights='frozen equal weights; no optimization')
    savejson(dest/'prespecified_protocol.json',protocol)
    data,expression,z,annotation,missing_prolif=load_frozen(root)
    data=add_clinical(data);outcomes,availability,conflicts=load_outcomes(root,data)
    data=pd.concat([data,outcomes],axis=1)
    data.to_csv(dest/'analysis_cohort.csv',index=False)
    overlap=[dict(outcome=o,clinical_n=len(data),rna_predictor_n=int(data.candidate_mean_z.notna().sum()),outcome_n=int(data[o].notna().sum()),
                  overlap_n=int(data[[o,'candidate_mean_z']].dropna().shape[0]),status='available' if data[o].notna().any() else 'unavailable') for o in OUTCOMES]
    pd.DataFrame(overlap).to_csv(dest/'outcome_overlap.csv',index=False)
    pd.DataFrame(availability).to_csv(dest/'outcome_resources.csv',index=False)
    pd.DataFrame(conflicts,columns=['patient_id','resource','status']).to_csv(dest/'outcome_matching_issues.csv',index=False)
    tables={};models={}
    primary,models['primary'],diagnostics=fit_model(data)
    primary['hypothesis_family']='primary; no secondary-family correction';tables['primary']=[primary]
    for row in diagnostics:row['patient_id']=data.loc[row['row_index'],'patient_id']
    pd.DataFrame(diagnostics,columns=['row_index','patient_id','fitted','residual','cooks_distance','leverage','influential']).to_csv(dest/'primary_diagnostics.csv',index=False)
    nested=[('Model 0',()),('Model 1',('age','stage')),('Model 2',('age','stage','purity')),('Model 3',('age','stage','purity','proliferation'))]
    adjusted=[fit_model(data,covariates=c,label=label)[0] for label,c in nested]
    # Complete-case companions
    for label,covs in nested[1:]:
        subset=data.dropna(subset=['HRD_total','candidate_mean_z',*covs])
        adjusted.append(fit_model(subset,label=label+' matched-case M0')[0])
    tables['adjusted']=adjusted
    tables['secondary']=add_fdr([fit_model(data,outcome=o,label='Secondary outcome',binary=o=='WGD')[0] for o in OUTCOMES[1:]])
    single=[]
    for g in CANDIDATE_GENES:
        r=fit_model(data,predictor=g+'_expr_z',label='Exploratory single candidate')[0];r['gene']=g;single.append(r)
    tables['single']=add_fdr(single)
    canonical=[]
    for g in COMPARISON_GENES:
        single=fit_model(data,predictor=g+'_expr_z',label='Canonical single gene')[0];single['gene']=g;canonical.append(single)
        full,fullfit,_=fit_model(data,covariates=(g+'_expr_z',),label='Candidate beyond canonical '+g);full['gene']=g
        if fullfit is not None:
            reduced_data=data.loc[fullfit.model.data.row_labels]
            reduced,rf,_=fit_model(reduced_data,predictor=g+'_expr_z',label='Reduced canonical model')
            if rf is not None:
                f,p,df=fullfit.compare_f_test(rf);full.update(incremental_r_squared=full['r_squared']-reduced['r_squared'],incremental_f_p=float(p))
        canonical.append(full)
    for label in ['Canonical single gene']:
        add_fdr([r for r in canonical if r['model']==label])
    add_fdr([r for r in canonical if r['model'].startswith('Candidate beyond')])
    tables['canonical']=canonical
    prolif,models['proliferation'],_=fit_model(data,outcome='proliferation',label='Candidate versus proliferation')
    pa=fit_model(data,covariates=('proliferation',),label='HRD adjusted for proliferation')[0]
    tables['proliferation']=[prolif,pa]
    logo=[]
    for gene,score in leave_one_out(z).items():
        d=data.copy();d['leave_one_out_score']=d.patient_id.map(score)
        r=fit_model(d,predictor='leave_one_out_score',label='Leave one gene out')[0];r['omitted_gene']=gene;logo.append(r)
    tables['logo']=logo
    sensitivity=[]
    for label,covs in [('Proliferation',('proliferation',)),('Purity',('purity',)),('BRCA1/2',('BRCA1_altered','BRCA2_altered')),('FGA',('FGA',)),('PAM50',('PAM50',))]:
        sensitivity.append(fit_model(data,covariates=covs,label='Adjust '+label)[0])
    for predictor in ['candidate_median_z','HORMAD1_z','STAG3_z','REC8_z']:
        sensitivity.append(fit_model(data,predictor=predictor,label='Prespecified alternative predictor')[0])
    for g in CANDIDATE_GENES:
        sensitivity.append(fit_model(data,covariates=(g+'_expr_z',),label='Exploratory candidate conditional on '+g)[0])
    basal=data[data.PAM50.astype(str).str.lower().isin(['basal','basal-like'])]
    sensitivity.append(fit_model(basal,label='PAM50 basal within clinical TNBC')[0])
    unaltered=data[(data.BRCA1_altered==0)&(data.BRCA2_altered==0)]
    sensitivity.append(fit_model(unaltered,label='Exclude BRCA1/2 altered (requires known negatives)')[0])
    excluded=[r['row_index'] for r in diagnostics if r['influential']]
    sensitivity.append(fit_model(data.drop(index=excluded),label="Exclude Cook's distance >4/N")[0])
    tables['sensitivity']=sensitivity
    null=pd.DataFrame(columns=['iteration','genes','beta','abs_t','mean_matching_distance','max_matching_distance'])
    null_summary=dict(status='not_run_missing_HRD',planned_iterations=5000,iterations=0,seed=SEED,empirical_p=None,percentile=None)
    perm_summary=dict(status='not_run_missing_HRD',planned_iterations=5000,iterations=0,seed=SEED+1,empirical_p=None)
    boot_summary=dict(status='not_run_missing_HRD',planned_iterations=2000,iterations=0,seed=SEED+2,ci_low=None,ci_high=None)
    if primary['status']=='estimated':
        d=data[['patient_id','candidate_mean_z','HRD_total']].dropna().sort_values('patient_id')
        ids=d.patient_id.tolist();x=d.candidate_mean_z.to_numpy();y=d.HRD_total.to_numpy()
        # Protein-coding background
        protein=set(annotation.loc[annotation.gene_type=='protein_coding','gene_symbol'])
        keep=[g for g in expression.index if g in protein]
        null,null_summary=matched_random(expression.loc[keep],z.loc[keep],ids,y,x)
        perm_summary,ts=permutation(x,y);pd.DataFrame({'t_statistic':ts}).to_csv(dest/'permutation_distribution.csv',index=False)
        boot_summary,betas=bootstrap(x,y);pd.DataFrame({'beta':betas}).to_csv(dest/'bootstrap_distribution.csv',index=False)
    tables['null']=null
    filenames={'primary':'primary_hrd_model.csv','adjusted':'adjusted_models.csv','secondary':'secondary_outcomes.csv','single':'single_gene_results.csv',
               'canonical':'canonical_comparison.csv','proliferation':'proliferation_analysis.csv','logo':'leave_one_gene_out.csv','sensitivity':'sensitivity_analysis.csv'}
    for key,filename in filenames.items():pd.DataFrame(tables[key]).to_csv(dest/filename,index=False)
    null.to_csv(dest/'random_gene_null.csv',index=False)
    savejson(dest/'random_gene_null_summary.json',null_summary);savejson(dest/'permutation_summary.json',perm_summary);savejson(dest/'bootstrap_summary.json',boot_summary)
    covariate_counts={c:int(data[c].notna().sum()) for c in ['age','stage','purity','proliferation','PAM50','BRCA1_altered','BRCA2_altered','FGA']}
    # Algorithmic evidence
    classification='WEAK';basis='No local independent HRD outcomes: hypothesis untested, not a demonstrated null association.'
    if primary['status']=='estimated':
        basis='Primary association absent or robustness insufficient.'
        if primary['p']<.05 and primary['hc3_p']<.05:
            classification='PROMISING';basis='Within-cohort association; external replication not assessed.'
            major=[pa]+[r for r in adjusted if r['model']=='Model 3']+[r for r in sensitivity if r['model']=="Exclude Cook's distance >4/N"]
            robust=all(r['status']=='estimated' and r['hc3_p']<.05 and np.sign(r['beta'])==np.sign(primary['beta']) and not r.get('coefficient_stability_warning',False) for r in major)
            if robust and primary['r_squared']>=.1 and null_summary.get('empirical_p',1)<.05 and perm_summary.get('empirical_p',1)<.05:
                classification='STRONG';basis='Meets prespecified internal association/robustness criteria; not causal or externally replicated.'
            if (pa['status']=='estimated' and pa['hc3_p']>=.05) or null_summary.get('empirical_p',1)>=.05:
                classification='WEAK';basis='Proliferation adjustment or matched random-set comparison does not support specificity.'
    estimated_genes=[r for r in tables['single'] if r['status']=='estimated']
    strongest=min(estimated_genes,key=lambda r:r['fdr']) if estimated_genes else None
    summary=dict(clinical_n=len(data),rna_n=int(data.candidate_mean_z.notna().sum()),HRD_overlap_n=primary['N'],primary=primary,
                 proliferation_adjusted=pa,most_adjusted=next(r for r in adjusted if r['model']=='Model 3'),
                 random_gene_null=null_summary,permutation=perm_summary,bootstrap=boot_summary,strongest_gene=strongest,
                 secondary_surviving_fdr=[r['outcome'] for r in tables['secondary'] if r.get('fdr',1)<.05],
                 classification=classification,classification_basis=basis,primary_test_performed=primary['status']=='estimated',
                 outcome_resources=availability,covariate_nonmissing=covariate_counts,missing_proliferation_genes=missing_prolif,
                 stage1_unchanged=True,source_hashes=frozen_hashes,python=platform.python_version(),
                 packages={name:version(name) for name in ['numpy','pandas','scipy','statsmodels','matplotlib']},
                 caveats=['BRCA1/2 alteration status cannot be inferred as wild-type from a reduced mutation list',
                          'No independent-cohort replication is tested in this stage','Missing results are not negative results'])
    savejson(dest/'stage2_summary.json',summary)
    render(root,data,tables,models,null_summary)
    write_docs(root,summary,prolif)
    for p,h in frozen_hashes.items():
        if sha256(root/p)!=h:raise RuntimeError('Stage 1 source changed: '+p)
    print(json.dumps(clean(dict(HRD_overlap_n=primary['N'],primary=primary,proliferation=prolif,classification=classification,covariate_counts=covariate_counts)),indent=2))
    return summary

def write_docs(root,s,prolif):
    docs=root/'docs'
    methods=f'''# V2 Stage 2 and core Stage 3 methods

Run `python -m src.stage2` from V2. Stage 1 is consumed unchanged: 114 clinically defined TNBC patients, 113 RNA matches, frozen equal-weight candidate predictors and log2(TPM+1) STAR matrix. No Stage 1 reconstruction, source edits or biological downloads occur. Runtime/package provenance and Stage 1 hashes are in the output manifest. Modules separate outcomes, covariates, statistics, null models, figures and orchestration.

## Independent outcomes and matching

Only four exact publication files registered in `config/data_sources.yaml` are accepted: `TCGA.HRD_withSampleID.txt` (HRD, hrd-loh, lst1, ai1); `ABSOLUTE_scores.tsv` (AS, LOH_frac_altered); `seg_based_scores.tsv` (frac_altered); `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (purity, ploidy, Genome doublings). Their source links are the [NCI PanCancer Atlas](https://gdc.cancer.gov/about-data/publications/pancanatlas) and [PanImmune resource](https://gdc.cancer.gov/about-data/publications/panimmune). Genome doublings >=1 is a binary WGD call derived from the curated field, never from a ploidy threshold. HRD-LOH is a scar count, not LOH fraction. Missing resources remain unavailable; RNA/CNA gene panels do not substitute for these outcomes.

Match the frozen patient and primary-tumor sample key, preferring an exact RNA aliquot when available. Record each genomic source sample ID separately for every outcome. Discordant genomic duplicate rows are excluded; identical duplicate values can collapse with all source IDs retained. The analysis cohort is left-joined to all frozen clinical patients. Models use endpoint-specific complete cases, not an intersection of all modalities. `outcome_overlap.csv`, `outcome_resources.csv` and `outcome_matching_issues.csv` expose attrition and conflicts.

## Predictors and proliferation

The primary predictor remains candidate_mean_z. Secondary alternatives are candidate_median_z, HORMAD1_z, STAG3_z and REC8_z. Individual z-scores for all nine candidates and the 13 separate canonical genes are computed from the frozen discovery matrix using population SD across all 113 RNA patients. The reconstructed candidate mean must numerically agree with Stage 1 within 1e-8; the source predictor is never overwritten or trained against HRD.

Proliferation uses the published 11-gene index: {', '.join(PROLIFERATION_GENES)}. Source: [Nielsen et al., 2010, proliferation-index gene table]({PROLIFERATION_SOURCE}). The score here is the mean of per-gene discovery-cohort z-scores of log2(TPM+1), an explicitly documented RNA-seq implementation of that gene set, not the proprietary PAM50 risk assay. All 11 genes are required; candidate genes are disjoint. Both candidate/proliferation correlation and HRD ~ candidate + proliferation are specified.

Age and pathologic stage use the original clinical file matched by exact clinical sample ID. Stage is categorical I–IV, with letter substages collapsed. PAM50 is a sensitivity covariate only and never defines TNBC. Purity comes only from the independent ABSOLUTE table. Reliable BRCA1/2 alteration/negative denominators, germline and biallelic status are not established locally; these covariates remain missing rather than labeling absent somatic rows wild-type.

## Models and diagnostics

Primary OLS: HRD_total ~ candidate_mean_z, continuous HRD. Report conventional beta, SE, Student-t 95% CI, P, R2, Pearson r/P and Spearman rho/P. HC3 sandwich SEs and t-based CIs/P are additional robustness estimates. Diagnostics include quadratic RESET (linearity/model form), Shapiro-Wilk residual test, Breusch–Pagan heteroscedasticity test, residual-versus-fitted/Q-Q plots, leverage and Cook's distance. These tests are descriptive and do not choose the predictor or primary estimator.

Nested M0 unadjusted; M1 age+stage; M2 adds purity; M3 adds proliferation. Models require >=10 complete cases, full-rank design and nonconstant outcome; adjusted models additionally require >=10 cases per estimated parameter including the intercept and stage dummy terms. VIF>5 is flagged, not automatically optimized away. Same-complete-case M0 companions distinguish attrition from adjustment. Missing/unstable models are not evidence of independence.

Secondary continuous outcomes use the same OLS/HC3/Spearman framework. Binary WGD uses logistic regression, with candidate_mean_z scaled by its SD in that model's complete-case cohort; report OR and 95% CI per 1 SD. At least 10 patients in each class are required; convergence/separation warnings block a result. Binary pseudo-R2 is distinct from OLS R2. BH correction covers the actually estimable eight-outcome secondary family, separate from the unadjusted primary hypothesis. Conventional model P is the prespecified FDR input; HC3 is reported separately. Missing tests retain missing q-values, never zero. The nine individual candidate tests form their own BH family. Canonical-gene single tests and 13 candidate-beyond-canonical tests are separate exploratory BH families. No multiplicity-controlled external replication claim is made.

For each canonical gene, fit a gene-only model and candidate+canonical model on identical complete cases; report candidate coefficient, incremental R2 and partial F-test when feasible. All nine candidate-conditioned models and nine leave-one-gene-out means are exploratory diagnostics of single-gene dependence, not an opportunity to select a new signature. Sensitivities include BRCA exclusion/adjustment when valid, purity, proliferation, FGA, PAM50, basal-only, Cook's distance >4/N exclusion and all prespecified alternative predictors. Report all, not the strongest subset.

## Empirical robustness

Random gene-set null: 5,000 nine-gene sets, seed {SEED}. Use unique-symbol protein-coding genes, excluding candidates and canonical comparisons. Require TPM>=0.1 in at least 20% of reference samples and log-expression variance>1e-8. For each candidate slot find up to 200 nearest eligible genes in mean/variance percentile space across the frozen RNA discovery reference; sample one without replacement within each set. Record selected genes and matching distances. This is approximate expression matching, not proof that all other genomic properties are balanced. Standardize random genes over the same discovery reference and fit the same M0 on the same HRD-complete patients. Compare absolute t, not raw betas whose scales differ. Empirical two-sided P=(1+null |t|>=observed |t|)/(B+1); percentile is the proportion of null |t| below observed. Sets may recur across iterations; sampling is not outcome-dependent.

Independently permute HRD patient labels 5,000 times (seed {SEED+1}) on fixed complete-case predictor/outcome pairs; use the same absolute-t statistic and plus-one P. This is an unadjusted exchangeability test, not a confounding-adjusted test. Bootstrap complete-case patients 2,000 times (seed {SEED+2}), sampling pairs with replacement and retaining the frozen predictor; report percentile 2.5–97.5% beta CI and invalid replicate counts. This conditional bootstrap does not re-estimate the discovery gene standardization.

No empirical procedure runs when the primary outcome is missing. Output JSONs report planned and actual iteration counts separately. Empty distribution CSVs do not represent null findings. Figures 1–5 are explicitly marked unavailable if HRD/CIN data are missing; Figure 6 can still show the observed expression/proliferation relationship. Outputs are vector PDF and 300-dpi PNG. Genomic forest effects are standardized for continuous outcome comparability; WGD OR stays separately tabulated.

## Evidence label and limitations

STRONG is reserved for conventional and HC3 primary P<0.05, primary R2>=0.10 (a declared working effect-size threshold, not a validated clinical threshold), stable same-direction M3/proliferation/Cook-exclusion HC3 P<0.05, empirical random/permutation P<0.05 and no VIF warning. PROMISING denotes a within-cohort signal with incomplete robustness, without claiming independent replication. WEAK covers absent/unsupported/surrogate/unexceptional evidence. If no HRD test is possible, WEAK is an administrative evidence-availability label; the hypothesis is **untested**, not disproven. The label does not modify the hypothesis or gene weights. Observational associations are not causal, and no external cohort is analyzed here.
'''
    (docs/'V2_METHODS_STAGE2.md').write_text(methods,encoding='utf-8')
    unavailable=s['primary']['status']!='estimated'
    primary_text='Not estimable: HRD-total overlap N=0. Beta, SE, CI, P, R2, HC3, adjusted effects, permutation and bootstrap results are unavailable.' if unavailable else json.dumps(clean(s['primary']),indent=2)
    result=f'''# V2 Stage 2 results

**Classification: {s['classification']}.** {s['classification_basis']}

{primary_text}

The frozen cohort remains {s['clinical_n']} patients with {s['rna_n']} complete expression predictors. Independent HRD/CIN outcomes are unavailable locally for all nine requested endpoints. This is not a negative association result. No primary HRD model, adjusted HRD model, single-gene HRD model or genomic robustness distribution was estimated in this run. No secondary outcome can be said to survive or fail FDR. No strongest HRD-associated candidate gene can be identified, and composite dependence on an individual gene cannot be evaluated.

The available expression-only proliferation analysis used N={prolif['N']}: Pearson r={prolif['pearson_r']:.6g}, P={prolif['pearson_p']:.6g}; Spearman rho={prolif['spearman_rho']:.6g}, P={prolif['spearman_p']:.6g}; proliferation-on-candidate beta={prolif['beta']:.6g}, 95% CI [{prolif['ci_low']:.6g}, {prolif['ci_high']:.6g}], R2={prolif['r_squared']:.6g}. These are expression relationships, not evidence for HRD or independence from proliferation. The HRD-adjusted result remains unavailable.

Four missing, already-configured public resources block the central tests: `tcga_hrd` (280,336 bytes), `tcga_aneuploidy` (456,006), `tcga_fga` (454,662), and `tcga_absolute` (901,812), totaling **2,092,816 bytes (~2.09 MB)**. They were reported rather than downloaded. Their stable source URLs and release metadata remain in `config/data_sources.yaml`. If acquired separately, place each in its configured `data/external/<resource>/` directory and rerun Stage 2. Original files and Stage 1 do not need to be rebuilt.

Additional material limitations: valid BRCA1/2 alteration-negative denominators and biallelic/germline status are unavailable; purity/FGA require the missing tables; one frozen clinical patient lacks RNA; many RNA aliquot barcodes are unresolved in Stage 1 although case and primary sample/file identities are retained. No clinical membership was changed. Covariate nonmissing counts: {s['covariate_nonmissing']}.

All requested output tables were written. Outcome-dependent table rows explicitly state unavailable with N=0 and missing estimates. Matched-random, permutation and bootstrap actual iteration counts are zero; their reusable implementations were tested on synthetic data, not run with fabricated HRD. Figures 1–5 are availability notices rather than fabricated scientific results. Figure 6 displays the observed proliferation relationship and identifies the unavailable adjusted test. No circular enrichment or retired composite analysis was performed.
'''
    (docs/'V2_RESULTS_STAGE2.md').write_text(result,encoding='utf-8')
    if not unavailable:
        result=f'''# V2 Stage 2 results

Classification: **{s['classification']}**. {s['classification_basis']}

Primary HRD overlap: {s['HRD_overlap_n']}. Primary estimates:

```json
{json.dumps(clean(s['primary']),indent=2)}
```

Proliferation-adjusted estimates:

```json
{json.dumps(clean(s['proliferation_adjusted']),indent=2)}
```

Matched-random empirical P: {s['random_gene_null'].get('empirical_p')}; percentile: {s['random_gene_null'].get('percentile')}.
Permutation empirical P: {s['permutation'].get('empirical_p')}. Bootstrap CI: [{s['bootstrap'].get('ci_low')}, {s['bootstrap'].get('ci_high')}].
Secondary outcomes surviving BH FDR: {s['secondary_surviving_fdr']}. Strongest single candidate by FDR: {s['strongest_gene']}.

All adjusted, canonical, individual-gene, leave-one-out and sensitivity estimates, including non-significant estimates, are in their named CSVs. P>=0.05 does not support an association at the specified threshold; unavailable estimates cannot establish independence. These are observational within-cohort results, not causal effects or independent-cohort replication.
Covariate nonmissing counts: {s['covariate_nonmissing']}. Source availability: {s['outcome_resources']}. Candidate definitions were unchanged.
'''
        (docs/'V2_RESULTS_STAGE2.md').write_text(result,encoding='utf-8')

if __name__=='__main__':run()
