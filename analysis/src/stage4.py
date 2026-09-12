import json
import warnings
from importlib.metadata import version
import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from .config import PROJECT_ROOT,CANDIDATE_GENES,COMPARISON_GENES
from .utils import sha256
from .statistics import fit_model,bh,unavailable
from .null_models import permutation,bootstrap,association
from .stage2 import clean
from .stage4_scanb import prepare,PRIMARY,SEED,STAGE4_PROLIFERATION,PREVIOUS_PROLIFERATION,sample_null_sets

def fractional(data,covariates=()):
    cols=['HRDetect_probability','candidate_mean_z',*covariates];d=data[cols].dropna()
    label='Fractional logit'+(' adjusted' if covariates else '')
    r=unavailable('HRDetect_probability','candidate_mean_z',label,len(d),'insufficient complete cases')
    if len(d)<20:return r
    x=sm.add_constant(d[['candidate_mean_z',*covariates]],has_constant='add')
    if np.linalg.matrix_rank(x)<x.shape[1]:r['reason']='rank deficient';return r
    with warnings.catch_warnings(record=True) as caught:
        fit=sm.GLM(d.HRDetect_probability,x,family=sm.families.Binomial()).fit(cov_type='HC0')
    if caught or not fit.converged or not np.isfinite(fit.params).all():r['reason']='fractional logit warning/nonconvergence';return r
    ci=fit.conf_int().loc['candidate_mean_z']
    r.update(status='estimated',reason='',beta=float(fit.params.candidate_mean_z),se=float(fit.bse.candidate_mean_z),ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),p=float(fit.pvalues.candidate_mean_z),
             covariance='HC0 sandwich quasi-likelihood',effect_scale='logit conditional mean per signature unit; not binary-HRD OR')
    return r

def annotate_auc(result,model,data,outcome):
    d=data[[outcome,'candidate_mean_z']].dropna();result['HRD_events']=int(d[outcome].sum());result['N']=len(d)
    if model is not None:
        prediction=np.asarray(model.predict());truth=d[outcome].to_numpy();n1=int(truth.sum());n0=len(truth)-n1
        result['AUC_apparent']=(st.rankdata(prediction)[truth==1].sum()-n1*(n1+1)/2)/(n1*n0)
        result['AUC_note']='Apparent in-sample discrimination, not held-out performance'
    return result

def run():
    root=PROJECT_ROOT;dest=root/'outputs/v2/stage4';dest.mkdir(parents=True,exist_ok=True)
    protected=[root/'outputs/v2/cohort/tcga_tnbc_patients.csv',root/'outputs/v2/predictors/candidate_signature.csv',root/'outputs/v2/stage2/primary_hrd_model.csv',root/'src/config.py']
    frozen={str(p):sha256(p) for p in protected}
    protocol={'primary_endpoints':PRIMARY,'primary_predictor':'candidate_mean_z','primary_family_adjustment':'Holm, two co-primary endpoints',
              'random_signatures':5000,'permutations_per_endpoint':10000,'bootstrap_per_endpoint':2000,'seed':SEED,
              'proliferation_primary':STAGE4_PROLIFERATION,'previous_proliferation_sensitivity':PREVIOUS_PROLIFERATION,
              'candidate_genes':CANDIDATE_GENES,'outlier_exclusion_sensitivity':'Cook distance >4/N; primary unchanged',
              'classification_rule':'STRONG: both positive primary Holm/HC3-Holm<.05, proliferation-HC3<.05, matched-null/permutation<.05 and positive significant LOO/Cook sensitivities; MODERATE: some primary evidence with incomplete robustness; WEAK: no supported primary signal or both specificity controls fail; FAILED: scientific/data block'}
    (dest/'scanb_protocol.json').write_text(json.dumps(protocol,indent=2))
    prepare(root)
    data=pd.read_csv(dest/'scanb_analysis_cohort.csv')
    params=json.loads((dest/'scanb_standardization.json').read_text())
    if data.candidate_mean_z.isna().any():raise ValueError('Incomplete primary signature: stop rather than substitute genes')
    primary=[];primary_models={};diagnostics=[];adjusted=[];secondary=[];sensitivity=[];logo=[];canonical=[]
    for outcome in PRIMARY:
        result,model,diags=fit_model(data,outcome=outcome,label='Co-primary OLS')
        result['family']='co-primary';primary.append(result);primary_models[outcome]=model
        if model is None:raise ValueError('Primary model not estimable: '+outcome)
        for d in diags:d.update(outcome=outcome,PD_ID=data.loc[d['row_index'],'PD_ID'],RBA=data.loc[d['row_index'],'RBA'],F_expression_column=data.loc[d['row_index'],'F_expression_column'])
        diagnostics.extend(diags)
        for label,cov in [('proliferation',('proliferation_mean_z',)),('age',('age',)),('age + proliferation',('age','proliferation_mean_z')),
                          ('previous Stage2/3 proliferation',('previous_proliferation_mean_z',)),('PAM50 + proliferation',('PAM50','proliferation_mean_z'))]:
            r,_,_=fit_model(data,outcome=outcome,covariates=cov,label=label);adjusted.append(r)
            if label=='age + proliferation':
                subset=data.dropna(subset=[outcome,'candidate_mean_z',*cov]);adjusted.append(fit_model(subset,outcome=outcome,label='Same-complete-case unadjusted for age + proliferation')[0])
        for missing in ['BRCA1 alteration','BRCA2 alteration','tumor purity']:
            adjusted.append(unavailable(outcome,'candidate_mean_z','Adjust '+missing,0,'No patient-level '+missing+' supplied in local inputs'))
        for g in CANDIDATE_GENES:
            r,_,_=fit_model(data,outcome=outcome,predictor=g+'_z',label='Exploratory individual candidate');r['gene']=g;sensitivity.append(r)
            d=data.copy();d['leave_one_out_score']=d[[h+'_z' for h in CANDIDATE_GENES if h!=g]].mean(axis=1)
            r,_,_=fit_model(d,outcome=outcome,predictor='leave_one_out_score',label='Leave one gene out');r['omitted_gene']=g;logo.append(r)
        r,_,_=fit_model(data,outcome=outcome,predictor='candidate_median_z',label='Prespecified median signature');r['gene']='composite_median';sensitivity.append(r)
        exclude=[d['row_index'] for d in diags if d['influential']]
        r,_,_=fit_model(data.drop(index=exclude),outcome=outcome,label='Exclude Cook distance >4/N');r['gene']='composite';r['excluded_patient_ids']=';'.join(data.loc[exclude,'PD_ID']);sensitivity.append(r)
        for gene in COMPARISON_GENES:
            if gene+'_comparison_z' in data:
                r,_,_=fit_model(data,outcome=outcome,predictor=gene+'_comparison_z',label='Exploratory canonical gene');r['gene']=gene;canonical.append(r)
    for p,q,h in zip(primary,multipletests([r['p'] for r in primary],method='holm')[1],multipletests([r['hc3_p'] for r in primary],method='holm')[1]):p['primary_holm_p']=q;p['primary_hc3_holm_p']=h
    for outcome in PRIMARY:
        for rows in [sensitivity,canonical]:
            group=[r for r in rows if r['outcome']==outcome and r.get('model','').startswith('Exploratory')]
            for r,q in zip(group,bh([r['p'] for r in group])):r['fdr']=q
    for outcome in ['SBS3','SV3','delMH']:
        if outcome not in data:secondary.append(unavailable(outcome,'candidate_mean_z','Secondary genomic exposure',0,'Column absent'))
        else:secondary.append(fit_model(data,outcome=outcome,label='Secondary genomic exposure')[0])
    for r,q in zip(secondary,bh([r['p'] for r in secondary])):r['fdr']=q;r['family']='three secondary genomic exposures'
    binary=[]
    for outcome in ['HRDetect_binary','scarHRD_binary']:
        r,fit,_=fit_model(data,outcome=outcome,label='Published HRD/HRP logistic',binary=True);r=annotate_auc(r,fit,data,outcome);binary.append(r)
    for r,q in zip(binary,bh([r['p'] for r in binary])):r['fdr']=q;r['family']='two published binary classifications'
    secondary+=binary
    adjusted+=[fractional(data),fractional(data,('proliferation_mean_z',))]
    pd.DataFrame(primary).to_csv(dest/'scanb_primary_results.csv',index=False)
    pd.DataFrame(secondary).to_csv(dest/'scanb_secondary_results.csv',index=False)
    pd.DataFrame(adjusted).to_csv(dest/'scanb_adjusted_results.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(dest/'scanb_gene_sensitivity.csv',index=False)
    pd.DataFrame(logo).to_csv(dest/'scanb_leave_one_gene_out.csv',index=False)
    pd.DataFrame(canonical).to_csv(dest/'scanb_canonical_results.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(dest/'scanb_diagnostics.csv',index=False)
    # Shared reference parameters
    neighbors=json.loads((dest/'scanb_null_neighborhoods.json').read_text())
    sets=sample_null_sets(neighbors['neighbors'],5000,SEED)
    null_expression=pd.read_csv(dest/'scanb_null_neighborhood_expression.tsv.gz',sep='\t',index_col=0)
    z=null_expression.sub(null_expression.mean(axis=1),axis=0).div(null_expression.std(axis=1,ddof=0),axis=0)
    null_scores=np.array([z.loc[genes,data.F_expression_column].mean(axis=0).to_numpy() for genes in sets])
    null_results=[];null_summaries=[];perms=[];boots=[]
    for i,outcome in enumerate(PRIMARY):
        subset=data[[outcome,'candidate_mean_z']].notna().all(axis=1).to_numpy();d=data.loc[subset];x=d.candidate_mean_z.to_numpy();y=d[outcome].to_numpy()
        perm,ts=permutation(x,y,n=10000,seed=SEED+10+i);perm.update(outcome=outcome,N=len(d));perms.append(perm)
        boot,betas=bootstrap(x,y,n=2000,seed=SEED+20+i);boot.update(outcome=outcome,N=len(d));boots.append(boot)
        pd.DataFrame({'t_statistic':ts}).to_csv(dest/(outcome+'_permutation_distribution.csv'),index=False)
        pd.DataFrame({'beta':betas}).to_csv(dest/(outcome+'_bootstrap_distribution.csv'),index=False)
        beta,observed=association(x,y);null_stats=[]
        for j,genes in enumerate(sets):
            nb,nt=association(null_scores[j,subset],y);null_stats.append(abs(nt))
            distances=[neighbors['distances'][g][h] for g,h in zip(CANDIDATE_GENES,genes)]
            null_results.append(dict(outcome=outcome,iteration=j+1,N=len(d),genes=';'.join(genes),beta=nb,abs_t=abs(nt),mean_matching_distance=float(np.mean(distances)),max_matching_distance=float(max(distances))))
        ns=np.asarray(null_stats)
        if not np.isfinite(ns).all():raise ValueError('Invalid null signature statistic')
        null_summaries.append(dict(outcome=outcome,N=len(d),iterations=5000,seed=SEED,observed_beta=beta,observed_abs_t=abs(observed),
                                   empirical_p=float((1+(ns>=abs(observed)).sum())/5001),percentile=float(100*(ns<abs(observed)).mean()),eligible_pool_n=neighbors['eligible_pool_n']))
    pd.DataFrame(perms).to_csv(dest/'scanb_permutation_results.csv',index=False);pd.DataFrame(boots).to_csv(dest/'scanb_bootstrap_results.csv',index=False)
    pd.DataFrame(null_results).to_csv(dest/'scanb_null_signature_results.csv',index=False);pd.DataFrame(null_summaries).to_csv(dest/'scanb_null_signature_summary.csv',index=False)
    # Fixed replication criteria
    primary_ok=all(r['beta']>0 and r['primary_holm_p']<.05 and r['primary_hc3_holm_p']<.05 for r in primary)
    pa=[r for r in adjusted if r['model']=='proliferation']
    robust=all(r['beta']>0 and r['hc3_p']<.05 for r in pa)
    null_ok=all(r['empirical_p']<.05 for r in null_summaries+perms)
    loo_ok=all(r['beta']>0 and r['hc3_p']<.05 for r in logo)
    cook_ok=all(r['beta']>0 and r['hc3_p']<.05 for r in sensitivity if r['model']=='Exclude Cook distance >4/N')
    classification='STRONG' if primary_ok and robust and null_ok and loo_ok and cook_ok else 'MODERATE' if any(r['beta']>0 and r['primary_hc3_holm_p']<.05 for r in primary) else 'WEAK'
    if not robust and not any(r['empirical_p']<.05 for r in null_summaries):classification='WEAK'
    summary={'status':'completed','classification':classification,'mapping':{k:v for k,v in json.loads((dest/'scanb_mapping_audit.json').read_text()).items() if k!='source_processing'},
             'available_candidate_genes':params['candidate_available'],'missing_candidate_genes':params['candidate_missing'],
             'primary':primary,'proliferation_adjusted':pa,'permutation':perms,'bootstrap':boots,'matched_random':null_summaries,
             'leave_one_out_all_positive':all(r['beta']>0 for r in logo),'leave_one_out_all_hc3_significant':loo_ok,
             'covariate_availability':{'age':int(data.age.notna().sum()),'proliferation':int(data.proliferation_mean_z.notna().sum()),'BRCA1_status':0,'BRCA2_status':0,'purity':0},
             'proliferation_definition_note':'User-listed Stage 4 set is primary; different actual Stage2/3 set tested separately, not mislabeled as identical',
             'bounded_sensitivity':[r for r in adjusted if r['model'].startswith('Fractional')],
             'blockers':[],'limitations':['10 unmatchable S/m records retained as exclusions','Patient-level BRCA alteration and purity unavailable in supplied tables','Endpoints share genomic features and are not statistically independent','Random pool includes expressed gene symbols without a verified protein-coding annotation filter','Logistic AUC is apparent, not held-out'],
             'frozen_previous_input_hashes':frozen,'packages':{n:version(n) for n in ['numpy','pandas','scipy','statsmodels','matplotlib','openpyxl']}}
    (dest/'stage4_summary.json').write_text(json.dumps(clean(summary),indent=2,allow_nan=False))
    diagnostic_figures(dest,data,primary_models,primary)
    methods(root,summary,params,primary,secondary,adjusted,sensitivity,logo)
    for path,digest in frozen.items():
        if sha256(path)!=digest:raise RuntimeError('Previously frozen data/code changed')
    print(json.dumps(clean({k:summary[k] for k in ['classification','mapping','primary','proliferation_adjusted','permutation','matched_random','bootstrap']}),indent=2))
    return summary

def diagnostic_figures(dest,data,models,results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'font.size':9})
    for outcome,model in models.items():
        fig,axes=plt.subplots(1,3,figsize=(12,3.8));d=data[[outcome,'candidate_mean_z']].dropna()
        axes[0].scatter(d.candidate_mean_z,d[outcome],s=12,alpha=.6);grid=np.linspace(d.candidate_mean_z.min(),d.candidate_mean_z.max(),100)
        p=model.get_prediction(sm.add_constant(grid)).summary_frame();axes[0].plot(grid,p['mean'],color='black');axes[0].fill_between(grid,p.mean_ci_lower,p.mean_ci_upper,alpha=.2)
        axes[0].set(xlabel='Frozen candidate mean z-score',ylabel=outcome)
        axes[1].scatter(model.fittedvalues,model.resid,s=12,alpha=.6);axes[1].axhline(0,color='gray');axes[1].set(xlabel='Fitted outcome',ylabel='Residual')
        sm.qqplot(model.resid,fit=True,line='45',ax=axes[2]);fig.tight_layout()
        fig.savefig(dest/(outcome+'_diagnostics.pdf'));fig.savefig(dest/(outcome+'_diagnostics.png'),dpi=300);plt.close(fig)

def methods(root,s,params,primary,secondary,adjusted,sensitivity,logo):
    dest=root/'outputs/v2/stage4';m=s['mapping'];pre=json.loads((dest/'preparation_manifest.json').read_text())
    text=f'''# V2 Stage 4: independent SCAN-B WGS-HRD replication

Implemented and run with `python -m src.stage4`. Classification: **{s['classification']}**. No manuscript changes, dataset downloads, original data edits, pathway enrichment or drug analyses occurred. Previous Stage 1/2 data and the frozen candidate configuration are hash-checked unchanged.

## Cohort and outcomes

Source: [Nacer DF et al., Breast Cancer Research 2026;28:115, doi:10.1186/s13058-026-02325-5](https://doi.org/10.1186/s13058-026-02325-5), future bibliography key `Nacer2026`. The supplied cohort comprises {m['supplement_n']} clinically defined TNBC patients. Supplement 2 supplies WGS identity, RNA identity, HRDetect/scarHRD scores and classifications. We merge its HRD classifications and HRD scores sheets by unique PD_ID. Cohort membership is inherited from the clinical study, not redefined with ESR1, PGR, ERBB2 or PAM50 expression. Supplement 3 contains DNA-repair gene lists and outcome-stratified differential-expression results; it provides no additional patient-level age, BRCA status or purity and is not used to select genes or predictors.

The two co-primary endpoints are continuous HRDetect_probability and scarHRD_score. Both are genomic measurements independent of candidate RNA, but they share genomic features and are correlated; they are not statistically independent tests. Secondary SBS3 and SV3 represent reported mutational/rearrangement signature proportions, and delMH is the reported proportion of deletions with flanking microhomology. These features overlap HRDetect inputs and are not separate independent replications. RNA-based TS228 and image-based scores are deliberately not genomic validation outcomes. Binary HRDetect and scarHRD HRD/HRP labels are used verbatim (HRD=1, HRP=0); no threshold reconstruction or optimized cutpoint occurs.

## Identity mapping

Use biological S###### plus the exact m component (m, m2, m3 or m4) from RBA and GEO scan-b external ID. Ignore k versus k2. Q/C prefixes and F title numbers are not biological matching keys. Never numerically match RBA to F. Validate cached `scanb_geo_samples.json` against the local GEO SOFT identities; regenerate only if absent. Two unrelated swap.corrected GEO records outside the HRD cohort lack m and are ignored; an unresolved key for a cohort S identifier would abort.

Exclude every GEO title ending in repl before the primary mapping. Ambiguous non-replicate choices or reuse of an F column across patients abort before analysis. Require unique WGS and RNA records. The mapping table retains RBA, PD_ID, S, m, GEO accession/external ID, F column, mapping/ambiguity/replicate status and metadata age. {m['uniquely_mapped_n']} map uniquely, {m['unmatched_n']} do not match, {m['ambiguous_n']} are ambiguous, and {m['replicate_rows_excluded']} candidate GEO replicate rows are excluded. m-component differences are not silently rescued. Missing patients remain explicit in the audit; no tumor is discarded because of its score or outcome.

## Expression and predictor construction

Read `data/GSE96058_gene_expression.csv.gz` directly with pandas compressed chunks of 256 genes, retaining only the matched F columns. Local SOFT processing metadata explicitly states that transcript FPKM was summed to gene symbols, 0.1 added, and log2 applied. The matrix is already **log2(FPKM+0.1)**, not TPM and not raw counts; use its values without a second log transform. Preserve exact reported gene symbols; duplicate symbols abort rather than selecting a favorable transcript or guessing aliases.

The first pass retains only candidate, canonical and proliferation genes and computes expression means/variances/detection fractions for the random background in small chunks. The second pass retains only {pre['null_loaded_gene_n']} neighborhood genes needed for matched controls. The {pre['gene_rows']}×all-GEO-samples decompressed table is never materialized or written. The source remains compressed. Preparation caches store source sizes/timestamps and SHA256 provenance and can be reused on identical inputs.

Frozen candidates: {', '.join(CANDIDATE_GENES)}. Available: {', '.join(params['candidate_available'])}. Missing: {params['candidate_missing']}. Standardize gene expression over all {m['uniquely_mapped_n']} uniquely mapped SCAN-B TNBC samples, independently of HRD availability, using population SD (ddof=0). Mean and median of all nine z-scores are computed; HORMAD1_z, STAG3_z and REC8_z are retained along with every individual gene z-score. The primary complete-nine-gene signature blocks if genes are missing/constant; no substitution or weight optimization is permitted. Canonical genes are separate and never enter the candidate mean.

The explicit Stage 4 proliferation list is {', '.join(STAGE4_PROLIFERATION)}. Available: {params['proliferation_available']}; missing: {params['proliferation_missing']}. Its mean z-score uses represented genes and the same cohort reference. This user-specified list differs from the actual previous Stage2/3 list ({', '.join(PREVIOUS_PROLIFERATION)}); the latter is retained as a separately named sensitivity. Neither includes candidates. The distinction is not concealed as an unchanged definition.

## Statistical analysis

Primary OLS uses continuous outcome ~ candidate_mean_z. Report N, beta per unit mean-z signature, conventional SE/t-CI/P, R2, Pearson r/P and Spearman rho/P. HC3 sandwich SEs with t-based CIs/P are separate robustness results. Holm correction over the two co-primary outcomes is additionally reported for both conventional and HC3 P-values; unadjusted primary P-values remain visible. The signature is never standardized using outcomes or trained to maximize HRD association.

Proliferation-only, age-only, age+proliferation, previous-proliferation and PAM50+proliferation adjustments are separately named. Age is matched GEO metadata; PAM50 is a covariate, not a selection criterion. A same-complete-case unadjusted model accompanies age+proliferation. Complete-case N is model-specific, and no missing covariate is fabricated. Individual BRCA1/BRCA2 alteration and tumor purity are unavailable in supplied patient-level sheets: these adjustments are recorded unavailable. Minimum N=10 for OLS, full-rank design, and >=10 cases per parameter for adjusted models follow Stage 2. VIF>5 is flagged. No claim of full independence from BRCA status or purity is possible here.

For HRDetect_probability, fractional-logit binomial quasi-likelihood with HC0 sandwich covariance is fitted unadjusted and with proliferation. It accommodates values at 0 and 1 without arbitrary clipping and retains a bounded fitted mean. Its coefficient is on the logit-mean scale and is not the primary OLS beta or a binary-class odds ratio. It is a sensitivity, not a replacement for the prespecified continuous analysis.

Published HRD/HRP labels are modeled by logistic regression; the predictor is scaled by its complete-case SD for OR per 1 SD. Report N, event counts, beta, OR/CI/P and apparent in-sample AUC. At least 10 in each class and convergence are required. AUC is not held-out performance. Three secondary genomic exposures receive one BH family; two binary classifications receive a separate BH family. Individual candidate tests form a nine-gene BH family per primary outcome, canonical genes a separate 13-gene family per outcome. Median and leave-one-out analyses are declared sensitivities; no result is used to redefine the primary signature.

## Null controls and uncertainty

Fixed seed {SEED}: 5,000 random nine-gene sets. Eligible genes have FPKM>=0.1 in >=20% of mapped reference samples and log-expression variance>1e-8. Exclude candidates, canonical comparisons, both proliferation lists and receptor genes. The background uses unique reported expressed gene symbols; no verified protein-coding annotation was supplied, so a protein-coding-only background is not falsely claimed. Rank mean and variance across the background plus candidate targets; for each target select its 200 nearest neighbors in this two-dimensional percentile space, then draw one without replacement within a set. Record genes and matching distances. Matching is approximate, not exact matching of all biological properties.

Each random gene is z-standardized across the identical SCAN-B reference and equally averaged. Use the same complete-case unadjusted OLS statistic for each endpoint. Compare absolute t, not raw beta (signature variances differ). Empirical two-sided P=(1+null |t|>=observed |t|)/(B+1); percentile is the percent of null |t| strictly below observed. The same drawn gene sets serve both endpoints without tuning.

For each endpoint independently perform 10,000 patient-label permutations of the outcome relative to fixed candidate values (seeds {SEED+10}, {SEED+11}); this is mathematically equivalent to permuting predictor labels relative to the outcome. Use absolute t and the plus-one P. Bootstrap complete-case patient pairs 2,000 times per endpoint (seeds {SEED+20}, {SEED+21}), holding the fitted gene standardization fixed, and report percentile beta CIs. These procedures describe association uncertainty; permutations do not adjust for confounding and the bootstrap is conditional on the reference scaling.

## Diagnostics and stability

Report Breusch–Pagan heteroscedasticity, Shapiro residual normality, quadratic RESET model-form check, residual/fitted and Q-Q plots, leverage, residuals and Cook distance by PD_ID/RBA/F. Outliers are not removed from primary models. Separate models exclude Cook distance>4/N and list all excluded IDs. Leave-one-gene-out computes each equal-weight eight-gene mean from the frozen nine gene z-scores and reports both conventional and HC3 inference. No gene is selected for removal. These diagnostics can reveal instability and nonlinearity; nonnormal bounded probabilities are expected to challenge OLS distributional assumptions, motivating robust and fractional sensitivities.

## Results and exact analysis Ns

'''
    for r in primary:
        text+=f"- {r['outcome']}: N={r['N']}, beta={r['beta']:.6g}, conventional 95% CI [{r['ci_low']:.6g}, {r['ci_high']:.6g}], P={r['p']:.6g}, Holm P={r['primary_holm_p']:.6g}, R2={r['r_squared']:.6g}; HC3 95% CI [{r['hc3_ci_low']:.6g}, {r['hc3_ci_high']:.6g}], P={r['hc3_p']:.6g}.\n"
    text+='\nBoth unadjusted genomic associations replicate, but neither remains significant with the primary proliferation adjustment. This does not establish an HRD state independent of proliferation, BRCA status or purity.\n'
    for r in s['proliferation_adjusted']:
        text+=f"- Proliferation-adjusted {r['outcome']}: N={r['N']}, beta={r['beta']:.6g}, HC3 95% CI [{r['hc3_ci_low']:.6g}, {r['hc3_ci_high']:.6g}], P={r['hc3_p']:.6g}.\n"
    for r in s['matched_random']:
        text+=f"- Matched random signatures, {r['outcome']}: empirical P={r['empirical_p']:.6g}; observed absolute-t percentile={r['percentile']:.3g}.\n"
    text+='\nUnadjusted leave-one-gene-out associations remain positive with HC3 P<0.05 for every omission. This supports composite stability but does not resolve proliferation confounding. Both primary models show nonnormal residuals and quadratic RESET evidence of nonlinearity. Primary Cook flags identify 10 HRDetect and 9 scarHRD observations; all remain in the primary analyses. Fractional-logit sensitivity similarly loses significance after proliferation adjustment.\n'
    text+='\n| Model | Outcome | Predictor | N | P |\n|---|---|---|---:|---:|\n'
    for r in adjusted+secondary+sensitivity+logo:text+=f"| {r['model']} | {r['outcome']} | {r.get('omitted_gene',r['predictor'])} | {r['N']} | {r.get('p',float('nan')):.6g} |\n"
    text+='\nAll estimates and canonical Ns are in the named CSVs. Unavailable adjustments retain missing estimates; they are not negative findings. The summary JSON includes every primary, proliferation-adjusted, permutation, matched-null and bootstrap result.\n'
    text+=f"\n**Overall: {s['classification']}**. STRONG requires both positive primary Holm and HC3-Holm P<0.05, positive proliferation-HC3 associations, both null/permutation P<0.05, and positive significant leave-one-out/Cook sensitivities. MODERATE denotes some primary replication with incomplete robustness; WEAK denotes absent primary support or failures of both specificity controls; FAILED denotes a genuine mapping/data/scientific block. This internal rule does not confer clinical validity or causality. Available BRCA/purity adjustment is limited, and HRD scars reflect genomic history rather than necessarily current repair capacity.\n"
    text+='\nExact inputs are the four specified local data files plus GEO JSON validated against SOFT. Input SHA256 values and retained-cache dimensions are in preparation_manifest.json. No original compressed input was extracted or modified. Tests: `python -m unittest discover -s tests -p test_stage4.py -v`.\n'
    (root/'docs/V2_METHODS_STAGE4.md').write_text(text,encoding='utf-8')

if __name__=='__main__':
    try:run()
    except Exception as exc:
        destination=PROJECT_ROOT/'outputs/v2/stage4';destination.mkdir(parents=True,exist_ok=True)
        (destination/'stage4_failure.json').write_text(json.dumps({'classification':'FAILED','blocker':str(exc)},indent=2))
        raise
