import json
from importlib.metadata import version
import numpy as np
import pandas as pd
from .config import PROJECT_ROOT,CANDIDATE_GENES
from .stage2 import clean
from .utils import sha256
from .stage3_patient import metabric,gse25066
from .stage3_single_cell import analyze as single_cell
from .stage3_functional import functional_results
from .figures_stage3 import make_figures

def evidence_row(gene,domain,endpoint,result=None,note=''):
    measure='linear regression coefficient per predictor unit'
    if domain=='single_cell':measure='paired mean difference in program z-score' if gene=='composite' else 'paired mean difference in log2(CPM+1)'
    elif endpoint in ['OS','RFS','DRFS']:measure='log hazard ratio per predictor SD'
    row=dict(gene=gene,domain=domain,endpoint=endpoint,N=np.nan,effect=np.nan,effect_measure=measure,ci_low=np.nan,ci_high=np.nan,p=np.nan,fdr=np.nan,support='unavailable',note=note)
    if result is None:return row
    for k in ['N','ci_low','ci_high','p','fdr']:
        if k in result:row[k]=result[k]
    row['effect']=result.get('effect',result.get('beta',np.nan))
    if result.get('status')=='estimated':
        significance=row['fdr'] if pd.notna(row['fdr']) else row['p']
        direction_ok=True if endpoint in ['OS','RFS','DRFS'] else row['effect']>0
        row['support']='supports' if pd.notna(significance) and significance<.05 and direction_ok else 'does_not_support'
    else:row['note']=result.get('reason',note)
    return row

def integrate(root,meta,geo,sc,dep):
    primary=pd.read_csv(root/'outputs/v2/stage2/primary_hrd_model.csv').iloc[0].to_dict()
    genes=pd.read_csv(root/'outputs/v2/stage2/single_gene_results.csv')
    secondary=pd.read_csv(root/'outputs/v2/stage2/secondary_outcomes.csv');rows=[]
    for g in ['composite',*CANDIDATE_GENES]:
        r=primary if g=='composite' else next((x for x in genes.to_dict('records') if x['gene']==g),None)
        rows.append(evidence_row(g,'TCGA_HRD','HRD_total',r,'Frozen TCGA discovery; not external replication'))
        for outcome in secondary.outcome:
            r=secondary[secondary.outcome==outcome].iloc[0].to_dict() if g=='composite' else None
            rows.append(evidence_row(g,'TCGA_CIN',outcome,r,'Individual-gene CIN tests not run' if g!='composite' else 'Frozen TCGA result'))
        for endpoint in ['altered_gene_fraction','altered_gene_fraction_adjusted','OS','RFS','HRD_scar']:
            r=next((x for x in meta if x['gene']==g and x['endpoint']==endpoint),None)
            rows.append(evidence_row(g,'METABRIC',endpoint,r,'Gene-weighted CNA proxy is not HRD scar'))
        for endpoint in ['pCR','DRFS']:
            r=next((x for x in geo if x['endpoint']==endpoint),None) if g=='composite' else None
            rows.append(evidence_row(g,'GSE25066',endpoint,r,'GPL96 mapping unavailable; no gene coverage inferred'))
        r=next((x for x in sc if x['gene']==g and x['population']=='TNBC' and x['endpoint']=='malignant_vs_nonmalignant_combined'),None)
        rows.append(evidence_row(g,'single_cell','malignant_vs_nonmalignant_combined',r,'Paired donor/sample analysis; RNA specificity is not HRD'))
        if g=='composite' and dep.get('status')=='estimated':
            strongest=dep.get('strongest_dependency') or {}
            rows.append(dict(gene=g,domain='DepMap',endpoint='genome_wide_CRISPR',N=dep.get('analysis_model_n',np.nan),effect=strongest.get('beta',np.nan),effect_measure='strongest genome-wide CRISPR gene-effect beta; negative means stronger dependency',ci_low=np.nan,ci_high=np.nan,p=strongest.get('p',np.nan),fdr=strongest.get('fdr',np.nan),support='supports' if dep.get('FDR_significant_dependencies',0)>0 else 'does_not_support',note='Genome-wide screen; strongest dependency='+str(strongest.get('gene',''))+'; cohort='+str(dep.get('analysis_cohort',''))))
        else:
            rows.append(evidence_row(g,'DepMap','genome_wide_CRISPR',note=dep.get('reason','Composite-only DepMap screen; individual candidate genes not tested') if g=='composite' else 'Composite-only DepMap screen; individual candidate genes not tested'))
        rows.append(evidence_row(g,'drug','compound_wide_response',note='Drug response resources missing; not tested'))
    return pd.DataFrame(rows)

def run():
    root=PROJECT_ROOT;dest=root/'outputs/v2/stage3';dest.mkdir(parents=True,exist_ok=True)
    frozen=['outputs/v2/cohort/tcga_tnbc_patients.csv','outputs/v2/predictors/candidate_signature.csv','outputs/v2/stage2/primary_hrd_model.csv',
            'outputs/v2/stage2/adjusted_models.csv','outputs/v2/stage2/single_gene_results.csv','outputs/v2/stage2/secondary_outcomes.csv']
    before={p:sha256(root/p) for p in frozen}
    protocol=dict(frozen_candidate_genes=CANDIDATE_GENES,metabric_duplicate_gene_rule='mean all reported symbol rows',
                  metabric_clinical_rule='ER IHC negative; no positive reported ER conflict; PR negative; HER2 negative; primary tumors',
                  metabric_cna_rule='fraction of unique assayed DNA gene calls with absolute state>=1',
                  metabric_primary_family=['CNA proxy','OS','RFS'],gse_platform='GPL96 only; no mapping inference',
                  sc_minimum_cells=20,sc_primary='TNBC paired malignant vs pooled nonmalignant',
                  sc_tests='two-sided donor-level Wilcoxon; bootstrap mean paired differences, 2000 resamples',
                  no_large_downloads=True,no_gene_weight_optimization=True)
    (dest/'protocol.json').write_text(json.dumps(clean(protocol),indent=2))
    meta,mdata,ms=metabric(root);geo,gdata,gs=gse25066(root);sc,sdata,ss=single_cell(root)
    dep,drug,resources=functional_results(root)
    evidence=integrate(root,meta,geo,sc,dep);evidence.to_csv(dest/'integrated_evidence.csv',index=False)
    make_figures(root,mdata,meta,gdata,gs,sdata,ss,evidence)
    adjusted=next(r for r in meta if r['endpoint']=='altered_gene_fraction_adjusted')
    # Proliferation-adjusted replication
    patient_support=ms['main_result'].get('fdr',1)<.05 and adjusted.get('hc3_p',1)<.05
    geo_support=any(r['status']=='estimated' and r.get('fdr',1)<.05 for r in geo)
    functional_support=dep.get('status')=='estimated' and dep.get('FDR_significant_dependencies',0)>0
    external_support=patient_support or geo_support
    if external_support and functional_support:
        classification='STRONG'
    elif external_support or functional_support:
        classification='PROMISING'
    else:
        classification='WEAK'
    basis=[]
    basis.append('external patient support present' if external_support else 'no robust external patient support after prespecified controls')
    basis.append('DepMap functional support present' if functional_support else ('DepMap screen negative at FDR<0.05' if dep.get('status')=='estimated' else 'DepMap unavailable'))
    summary=dict(classification=classification,classification_basis='; '.join(basis),
                 metabric=ms,metabric_adjusted=adjusted,gse25066=gs,single_cell=ss,
                 depmap=dep,
                 drugs=drug,
                 missing_resources=resources[resources.status=='missing'].to_dict('records'),
                 direct_METABRIC_HRD='Not identified locally; exact resource/size unknown; not approximated',
                 stage1_stage2_source_hashes=before,packages={n:version(n) for n in ['numpy','pandas','scipy','statsmodels','matplotlib']},
                 functional_analysis_implementation='Genome-wide DepMap CRISPR gene-effect screen implemented when local resources are present; drug analysis remains gated until separately implemented',
                 no_new_downloads=True)
    (dest/'stage3_summary.json').write_text(json.dumps(clean(summary),indent=2,allow_nan=False))
    write_docs(root,summary,meta,sc,resources)
    for path,digest in before.items():
        if sha256(root/path)!=digest:raise RuntimeError('Frozen input changed: '+path)
    print(json.dumps(clean({'METABRIC_N':ms['clinical_tnbc_n'],'METABRIC_CNA':ms['main_result'],
                            'GSE25066':gs,'single_cell':ss['primary_result'],'DepMap':dep,'classification':classification}),indent=2))
    return summary

def write_docs(root,s,meta,sc,resources):
    m=s['metabric'];g=s['gse25066'];single=s['single_cell'];r=m['main_result'];a=s['metabric_adjusted'];sr=single['primary_result']
    methods=f'''# V2 Stage 3 methods: external and functional validation

Run `python -m src.stage3` from V2. Stage 1 and Stage 2 results are read, not rebuilt. Source hashes before/after are checked. Original notebooks and datasets are read-only. No new datasets were downloaded. Candidate genes remain {', '.join(CANDIDATE_GENES)}; no weights, genes or cutpoints were optimized against validation outcomes.

## METABRIC

Sources are the existing `dataset/brca_metabric_clinical_data.tsv` and `results/metabric_extracted/brca_metabric/data_mrna_illumina_microarray.txt`, `meta_mrna_illumina_microarray.txt`, `data_cna.txt` under the frozen original root. Source metadata explicitly identifies log2 Illumina HT-12 v3 intensities; no second log transform is applied. Only needed TNBC expression/CNA columns are loaded.

Require primary tumors, Negative ER measured by IHC, Negative PR Status and Negative HER2 Status, and no conflicting positive reported ER Status. The source typo `Positve` is normalized to positive. The generic reported three-negative fields yield {m['source_reported_three_negative_n']} records, but {m['reported_negative_but_IHC_positive_n']} have positive ER IHC and 11 lack resolved ER IHC. These are excluded before modeling, yielding {m['clinical_tnbc_n']} unique patients. IHC is used rather than the discordant generic ER field; PR/HER2 are the source-reported receptor annotations, whose individual assay details are not available in this table. PAM50 basal and the SNP6 HER2 field do not define TNBC. A complete receptor decision ledger is retained. No expression medians or inferred receptor rescue is used.

Retain all reported expression rows for each gene symbol and average log intensities equally for duplicates, including both STAG3 rows. This is a fixed probe/row aggregation rule; neither row is selected by expression variance or association. Z-standardize each of the nine genes using population SD within the {m['rna_n']} METABRIC clinical-TNBC RNA samples, then take their equal-weight mean. All nine are required. This is cross-platform cohort restandardization, not application of a TCGA absolute assay cutoff.

Independent DNA endpoint: proportion of nonmissing unique assayed genes with absolute discrete CNA state >=1, with median aggregation of duplicate gene DNA rows. This is a **gene-weighted CNA burden proxy**, not length-weighted FGA, allele-specific HRD scars, CIN rate or direct HRD replication. It uses genome-wide DNA measurements and never candidate RNA. Primary validation association uses continuous signature OLS with conventional and HC3 CIs/P; an age+published 11-gene proliferation adjustment is a prespecified sensitivity. The proliferation implementation reuses the disjoint Stage 2 gene set.

OS and RFS are secondary clinical endpoints. Source months and explicit event-status labels are used, with positive follow-up and known 0/1 status required. Cox proportional hazards models use Efron ties and continuous signature per one reference-cohort SD, with HR, 95% CI, P, N and event counts. No optimized cutpoints or Kaplan–Meier high/low groups are created. A Schoenfeld-residual/log-time Spearman diagnostic is reported, not a formal multivariate PH test. The composite CNA/OS/RFS family receives BH correction. Individual-gene CNA and survival tests are exploratory, with separate nine-gene BH families per endpoint; they do not alter the composite.

## GSE25066

Read embedded metadata in the local `GSE25066_series_matrix.txt.gz`. Every sample declares GPL96. Clinical TNBC requires `er_status_ihc=N`, `pr_status_ihc=N`, `her2_status=N`; the fields using ESR1 to resolve indeterminate ER are explicitly unused. Preserve GSM and embedded sample identity. pCR is observed `pathologic_response_pcr_rd`, not any expression-derived response prediction; DRFS uses the source event flag and `drfs_even_time_years` field.

The authoritative GPL96 annotation is absent. Therefore gene coverage cannot be verified and neither reduced signature nor pCR/DRFS association was estimated. Missing annotation is not zero gene expression or evidence that all nine genes are absent. No GPL570 mapping was substituted. The implemented adapter, if the configured GPL96 annotation is later supplied, retains only unambiguous candidate probe mappings, averages all probes for each represented candidate, z-standardizes within clinical TNBC and uses all represented frozen genes without outcome-based selection. Logistic pCR and Cox DRFS remain continuous-predictor tests with a two-endpoint BH family. That branch is not claimed as executed here.

## GSE176078

Use the local `GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz` and existing matching `metadata.csv`, `count_matrix_genes.tsv`, `count_matrix_barcodes.tsv` in `results/sc_extracted/Wu_etal_2021_BRCA_scRNASeq`. The [original study](https://www.nature.com/articles/s41588-021-00911-1) supplies the Cancer Epithelial and other cell annotations; they are not inferred from the candidate genes. The sparse archive contains 29,733 genes, 100,064 cells and 177,994,136 nonzero entries. It is streamed in 200,000-triplet chunks, accumulating library counts and the nine candidate counts by exact `orig.ident` × class. No dense gene-by-cell matrix or original-directory extraction is created. Barcode prefixes and metadata identities are checked; original CID suffixes are retained rather than guessed away.

There are {single['all_donor_sample_n']} published donor/sample units, including {single['tnbc_donor_sample_n']} labeled TNBC by the source metadata. We treat exact orig.ident as the donor/sample unit; an additional independent patient-identity crosswalk is not available, so we do not claim identity resolution beyond the published grouping. Cell classes: Cancer Epithelial→malignant epithelial; Normal Epithelial→normal epithelial; T/B/Myeloid/Plasmablast→immune; CAF/PVL/Endothelial→stromal. Single-cell TNBC labels are source-provided, not reconstructed from RNA receptor thresholds. No independent HRD annotation is available.

Require >=20 cells and positive library counts per donor/class. Sum integer counts and normalize to counts per million of all gene counts, then log2(CPM+1). This is pseudobulk abundance, not a mean of log-normalized cells. The gene z-score reference consists of eligible TNBC donor-by-four-disjoint-class pseudobulks; fit all nine genes and equally average their z-scores. Apply those same reference parameters to the independently pooled nonmalignant compartment and all-breast sensitivity. True observed zero counts remain zeros before transformation; genes with no reference variance would block a complete score rather than be removed.

Primary contrast is within-donor malignant minus combined nonmalignant program activity in TNBC: {sr['N']} eligible paired units. Use two-sided paired Wilcoxon and a 2,000-resample paired-donor bootstrap CI for the mean difference, fixed seed 20260908. Compare malignant to each separate class as descriptive secondary contrasts, with BH across available four composite contrasts per population. For each of nine genes, compare paired log2(CPM+1) abundance, with BH across nine genes; these identify expression contributors, not reweighted signatures. A minimum of five pairs is required for an inferential test. All-breast comparisons are separately labeled sensitivities. Cells are never independent replicates.

## DepMap and drugs

When the local DepMap model, expression and CRISPR gene-effect files are present, Stage 3 performs a genome-wide dependency screen. Breast lineage and TNBC eligibility are derived only from curated non-expression `Model.csv` metadata. A TNBC model must be explicitly annotated with `TNBC` or `triple-negative`; ESR1, PGR and ERBB2 expression are never used for eligibility. If fewer than eight curated TNBC models are available, the code performs a clearly labelled breast-lineage exploratory screen instead of inventing TNBC labels.

All nine frozen candidate genes are z-standardized within the analysis model set and equally averaged. A disjoint 11-gene proliferation score is used as a covariate when sufficiently complete. For each CRISPR dependency, gene effect is regressed on the continuous candidate signature, both unadjusted and proliferation-adjusted, followed by BH FDR across the genome. Negative beta means stronger dependency with increasing candidate-signature activity. DDR/replication-stress genes are annotated only after the genome-wide test. The current implementation does not reweight the candidate signature and does not use known DDR genes to prefilter the screen.

Drug response remains unavailable until a separately reviewed PRISM/GDSC model and compound crosswalk is implemented. No therapeutic efficacy inference is made.

## Integration and figures

`integrated_evidence.csv` is a long table for the composite and all nine genes. It preserves actual effect, CI, P/FDR, N and whether the prespecified evidence supports, does not support at threshold, or is unavailable. Positive genomic/malignant-specificity effects with q<0.05 (or primary P<0.05) count as directional support. Survival entries describe two-sided association support without asserting a favorable direction. Untested gene-by-domain combinations remain unavailable. No points are assigned or summed. The heatmap encodes categories only; numeric colors are not an evidence score. TCGA entries read the frozen Stage 2 tables, including their actual N=103 primary model, and are not external replication.

Figures are PDF plus 300-dpi PNG. METABRIC shows the CNA proxy and continuous-score survival HRs; GSE25066 shows observed response counts only, with the mapping blocker explicit; single-cell shows paired donor values. DepMap volcano/forest figures are explicit unavailability notices, not empty discovery claims. The evidence heatmap includes both unadjusted and adjusted METABRIC CNA results.

Current classification uses the frozen independence question conservatively: unadjusted CNA association alone that loses support after age/proliferation control is not robust external replication. STRONG requires an independent patient result plus coherent orthogonal functional support; PROMISING requires supporting external or functional evidence with incomplete replication; WEAK reflects insufficient meaningful replication in the available components. Unavailable functional data are not negative experiments. No signature alteration follows the label.
'''
    (root/'docs/V2_METHODS_STAGE3.md').write_text(methods,encoding='utf-8')
    osr=next(x for x in meta if x['endpoint']=='OS' and x['gene']=='composite');rfs=next(x for x in meta if x['endpoint']=='RFS' and x['gene']=='composite')
    dependency=resources[resources.resource_id.isin(['depmap_effect','depmap_expression','depmap_model','depmap_default','depmap_profiles','depmap_condition','depmap_readme'])]
    res=f'''# V2 Stage 3 results

**Classification: {s['classification']}.** The unadjusted METABRIC CNA proxy association is positive but is not supported after age/proliferation adjustment; the TNBC single-cell composite is not enriched in malignant cells; functional validation is unavailable. This does not refute the TCGA HRD result or constitute negative DepMap/drug experiments.

## METABRIC

Clinical TNBC and RNA N={m['clinical_tnbc_n']}; all nine genes available. The 320 generic reported triple-negative records were reduced to 258 by requiring resolved negative ER IHC and excluding the 51 positive-IHC conflicts and 11 missing IHC records. Gene-row aggregation retained both STAG3 rows.

Independent gene-weighted CNA burden: beta={r['beta']:.6g}, 95% CI [{r['ci_low']:.6g}, {r['ci_high']:.6g}], P={r['p']:.6g}, q={r['fdr']:.6g}, R2={r['r_squared']:.6g}. HC3 P={r['hc3_p']:.6g}. This is partial genomic replication, not direct HRD replication.

Age/proliferation-adjusted beta={a['beta']:.6g}, 95% CI [{a['ci_low']:.6g}, {a['ci_high']:.6g}], P={a['p']:.6g}; HC3 P={a['hc3_p']:.6g}. This adjusted result does not support an independent association at P<0.05.

OS: N={osr['N']}, events={osr['events']}, HR per signature SD={osr['hazard_ratio']:.6g}, 95% CI [{osr['hr_ci_low']:.6g}, {osr['hr_ci_high']:.6g}], P={osr['p']:.6g}, q={osr['fdr']:.6g}. RFS: N={rfs['N']}, events={rfs['events']}, HR={rfs['hazard_ratio']:.6g}, 95% CI [{rfs['hr_ci_low']:.6g}, {rfs['hr_ci_high']:.6g}], P={rfs['p']:.6g}, q={rfs['fdr']:.6g}. Neither survival endpoint supports an association.

## GSE25066

Clinical TNBC N={g['clinical_tnbc_n']}; observed pCR status N={g['pcr_observed_n']}, including {g['pcr_events']} pCR events. All samples declare GPL96. The 4,522,748-byte configured GPL96 annotation is missing, so candidate platform coverage and pCR/DRFS signature associations are unavailable. No missing gene was assigned zero and no GPL570 map was used.

## Single-cell

{single['all_donor_sample_n']} donor/sample units overall, {single['tnbc_donor_sample_n']} source-labeled TNBC; {sr['N']} TNBC units have >=20 malignant and nonmalignant cells. Mean paired malignant-minus-nonmalignant composite difference={sr['effect']:.6g}, bootstrap 95% CI [{sr['ci_low']:.6g}, {sr['ci_high']:.6g}], two-sided Wilcoxon P={sr['p']:.6g}, q={sr['fdr']:.6g}. The frozen composite is **not supported as malignant-specific** in this analysis.

HORMAD1, SMC1B and SYCP2 show positive descriptive malignant-cell differences, while STAG3 and REC8 are lower in malignant than pooled nonmalignant pseudobulk. No individual candidate survives the nine-gene TNBC FDR family. All gene-level estimates and all-breast sensitivities are retained in `single_cell_results.csv`; these findings do not justify removing genes or changing weights. Low expression, only eight paired TNBC units and source-defined donor/sample grouping limit precision. No HRD was inferred from single-cell RNA.

## Functional validation and missing resources

DepMap status: **{s['depmap'].get('status')}**. Breast models={s['depmap'].get('breast_model_n')}; curated TNBC models={s['depmap'].get('tnbc_model_n')}; models in the functional screen={s['depmap'].get('analysis_model_n')}. Analysis cohort={s['depmap'].get('analysis_cohort')}. FDR-significant stronger dependencies={s['depmap'].get('FDR_significant_dependencies')}. Strongest dependency={s['depmap'].get('strongest_dependency')}. TNBC eligibility used curated non-expression metadata only; ESR1/PGR/ERBB2 expression thresholds were not used. Drug response status is **{s['drugs'].get('status')}** and no drug efficacy claim is made. Direct METABRIC HRD scars remain unavailable.

| Resource | Provider file ID | Bytes | Filename |
|---|---|---:|---|
'''
    for row in resources[resources.status=='missing'].to_dict('records'):
        res+=f"| {row['resource_id']} | {row['file_id']} | {row['size_bytes']:,} | {row['filename']} |\n"
    res+='\nNo new datasets were downloaded. Original data and frozen Stage 1/2 inputs were unchanged. Missing resources block only their components; results above use the available patient and single-cell data.\n'
    (root/'docs/V2_RESULTS_STAGE3.md').write_text(res,encoding='utf-8')

if __name__=='__main__':run()
