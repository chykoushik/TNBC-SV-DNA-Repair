import csv
import datetime
import gzip
import hashlib
import json
from collections import Counter, defaultdict
import numpy as np
from . import config
from .utils import open_text, output_path, sha256, tcga_ids, write_csv, write_json
from .signatures import fit_standardization, score_predictors, PREDICTOR_COLUMNS

PREFIX = 'outputs/v2/expression/'

def parse_star(path):
    genes, tpm, counts, seen, counters, comments = [], [], [], set(), {}, []
    with open_text(path) as f:
        line = f.readline()
        while line.startswith('#'):
            comments.append(line.strip()); line = f.readline()
        reader = csv.DictReader(f, fieldnames=line.rstrip('\r\n').split('\t'), delimiter='\t')
        required = {'gene_id','gene_name','gene_type','unstranded','tpm_unstranded'}
        if not required.issubset(reader.fieldnames): raise ValueError('Missing augmented STAR columns')
        for r in reader:
            if None in r or any(v is None for v in r.values()): raise ValueError('Malformed STAR row')
            gid = r['gene_id']
            if gid.startswith('N_'):
                counters[gid] = int(r['unstranded']); continue
            if not gid.startswith('ENSG') or gid in seen: raise ValueError('Invalid/duplicate Ensembl gene ID')
            seen.add(gid)
            v, c = float(r['tpm_unstranded']), float(r['unstranded'])
            if not np.isfinite(v) or not np.isfinite(c) or min(v,c) < 0 or c != int(c):
                raise ValueError('Invalid TPM/count')
            genes.append((gid,r['gene_name'],r['gene_type']));tpm.append(v);counts.append(c)
    if not genes: raise ValueError('No gene rows')
    return genes, np.asarray(tpm), np.asarray(counts), counters, comments

def resolve_metadata(record, clinical, detail=None):
    cases = record['cases']
    if len(cases) != 1 or cases[0]['submitter_id'] != clinical['patient_id']: raise ValueError('Case ID conflict')
    samples = cases[0]['samples']
    if len(samples) != 1: raise ValueError('Ambiguous sample metadata')
    sample = samples[0]['submitter_id']; ids = tcga_ids(sample); frozen = tcga_ids(clinical['sample_id'])
    if samples[0]['sample_type'] != 'Primary Tumor' or ids['sample_type_code'] != '01': raise ValueError('Not primary tumor')
    if ids['sample_key'] != frozen['sample_key']: raise ValueError('Frozen clinical sample key mismatch')
    if frozen['vial'] and frozen['normalized_sample_id'] != ids['normalized_sample_id']: raise ValueError('Vial mismatch')
    aliquots = set()
    if detail:
        d = detail['data']
        if d['file_id'] != record['file_id'] or d['md5sum'] != record['md5sum']: raise ValueError('GDC detail snapshot conflict')
        for case in d['cases']:
            if case['submitter_id'] != clinical['patient_id']: raise ValueError('Detailed case conflict')
            for s in case['samples']:
                if s['submitter_id'] != sample: raise ValueError('Detailed sample conflict')
                for portion in s.get('portions',[]):
                    for analyte in portion.get('analytes',[]):
                        for a in analyte.get('aliquots',[]): aliquots.add(a['submitter_id'])
    if len(aliquots) > 1: raise ValueError('Multiple aliquots linked to one file')
    aliquot = next(iter(aliquots), '')
    if aliquot and not aliquot.startswith(sample+'-'): raise ValueError('Aliquot/sample conflict')
    return sample, aliquot

def select_replicate(records):
    return sorted(records, key=lambda x:(-x['assigned_reads'],x['gdc_sample_id'],x['aliquot_id'],x['file_id']))[0]

def run():
    root = config.PROJECT_ROOT
    frozen_path = root/'outputs/v2/cohort/tcga_tnbc_patients.csv'
    frozen_hash = sha256(frozen_path)
    with frozen_path.open(newline='',encoding='utf-8') as f: cohort = list(csv.DictReader(f))
    if len({r['patient_id'] for r in cohort}) != len(cohort): raise ValueError('Frozen cohort has duplicate patients')
    clinical = {r['patient_id']:r for r in cohort}
    inventory_path = root/'outputs/acquisition/tcga_rna_files.json'
    inventory = json.loads(inventory_path.read_text())
    files = {r['file_id']:r for r in inventory['files']}
    detail_path = root/(PREFIX+'gdc_aliquot_metadata.json')
    details = json.loads(detail_path.read_text()) if detail_path.exists() else {}
    # Frozen file inventory
    local_paths = sorted(config.STAR_DIR.rglob('*.tsv'))
    manifest, parsed, groups = [], {}, defaultdict(list)
    for path in local_paths:
        fid = path.parent.name; meta = files.get(fid)
        if not meta: raise ValueError('Local STAR UUID absent from pinned GDC inventory: '+fid)
        patients = [c['submitter_id'] for c in meta['cases']]
        if not any(p in clinical for p in patients): continue
        row = dict(file_id=fid,patient_id=patients[0],clinical_sample_id=clinical[patients[0]]['sample_id'],
                   gdc_sample_id='',aliquot_id='',expression_sample_id='',source_path=str(path),source_md5=meta['md5sum'],
                   sha256='',assigned_reads=0,selected=False,status='',selection_rule='highest assigned reads; lexical sample/aliquot/file tie-break')
        try:
            if meta.get('analysis',{}).get('workflow_type') != 'STAR - Counts': raise ValueError('Not a GDC STAR - Counts file')
            if path.name != meta['file_name'] or path.stat().st_size != meta['file_size']: raise ValueError('Filename/size differs from GDC inventory')
            md5 = hashlib.md5()
            with path.open('rb') as f:
                for chunk in iter(lambda:f.read(1048576),b''): md5.update(chunk)
            if md5.hexdigest() != meta['md5sum']: raise ValueError('Provider MD5 mismatch')
            row['sha256'] = sha256(path)
            row['gdc_sample_id'],row['aliquot_id'] = resolve_metadata(meta,clinical[row['patient_id']],details.get(fid))
            # UUID identity preservation
            row['expression_sample_id'] = row['aliquot_id'] or row['gdc_sample_id']+'|'+fid
            parsed[fid] = parse_star(path)
            row['assigned_reads'] = int(parsed[fid][2].sum())
            if row['assigned_reads'] <= 0: raise ValueError('Empty assigned-read library')
            row['status'] = 'eligible' if row['aliquot_id'] else 'eligible_aliquot_barcode_unresolved'
            groups[row['patient_id']].append(row)
        except (ValueError, OSError) as exc: row['status']='excluded: '+str(exc)
        manifest.append(row)
    selected = []
    for patient in sorted(groups):
        chosen=select_replicate(groups[patient]);chosen['selected']=True;selected.append(chosen)
        for r in groups[patient]:
            if r is not chosen: r['status']='unselected_duplicate'
    if not selected: raise ValueError('No valid frozen-cohort STAR samples; no predictor outputs replaced')
    genes = parsed[selected[0]['file_id']][0]
    comments = parsed[selected[0]['file_id']][4]
    for r in selected:
        if parsed[r['file_id']][0] != genes or parsed[r['file_id']][4] != comments:
            raise ValueError('Inconsistent gene annotation/order; explicit reconciliation required')
    tpm = np.column_stack([parsed[r['file_id']][1] for r in selected]); logged=np.log2(tpm+1)
    ids=[r['expression_sample_id'] for r in selected]
    matrix_path=output_path(PREFIX+'tcga_tnbc_expression.tsv.gz')
    with gzip.open(matrix_path,'wt',encoding='utf-8',newline='') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(['ensembl_gene_id','gene_symbol','gene_type']+ids)
        for gene,values in zip(genes,logged): writer.writerow(list(gene)+[format(v,'.12g') for v in values])
    symbols=defaultdict(list)
    for i,g in enumerate(genes):symbols[g[1]].append(i)
    expression={}
    requested=config.CANDIDATE_GENES+config.COMPARISON_GENES
    for symbol in requested:
        if len(symbols[symbol])==1: expression[symbol]=dict(zip(ids,logged[symbols[symbol][0]].tolist()))
    parameters=fit_standardization(expression,ids)
    scored={r['expression_sample_id']:r for r in score_predictors(expression,ids,parameters)}
    selected_patients={r['patient_id']:r for r in selected};predictions=[];missing=[]
    for c in cohort:
        p=dict(patient_id=c['patient_id'],sample_id=c['sample_id'],expression_sample_id='',predictor_role='expression_predictor',
               n_genes_scored=0,predictor_status='no_complete_valid_local_STAR',**{k:None for k in PREDICTOR_COLUMNS})
        if c['patient_id'] in selected_patients:p.update(scored[selected_patients[c['patient_id']]['expression_sample_id']])
        else:missing.append(dict(patient_id=c['patient_id'],clinical_sample_id=c['sample_id'],reason='no complete valid matching local STAR file'))
        predictions.append(p)
    gene_qc=[]
    for i,(gid,symbol,kind) in enumerate(genes):
        values=tpm[i];lv=logged[i]
        gene_qc.append(dict(ensembl_gene_id=gid,gene_symbol=symbol,gene_type=kind,
            gene_set='candidate' if symbol in config.CANDIDATE_GENES else 'comparison' if symbol in config.COMPARISON_GENES else 'other',
            symbol_mapping_count=len(symbols[symbol]),n_samples=len(ids),zero_n=int((values==0).sum()),
            below_0_1_tpm_n=int((values<.1).sum()),all_zero=bool((values==0).all()),
            near_zero_expression=bool((values<.1).mean()>=.9),near_zero_log_variance=bool(lv.std()<=1e-8),
            mean_tpm=float(values.mean()),median_tpm=float(np.median(values)),max_tpm=float(values.max()),
            mean_log2_tpm1=float(lv.mean()),sd_log2_tpm1=float(lv.std()),
            p05_log2_tpm1=float(np.quantile(lv,.05)),p95_log2_tpm1=float(np.quantile(lv,.95))))
    mask=np.array([g[2]=='protein_coding' for g in genes])&(tpm.max(axis=1)>=1)&(logged.std(axis=1)>1e-8)
    if mask.sum()<2:raise ValueError('Insufficient variable protein-coding genes for correlation QC')
    correlations=np.corrcoef(logged[mask].T) if len(ids)>1 else np.ones((1,1))
    medians=np.array([np.median(np.delete(correlations[i],i)) if len(ids)>1 else np.nan for i in range(len(ids))])
    center=np.nanmedian(medians);mad=np.nanmedian(abs(medians-center));sample_qc=[]
    for i,r in enumerate(selected):
        robust=float((medians[i]-center)/(1.4826*mad)) if len(ids)>=5 and mad>0 else None
        counters=parsed[r['file_id']][3];assigned=r['assigned_reads'];total=assigned+sum(counters.values())
        sample_qc.append(dict(patient_id=r['patient_id'],clinical_sample_id=r['clinical_sample_id'],gdc_sample_id=r['gdc_sample_id'],aliquot_id=r['aliquot_id'],
            expression_sample_id=ids[i],file_id=r['file_id'],n_gene_rows=len(genes),assigned_reads=assigned,
            assigned_fraction_of_STAR_accounting=assigned/total if total else None,
            detected_genes=int((tpm[:,i]>0).sum()),genes_ge_1_tpm=int((tpm[:,i]>=1).sum()),tpm_sum=float(tpm[:,i].sum()),
            mean_log2_tpm1=float(logged[:,i].mean()),median_log2_tpm1=float(np.median(logged[:,i])),
            p05_log2_tpm1=float(np.quantile(logged[:,i],.05)),p95_log2_tpm1=float(np.quantile(logged[:,i],.95)),p99_log2_tpm1=float(np.quantile(logged[:,i],.99)),
            correlation_gene_n=int(mask.sum()),median_pearson_other_samples=float(medians[i]) if len(ids)>1 else None,
            correlation_robust_z=robust,low_correlation_flag=bool(medians[i]<.8) if len(ids)>1 else False,
            robust_outlier_flag=bool(robust is not None and robust < -3),
            outlier_assessment='descriptive_only_n_less_than_5' if len(ids)<5 else 'review_flags_only_no_automatic_exclusion'))
    write_csv(PREFIX+'expression_sample_manifest.csv',manifest,list(manifest[0]))
    write_csv(PREFIX+'expression_gene_qc.csv',gene_qc,list(gene_qc[0]))
    write_csv(PREFIX+'expression_sample_qc.csv',sample_qc,list(sample_qc[0]))
    write_csv(PREFIX+'missing_patients.csv',missing,['patient_id','clinical_sample_id','reason'])
    write_csv(PREFIX+'sample_correlations.csv',[dict(expression_sample_id=s,**dict(zip(ids,correlations[i].tolist()))) for i,s in enumerate(ids)],['expression_sample_id']+ids)
    write_csv('outputs/v2/predictors/candidate_signature.csv',predictions,['patient_id','sample_id','expression_sample_id','predictor_role','n_genes_scored','predictor_status']+PREDICTOR_COLUMNS)
    write_json('outputs/v2/predictors/standardization_parameters.json',dict(status='fitted_provisional_partial_download' if missing else 'fitted',reference_samples=ids,reference_patient_ids=[r['patient_id'] for r in selected],ddof=0,expression_field='tpm_unstranded',transform='log2(TPM+1)',genes=parameters))
    summary_path=root/'outputs/v2/cohort/cohort_summary.json';summary=json.loads(summary_path.read_text())
    problems=[]
    if missing:problems.append(f'{len(missing)} frozen patients have no complete valid matching local RNA file')
    unresolved=[r['file_id'] for r in selected if not r['aliquot_id']]
    if unresolved:problems.append('Aliquot barcode unavailable for: '+','.join(unresolved))
    summary.update(status='predictors_computed_partial_RNA' if missing else 'complete',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        clinical_tnbc_n=len(cohort),expression_available_n=len(ids),predictor_complete_n=sum(p['predictor_status']=='complete' for p in predictions),
        available_candidate_genes=[g for g in config.CANDIDATE_GENES if g in expression],missing_from_verified_expression=[g for g in config.CANDIDATE_GENES if g not in expression],
        unassessable_candidate_genes=[],available_comparison_genes=[g for g in config.COMPARISON_GENES if g in expression],
        missing_comparison_genes=[g for g in config.COMPARISON_GENES if g not in expression],
        expression_file=str(matrix_path),expression_sha256=sha256(matrix_path),expression_field='tpm_unstranded',expression_transform='log2(TPM+1)',
        expression_provenance='GDC STAR - Counts; GRCh38; '+str(comments),blocking_problems=problems,
        frozen_cohort_sha256=frozen_hash,local_complete_STAR_n=len(local_paths),ignored_partial_n=len(list(config.STAR_DIR.rglob('*.part'))),
        duplicate_patient_groups={k:len(v) for k,v in groups.items() if len(v)>1},numpy_version=np.__version__,
        original_cohort_csv_unchanged=True,clinical_csv_expression_columns='historical Stage 1 availability; authoritative current match is expression_sample_manifest.csv',
        genomic_outcomes_defined=False,downloads_performed=False)
    write_json('outputs/v2/cohort/cohort_summary.json',summary)
    write_json(PREFIX+'expression_run_manifest.json',dict(cohort_sha256=frozen_hash,inventory_sha256=sha256(inventory_path),
        metadata_sha256=sha256(detail_path) if detail_path.exists() else None,complete_file_paths_at_start=[str(p) for p in local_paths],
        source_hashes={r['source_path']:r['sha256'] for r in selected},source_release=inventory['release'],
        metadata_API='https://api.gdc.cancer.gov/files/{file_id}?expand=cases.samples.portions.analytes.aliquots',
        metadata_only_lookup=True,biological_data_downloaded=False,annotation=comments))
    write_methods(summary,selected,missing,mask.sum(),gene_qc)
    if sha256(frozen_path)!=frozen_hash:raise RuntimeError('Frozen clinical cohort changed')
    print(json.dumps({k:summary[k] for k in ['clinical_tnbc_n','expression_available_n','predictor_complete_n','available_candidate_genes','missing_from_verified_expression','blocking_problems']},indent=2))
    return summary

def write_methods(s,selected,missing,corr_n,gene_qc):
    text=f'''# V2 Stage 1: clinical cohort and GDC STAR expression predictors

Run `python -m src.stage1` (or `python -m src.star`) from V2. This continuation reads the frozen clinical cohort directly; it never re-runs clinical selection or the local dataset search. The original source notebooks and data are unchanged.

## Frozen clinical definition and identifiers

`outputs/v2/cohort/tcga_tnbc_patients.csv` retains **{s['clinical_tnbc_n']} patients** and its pre-run SHA256 `{s['frozen_cohort_sha256']}`. It came from `{s['clinical_file']}` using Negative ER (`breast_carcinoma_estrogen_receptor_status`), Negative PR (`breast_carcinoma_progesterone_receptor_status`) and Negative HER2 IHC (`lab_proc_her2_neu_immunohistochemistry_receptor_status`). Unknown, equivocal and indeterminate values were not counted negative. Positive HER2 ISH excluded two discordant cases from 116 triple-negative primary records. Missing/equivocal ISH did not invalidate an otherwise negative IHC. Primary tumor status was resolved using code 01 and concordant clinical sample-type annotations. No receptor RNA or PAM50 determined membership. This definition has not changed.

The frozen CSV's earlier expression-availability columns remain historical to preserve the file byte-for-byte. Current matches are authoritative in `expression_sample_manifest.csv` and `cohort_summary.json`, not those old columns.

## Actual STAR input and field

Only completed files under `data/external/tcga_rna/<GDC UUID>/` are considered. At run start there were {s['local_complete_STAR_n']} complete TSVs and {s['ignored_partial_n']} ignored partials. Only frozen-cohort matches are parsed. The acquisition inventory `outputs/acquisition/tcga_rna_files.json` supplies case/sample IDs, workflow, filenames, sizes and MD5s. Each matched input must pass filename/byte-count/MD5 validation. SHA256s and exact paths are recorded in `expression_sample_manifest.csv` and `expression_run_manifest.json`. A bounded GDC metadata-only request supplied the full aliquot barcodes in `gdc_aliquot_metadata.json`; no biological data were downloaded in this continuation.

The [GDC mRNA pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/) quantifies gene-level RNA on GRCh38, treats reads as unstranded for harmonization and supplies normalized TPM alongside counts. The observed files declare GENCODE v36. We use **`tpm_unstranded`**, then **log2(TPM+1)** for expression predictors and correlations. `unstranded` is retained for library/count QC, not standardized as abundance. Neither stranded-count column nor FPKM is substituted. `N_*` rows are STAR read-accounting summaries, excluded from the gene matrix and retained for QC. Gene Ensembl IDs (including versions), symbols and biotypes are preserved exactly; no exon reconstruction, gene-coordinate inference or alias guessing occurs.

## Matching and replicate rule

Require a single GDC case matching a frozen patient, a unique primary-tumor sample and agreement of the patient-plus-sample-type key with the frozen clinical sample. If the clinical ID includes a vial, the vial must match. Full GDC aliquot IDs are preserved and validated against the sample; ambiguous identity or mismatched metadata excludes the file for review. A missing detailed aliquot barcode is explicitly flagged and uses sample plus file UUID rather than an invented barcode.

After identity/integrity checks, select at most one file per patient: highest total assigned `unstranded` gene reads, then lexical sample ID, aliquot ID and file UUID. This is a documented project technical tie-break, not an official TCGA biological-quality ranking. It does not inspect candidate activity or outcomes. No duplicate groups occurred in this run: {s['duplicate_patient_groups']}. The local download contains no external pathology exclusion resource or RIN metrics, so those checks are unavailable; an MD5 pass is not a claim of comprehensive biological QC. Correlation flags are review-only and do not silently change clinical membership.

Matched **{s['expression_available_n']} of {s['clinical_tnbc_n']}** patients. Missing **{len(missing)}**, listed individually in `missing_patients.csv`. A missing local file is not evidence that the patient lacks a GDC assay. Only these matched discovery samples contribute to standardization. No non-TNBC sample is used to fit parameters.

## Matrix and predictor computation

`tcga_tnbc_expression.tsv.gz` is genome-wide gene-by-sample **log2(TPM+1)**, with leading Ensembl ID, gene symbol and gene-type columns; expression columns use aliquot barcodes (or an explicitly flagged sample/file key). The sample manifest links each column to the original clinical sample and patient. All selected files must have identical annotation/order; the runner fails instead of silently intersecting incompatible annotations.

Candidate genes: {', '.join(config.CANDIDATE_GENES)}. Available: {', '.join(s['available_candidate_genes'])}. Missing or ambiguous: {s['missing_from_verified_expression']}.

Separate canonical comparisons: {', '.join(config.COMPARISON_GENES)}. Missing or ambiguous: {s['missing_comparison_genes']}. Duplicate symbols are not summed; any nonunique requested symbol is unscorable and reported in gene QC. Comparisons do not enter the nine-gene predictor.

The existing `fit_standardization` and `score_predictors` functions are reused unchanged. Fit each gene's mean and population SD (ddof=0) across the {s['expression_available_n']} matched TNBC discovery samples, then calculate z=(log2(TPM+1)-mean)/SD. Require at least two observed values and nonzero SD. Compute the equal-weight mean and median across all nine candidate z-scores, plus HORMAD1, STAG3 and REC8 single-gene z-scores. All nine must be usable for a composite. Missing values and constant genes are not set to zero. Frozen fitted parameters and exact reference sample/patient lists are recorded in `standardization_parameters.json`.

Complete nine-gene predictors: **{s['predictor_complete_n']}**. All 114 clinical patients remain in `candidate_signature.csv`; unmatched rows have blank values and a missing-file status. With only {s['expression_available_n']} RNA matches these fitted scores are **provisional and unsuitable for inference**. As additional discovery files become available, a new explicitly recorded fit will change the reference parameters; these partial-reference scores must not be mixed with a later complete-cohort fit.

## Expression QC

`expression_gene_qc.csv` reports every Ensembl row, symbol mapping count, candidate/comparison membership, zero counts, mean/median/max TPM, log-expression mean/SD and 5th/95th percentiles. Near-zero expression means TPM below 0.1 in at least 90% of matched samples. Near-zero log variation means SD at most 1e-8. These are descriptive flags; only exactly zero SD blocks the existing signature implementation. TPM zero is an observed zero, not a missing gene.

`expression_sample_qc.csv` reports detected genes, genes at TPM>=1, assigned reads, assigned fraction of the available STAR accounting categories, TPM sums, mean/median and 5th/95th/99th percentile log expression. The STAR accounting fraction is not a complete alignment-quality or purity estimate.

Pearson correlations use {int(corr_n)} protein-coding genes with TPM>=1 in at least one matched sample and log SD>1e-8. `sample_correlations.csv` contains the full symmetric matrix. Median correlation to other samples below 0.8 is flagged for review. A robust median/MAD z-score below -3 is additionally flagged only when N>=5 and MAD>0; it is unavailable at N=3. These are prespecified descriptive checks, not validated exclusion thresholds. No correlation outlier was automatically removed.

The package does not construct genomic outcomes from candidates. No association testing, pathway enrichment, survival analysis, immune analysis or retired composite calculation is performed. Stage 1 tests are run with `python -m unittest discover -s tests -p test_stage1*.py -v`; `outputs/v2/stage1_test_report.json` records results.
'''
    output_path('docs/V2_METHODS_STAGE1.md').write_text(text,encoding='utf-8')

if __name__=='__main__':run()
