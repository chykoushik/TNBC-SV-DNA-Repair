import datetime
import platform
from collections import Counter
from . import config
from .cohort import build_cohort
from .data_loading import load_clinical
from .expression import discover_expression, load_gene_expression, match_samples
from .signatures import fit_standardization, score_predictors, PREDICTOR_COLUMNS
from .utils import write_csv, write_json, sha256, output_path

def build_initial_cohort():
    clinical_hash = sha256(config.CLINICAL_FILE)
    discovery = discover_expression(config.DATASET_DIR)
    write_csv('outputs/v2/cohort/local_expression_search.csv', discovery, list(discovery[0]))
    clinical = load_clinical(config.CLINICAL_FILE)
    cohort, ledger = build_cohort(clinical)
    expression = {}; matches = {}; expression_hash = None
    if config.EXPRESSION_FILE is not None:
        resolved = config.EXPRESSION_FILE.resolve()
        if not resolved.is_relative_to(config.DATASET_DIR.resolve()):
            raise ValueError('Stage 1 requires an existing expression file in the original dataset directory')
        expression_hash = sha256(resolved)
        samples, expression = load_gene_expression(resolved, config.CANDIDATE_GENES + config.COMPARISON_GENES,
                                                  provenance=config.EXPRESSION_PROVENANCE, transform=config.EXPRESSION_TRANSFORM)
        matches = match_samples(cohort, samples)
    parameters = fit_standardization(expression, matches.values())
    scores = {r['expression_sample_id']:r for r in score_predictors(expression, matches.values(), parameters)}
    prediction_rows = []
    for row in cohort:
        expression_id = matches.get(row['sample_id'], '')
        row.update(expression_available=bool(expression_id), expression_sample_id=expression_id)
        predicted = dict(patient_id=row['patient_id'], sample_id=row['sample_id'], expression_sample_id=expression_id,
                         predictor_role='expression_predictor', n_genes_scored=0,
                         predictor_status='blocked_no_trustworthy_local_gene_expression' if config.EXPRESSION_FILE is None else 'no_unambiguous_expression_match')
        predicted.update({column:None for column in PREDICTOR_COLUMNS})
        if expression_id: predicted.update(scores[expression_id])
        prediction_rows.append(predicted)
    available = [g for g in config.CANDIDATE_GENES if g in expression]
    missing = [g for g in config.CANDIDATE_GENES if g not in expression] if config.EXPRESSION_FILE else []
    blocked = []
    if config.EXPRESSION_FILE is None:
        blocked.append('No provenance-verified local TCGA-BRCA gene-level RNA matrix. The local BRCA matrix is exon-level. Predictor values cannot be calculated without prohibited reconstruction or new data.')
    if config.EXPRESSION_FILE and not matches: blocked.append('No unambiguous primary-tumor expression matches')
    summary = dict(stage='1', status='cohort_complete_predictors_blocked' if blocked else 'complete',
                   created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), python_version=platform.python_version(),
                   clinical_file=str(config.CLINICAL_FILE), clinical_sha256=clinical_hash,
                   clinical_records=len(clinical), primary_records=sum(x['primary_status']=='primary' for x in ledger),
                   clinical_tnbc_n=len(cohort), expression_available_n=len(matches),
                   strict_three_negative_primary_before_ish_conflict=sum(x['primary_status']=='primary' and all(x[k]=='negative' for k in ['er_status','pr_status','her2_ihc_status']) for x in ledger),
                   exclusion_counts=dict(sorted(Counter(x['exclusion_reason'] for x in ledger if not x['included']).items())),
                   candidate_genes=list(config.CANDIDATE_GENES), comparison_genes=list(config.COMPARISON_GENES),
                   available_candidate_genes=available, missing_from_verified_expression=missing,
                   unassessable_candidate_genes=list(config.CANDIDATE_GENES) if config.EXPRESSION_FILE is None else [],
                   expression_file=str(config.EXPRESSION_FILE) if config.EXPRESSION_FILE else None,
                   expression_sha256=expression_hash, expression_transform=config.EXPRESSION_TRANSFORM,
                   expression_provenance=config.EXPRESSION_PROVENANCE,
                   predictor_complete_n=sum(x['predictor_status']=='complete' for x in prediction_rows),
                   blocking_problems=blocked, genomic_outcomes_defined=False,
                   local_files_searched=len(discovery), downloads_performed=False)
    write_csv('outputs/v2/cohort/tcga_tnbc_patients.csv', cohort, list(cohort[0]) if cohort else list(ledger[0])+['expression_available','expression_sample_id'])
    write_csv('outputs/v2/cohort/cohort_decisions.csv', ledger, list(ledger[0]))
    write_json('outputs/v2/cohort/cohort_summary.json', summary)
    write_csv('outputs/v2/predictors/candidate_signature.csv', prediction_rows,
              ['patient_id','sample_id','expression_sample_id','predictor_role','n_genes_scored','predictor_status']+PREDICTOR_COLUMNS)
    write_json('outputs/v2/predictors/standardization_parameters.json',
               dict(status='not_fitted_no_expression' if not matches else 'fitted', reference_samples=sorted(matches.values()), ddof=0, genes=parameters))
    write_json('outputs/v2/predictors/gene_sets.json',dict(candidate=list(config.CANDIDATE_GENES),canonical_comparison=list(config.COMPARISON_GENES)))
    assert clinical_hash == sha256(config.CLINICAL_FILE), 'Clinical source changed during run'
    if config.EXPRESSION_FILE: assert expression_hash == sha256(config.EXPRESSION_FILE), 'Expression source changed during run'
    methods = f'''# V2 methods: Stage 1

This is an implementation run against existing local data, not a new audit or a biological association analysis. Run from the V2 root with `python -m src.stage1`. The six requested modules, `src/__init__.py` and the Stage 1 entry point form an importable Python package. Only the Python standard library is needed for this stage; the recorded runtime was Python {platform.python_version()}. Existing environment pins were not installed or modified.

## Clinical TNBC definition

Input: `{config.CLINICAL_FILE}`. SHA256: `{clinical_hash}`. Despite its `.gz` suffix this file is plain text; loading selects compression by file magic. The input contains {len(clinical)} rows; {summary['primary_records']} are resolved primary tumor records.

Inclusion requires source-reported **Negative** in each of:

- ER: `{config.ER_COLUMN}`.
- PR: `{config.PR_COLUMN}`.
- HER2 IHC: `{config.HER2_COLUMN}`.

Values are stripped and compared case-insensitively. Blank, unknown, indeterminate and equivocal receptor values are not negative. This strict stage does not rescue equivocal or missing IHC with negative ISH. A positive `{config.HER2_ISH_COLUMN}` contradicts negative IHC and excludes the record for adjudication. A negative IHC with missing/equivocal ISH remains eligible because its qualifying negative evidence is IHC, not an assumed negative ISH. This conservative rule may exclude clinically adjudicable cases; future adjudication must be documented separately.

Primary tumors are resolved from TCGA sample-type code 01, `sample_type_id` and/or the primary-tumor label. At least one source must resolve primary status, and all available resolved evidence must agree. Non-primary, unresolved and conflicting records are excluded. Normal code 11 and metastatic code 06 are not included. The code preserves original sample IDs and normalized IDs separately and extracts the 12-character TCGA patient identifier.

There were **{summary['strict_three_negative_primary_before_ish_conflict']}** primary records negative for all three receptors before ISH conflict screening. **{summary['exclusion_counts'].get('her2_ihc_ish_conflict',0)}** were excluded for positive ISH despite negative IHC. Final clinical TNBC: **{len(cohort)} unique patients**. The full ledger `outputs/v2/cohort/cohort_decisions.csv` records every inclusion/exclusion. For multiple primary records per patient, discordant receptor records are excluded; otherwise the lexicographically first normalized eligible sample is selected. The local result has one included sample per patient.

ESR1, PGR and ERBB2 expression, receptor-expression medians and PAM50 are never used for selection. No GIPS calculation or outcome construction is present in this package.

## Local expression search and actual expression used

`config.py` references `{config.DATASET_DIR}` directly. No source data were copied or edited. The recursive search inspected {len(discovery)} files using filenames and bounded text prefixes; archives were not extracted. Its evidence is in `outputs/v2/cohort/local_expression_search.csv`.

**Actual expression input: {str(config.EXPRESSION_FILE) if config.EXPRESSION_FILE else 'none; no trustworthy TCGA-BRCA gene-level expression matrix found locally'}.**

`TCGA.BRCA.sampleMap_HiSeqV2_exon.gz` is exon-level and was rejected. `TCGA.OV.sampleMap_HiSeqV2.gz` belongs to ovarian cancer. The BRCA GISTIC matrix measures copy number, not RNA. METABRIC and GEO files represent other cohorts and cannot provide expression for TCGA-BRCA patients. No exon-coordinate aggregation, coordinate hardcoding, surrogate RNA reconstruction or download was performed. **Expression-available clinical TNBC N: {len(matches)}.** This is unavailable suitable local RNA, not proof that these patients lack RNA assays elsewhere.

Candidate set: {', '.join(config.CANDIDATE_GENES)}.

Canonical comparison genes are kept in a distinct configuration tuple and do not enter the candidate signature: {', '.join(config.COMPARISON_GENES)}.

Available candidate genes in verified expression: {', '.join(available) or 'none assessable'}. Missing from a verified matrix: {', '.join(missing) or 'not assessable without an expression matrix'}. All nine genes are currently **unassessable**, not demonstrated biologically absent. The output CSV includes each clinical TNBC patient, blank predictor values and an explicit blocking status; blanks are not zeros or calculated scores.

## Implemented expression matching and predictor calculation

The loader requires an explicitly configured existing source, documented provenance and transform. It supports a symbol-indexed gene-by-sample TSV, streams rows and retains only candidate/comparison genes. It rejects exon files, duplicate genes/samples, malformed rows, and coordinate/Ensembl identifiers requiring an unimplemented verified mapping. No automatic gene alias or coordinate inference occurs. A documented already-logged matrix uses identity transformation; nonnegative abundance can use log2(value+1). The actual run selected neither because no valid matrix exists.

Matching uses an exact normalized sample ID first, then a unique patient-plus-two-digit-sample-type key. Only primary RNA samples are considered; a patient-only join is not allowed. Multiple aliquots at a fallback key are ambiguous and remain unmatched. This prevents row-order-dependent first-match selection. Both clinical and expression IDs are retained. A patient-only clinical record can enter a resolved clinical cohort but cannot be RNA-matched without a sample key.

The implemented reference set is the expression-matched clinical TNBC discovery cohort, one tumor per patient, with IDs sorted and stored in `standardization_parameters.json`. For each candidate gene, fit the reference mean and population SD (ddof=0) using observed values. At least two observations and positive SD are required. Then z=(transformed expression - reference mean)/reference SD. Fit and scoring are separate functions so held-out cohorts can reuse frozen parameters without refitting; cross-platform validation will need a separately prespecified calibration strategy.

The five implemented predictor columns are the equal-weight **mean of all nine gene z-scores**, **median of all nine gene z-scores**, and single-gene **HORMAD1_z**, **STAG3_z**, **REC8_z**. All nine usable values are required for either composite; missing/constant genes are not replaced with zeros or silently dropped. Single-gene predictors can remain available independently. Coverage and calculation status are recorded. These are expression predictors only. Candidate genes do not define HRD, CIN, SV or other outcomes. No standardization parameters were actually fitted and no predictor values were computed in this blocked local run.

## Results, checks and scope

Clinical TNBC N: **{len(cohort)}**. Expression-available N: **{len(matches)}**. Complete candidate predictor N: **{summary['predictor_complete_n']}**. Blocking problem: {'; '.join(blocked) or 'none'}.

Required artifacts are in `outputs/v2/cohort/` and `outputs/v2/predictors/`. Additional decision, search, gene-set and standardization files make all exclusions and missing scores explicit. Tests in `tests/test_stage1.py` verify clinical-only selection under changed receptor-expression values, unknown/equivocal handling, primary-only selection, positive-ISH conflict exclusion, deterministic IDs/duplicates, predictor arithmetic, strict missing coverage, separation from genomic outcomes and absence of the retired composite framework from executable source. Run with `python -m unittest discover -s tests -p test_stage1.py -v`; the test report is saved separately after execution.

No original notebook or dataset was modified. Clinical SHA256 was checked before and after execution. No HRD association, enrichment, immune deconvolution, survival modeling, downloading or GIPS was run. Stage 1 cohort construction is complete; numerical candidate predictors are blocked on suitable gene-level RNA, which cannot be supplied under the current no-download constraint.
'''
    output_path('docs/V2_METHODS_STAGE1.md').write_text(methods, encoding='utf-8')
    print(__import__('json').dumps(summary, indent=2))
    return summary

def run():
    from .star import run as continue_star
    return continue_star()

if __name__ == '__main__': run()
