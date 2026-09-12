from collections import defaultdict
from .config import ER_COLUMN, PR_COLUMN, HER2_COLUMN, HER2_ISH_COLUMN
from .utils import tcga_ids

def receptor_status(value):
    value = str(value or '').strip().lower()
    return value if value in {'negative', 'positive', 'equivocal', 'indeterminate'} else 'unknown'

def primary_status(row, ids):
    evidence = []
    if ids['sample_type_code']: evidence.append(ids['sample_type_code'] == '01')
    explicit = str(row.get('sample_type_id', '')).strip()
    if explicit.isdigit(): evidence.append(explicit.zfill(2) == '01')
    label = row.get('sample_type', '').strip().lower()
    if label in {'primary tumor', 'primary solid tumor'}: evidence.append(True)
    elif label and label not in {'na', 'nan', 'unknown', 'not reported'}: evidence.append(False)
    if not evidence: return 'unresolved'
    if any(evidence) and not all(evidence): return 'conflicting'
    return 'primary' if all(evidence) else 'non_primary'

def build_cohort(clinical_rows):
    ledger = []
    for source in clinical_rows:
        raw_id = source['sampleID']
        record = dict(sample_id=raw_id, patient_id='', normalized_sample_id='', sample_key='',
                      sample_type_code='', primary_status='', er_raw=source[ER_COLUMN], pr_raw=source[PR_COLUMN],
                      her2_ihc_raw=source[HER2_COLUMN], her2_ish_raw=source.get(HER2_ISH_COLUMN, ''),
                      er_status=receptor_status(source[ER_COLUMN]), pr_status=receptor_status(source[PR_COLUMN]),
                      her2_ihc_status=receptor_status(source[HER2_COLUMN]),
                      her2_ish_status=receptor_status(source.get(HER2_ISH_COLUMN)), included=False, exclusion_reason='')
        try:
            ids = tcga_ids(raw_id)
            record.update({k:v for k,v in ids.items() if k != 'vial'})
            record['primary_status'] = primary_status(source, ids)
            if record['primary_status'] != 'primary': record['exclusion_reason'] = record['primary_status']
            elif not all(record[c] == 'negative' for c in ['er_status', 'pr_status', 'her2_ihc_status']):
                record['exclusion_reason'] = 'not_all_three_receptors_negative'
            elif record['her2_ish_status'] == 'positive': record['exclusion_reason'] = 'her2_ihc_ish_conflict'
            else: record['included'] = True
        except ValueError:
            record['exclusion_reason'] = 'invalid_sample_id'
        ledger.append(record)
    groups = defaultdict(list)
    for r in ledger:
        if r['primary_status'] == 'primary': groups[r['patient_id']].append(r)
    for group in groups.values():
        statuses = {tuple(r[k] for k in ['er_status', 'pr_status', 'her2_ihc_status', 'her2_ish_status']) for r in group}
        if len(statuses) > 1:
            for r in group: r.update(included=False, exclusion_reason='discordant_primary_receptor_records')
        else:
            eligible = sorted((r for r in group if r['included']), key=lambda r:(r['normalized_sample_id'], r['sample_id']))
            for r in eligible[1:]: r.update(included=False, exclusion_reason='duplicate_patient_primary_sample')
    ledger.sort(key=lambda r:(r['patient_id'],r['normalized_sample_id'],r['sample_id'],r['exclusion_reason']))
    return [dict(r) for r in ledger if r['included']], ledger
