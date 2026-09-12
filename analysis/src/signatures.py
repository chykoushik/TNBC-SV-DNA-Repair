import statistics
from .config import CANDIDATE_GENES

PREDICTOR_COLUMNS = ['candidate_mean_z', 'candidate_median_z', 'HORMAD1_z', 'STAG3_z', 'REC8_z']

def fit_standardization(expression, reference_samples):
    samples = sorted(set(reference_samples))
    parameters = {}
    for gene in CANDIDATE_GENES:
        observed = [expression.get(gene, {}).get(s) for s in samples]
        observed = [x for x in observed if x is not None]
        sd = statistics.pstdev(observed) if len(observed) >= 2 else 0
        parameters[gene] = dict(n=len(observed), mean=statistics.mean(observed) if observed else None,
                                sd=sd, usable=len(observed) >= 2 and sd > 0)
    return parameters

def score_predictors(expression, samples, parameters):
    rows = []
    for sample in sorted(set(samples)):
        zs = {}
        for gene in CANDIDATE_GENES:
            p = parameters[gene]; value = expression.get(gene, {}).get(sample)
            zs[gene] = (value-p['mean'])/p['sd'] if value is not None and p['usable'] else None
        complete = all(x is not None for x in zs.values())
        row = dict(expression_sample_id=sample, predictor_role='expression_predictor', n_genes_scored=sum(x is not None for x in zs.values()),
                   candidate_mean_z=statistics.mean(zs.values()) if complete else None,
                   candidate_median_z=statistics.median(zs.values()) if complete else None,
                   predictor_status='complete' if complete else 'incomplete_candidate_coverage')
        row.update({g+'_z':zs[g] for g in ['HORMAD1','STAG3','REC8']})
        rows.append(row)
    return rows
