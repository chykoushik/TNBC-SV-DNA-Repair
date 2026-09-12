import csv
import math
from pathlib import Path
from .utils import open_text, tcga_ids

def discover_expression(dataset_dir):
    report = []
    for path in sorted(Path(dataset_dir).rglob('*')):
        if not path.is_file(): continue
        name = path.name.lower()
        record = dict(path=str(path), size_bytes=path.stat().st_size, classification='', header_preview='', gene_preview='')
        if any(name.endswith(s) for s in ('.tar', '.tar.gz', '.zip')):
            record['classification'] = 'archive; cohort named in filename; not verified TCGA-BRCA gene-level RNA'
        elif any(t in name for t in ['vcf', '.bed', 'methyl', 'hm450']):
            record['classification'] = 'non-expression genomic data'
        else:
            try:
                with open_text(path) as f:
                    header = f.readline(65536)
                    second = f.readline(65536)
                record['header_preview'] = header[:300].strip()
                record['gene_preview'] = second[:150].strip()
                if 'exon' in name:
                    record['classification'] = 'rejected: exon-level expression; no gene reconstruction permitted'
                elif 'gistic' in name or 'copynumber' in name:
                    record['classification'] = 'rejected: copy number is not RNA expression'
                elif 'tcga.ov' in name:
                    record['classification'] = 'rejected: ovarian cohort, not TCGA-BRCA'
                elif 'gse' in name or 'metabric' in name:
                    record['classification'] = 'external cohort; not TCGA-BRCA RNA'
                elif ('brca' in name or 'pancan' in name) and any(s in name for s in ['rnaseq', 'hiseq', 'geneexp', 'expression', 'tpm']):
                    record['classification'] = 'possible gene-level RNA: requires explicit provenance, identifier and transform verification'
                else:
                    record['classification'] = 'not identified as TCGA-BRCA gene-level RNA'
            except (UnicodeError, OSError) as exc:
                record['classification'] = 'unreadable text prefix: ' + str(exc)
        report.append(record)
    return report

def match_samples(cohort, expression_ids):
    exact, key = {}, {}
    for raw in sorted(expression_ids):
        ids = tcga_ids(raw)
        if ids['sample_type_code'] != '01': continue
        exact.setdefault(ids['normalized_sample_id'], []).append(raw)
        key.setdefault(ids['sample_key'], []).append(raw)
    result = {}
    for row in cohort:
        ids = tcga_ids(row['sample_id'])
        choices = exact.get(ids['normalized_sample_id'], []) or key.get(ids['sample_key'], [])
        if len(choices) == 1: result[row['sample_id']] = choices[0]
    return result

def load_gene_expression(path, genes, *, provenance, transform):
    if not provenance or transform not in {'identity', 'log2p1'}:
        raise ValueError('Explicit expression provenance and transform required')
    if 'exon' in Path(path).name.lower(): raise ValueError('Exon matrices are unsupported')
    wanted = set(genes); values = {}; seen = set()
    with open_text(path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader); samples = header[1:]
        if not samples or len(samples) != len(set(samples)): raise ValueError('Missing/duplicate expression sample IDs')
        for s in samples: tcga_ids(s)
        for row in reader:
            if len(row) != len(header): raise ValueError('Expression row width mismatch')
            gene = row[0].strip()
            if ':' in gene or gene.startswith('ENSG'):
                raise ValueError('Coordinate/Ensembl matrix requires a verified gene annotation adapter')
            if gene in seen: raise ValueError('Duplicate gene identifier: ' + gene)
            seen.add(gene)
            if gene not in wanted: continue
            vector = []
            for cell in row[1:]:
                if cell.strip().lower() in {'', 'na', 'nan', 'null'}: vector.append(None); continue
                v = float(cell)
                if not math.isfinite(v): raise ValueError('Nonfinite expression')
                if transform == 'log2p1':
                    if v < 0: raise ValueError('Negative abundance cannot use log2p1')
                    v = math.log2(v + 1)
                vector.append(v)
            values[gene] = dict(zip(samples, vector))
    return samples, values
