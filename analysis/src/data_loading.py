import csv
from .config import CLINICAL_FILE, ER_COLUMN, PR_COLUMN, HER2_COLUMN
from .utils import open_text

def load_clinical(path=CLINICAL_FILE):
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'sampleID', ER_COLUMN, PR_COLUMN, HER2_COLUMN}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Missing required clinical columns: {required-set(reader.fieldnames or [])}')
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'Malformed clinical row {line}')
            rows.append(row)
    return rows
