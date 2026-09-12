import csv
import gzip
import hashlib
import json
import re
from pathlib import Path
from .config import PROJECT_ROOT

def open_text(path):
    # Compression magic detection
    with open(path, 'rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt', encoding='utf-8-sig', newline='') if compressed else open(path, encoding='utf-8-sig', newline='')

def tcga_ids(value):
    sample = value.strip().upper()
    match = re.fullmatch(r'(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})(?:-(\d{2})([A-Z])?(?:-[A-Z0-9-]+)?)?', sample)
    if not match:
        raise ValueError(f'Invalid TCGA identifier: {value!r}')
    patient, code, vial = match.groups()
    return dict(patient_id=patient, normalized_sample_id=sample, sample_type_code=code or '',
                sample_key=patient + '-' + code if code else '', vial=vial or '')

def sha256(path):
    result = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): result.update(block)
    return result.hexdigest()

def output_path(relative):
    root = PROJECT_ROOT.resolve()
    p = (root / relative).resolve()
    if not p.is_relative_to(root) or p == root:
        raise ValueError('Output must remain inside V2')
    p.parent.mkdir(parents=True, exist_ok=True)
    return p

def write_json(relative, data):
    output_path(relative).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')

def write_csv(relative, rows, columns):
    with output_path(relative).open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader(); writer.writerows(rows)
