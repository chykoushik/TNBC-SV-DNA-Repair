from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATASET_DIR = PROJECT_ROOT.parent.parent / 'dataset'
CLINICAL_FILE = DATASET_DIR / 'BRCA_clinicalMatrix.tsv.gz'
OUTPUT_DIR = PROJECT_ROOT / 'outputs' / 'v2'
STAR_DIR = PROJECT_ROOT / 'data' / 'external' / 'tcga_rna'
STAR_EXPRESSION_FIELD = 'tpm_unstranded'
# Verified gene matrix
EXPRESSION_FILE = None
EXPRESSION_PROVENANCE = None
EXPRESSION_TRANSFORM = None  # Expression transforms
ER_COLUMN = 'breast_carcinoma_estrogen_receptor_status'
PR_COLUMN = 'breast_carcinoma_progesterone_receptor_status'
HER2_COLUMN = 'lab_proc_her2_neu_immunohistochemistry_receptor_status'
HER2_ISH_COLUMN = 'lab_procedure_her2_neu_in_situ_hybrid_outcome_type'
CANDIDATE_GENES = ('HORMAD1', 'HORMAD2', 'STAG3', 'REC8', 'SMC1B', 'SYCP2', 'SYCP3', 'MSH4', 'MSH5')
COMPARISON_GENES = ('BRCA1', 'BRCA2', 'PALB2', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'BRIP1', 'ATM', 'CHEK2', 'STAG2', 'RAD21', 'SMC1A')
