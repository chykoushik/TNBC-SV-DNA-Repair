# Structural Variation and DNA Repair Gene Expression in Triple-Negative Breast Cancer

Code and reproducible analysis outputs for the computational study:

**A Computational Multi-omic Screen of Structural Variation and Gene Disruption at DNA Repair and Cohesin Loci in Triple-Negative Breast Cancer**

## Overview

This repository contains the analysis pipeline used to evaluate a fixed nine-gene meiosis-related expression score in triple-negative breast cancer (TNBC) and its association with genomic homologous recombination deficiency (HRD) and related genomic phenotypes.

The candidate gene set is:

- HORMAD1
- HORMAD2
- STAG3
- REC8
- SMC1B
- SYCP2
- SYCP3
- MSH4
- MSH5

The score is calculated by z-standardizing expression of each candidate gene within the specified cohort and taking the equal-weight mean of the nine standardized values.

Established homologous recombination and somatic cohesin genes are analysed separately as comparison genes and are not included in the candidate score.

The repository also contains the locus-level structural variation analysis across the broader DNA repair and cohesin gene panel.

## Analyses

The final analysis includes:

- Population structural variation at DNA repair and cohesin loci using gnomAD v4.1 and HGSVC2.
- TCGA-BRCA TNBC analysis using clinically defined ER-negative, PR-negative, and HER2-negative primary tumours.
- Association of the nine-gene candidate expression score with TCGA genomic HRD.
- Independent analysis in SCAN-B using HRDetect probability and scarHRD score.
- Proliferation-adjusted sensitivity analyses.
- Leave-one-gene-out and null-signature analyses.
- METABRIC analysis of copy-number alteration burden.
- GSE25066 secondary analyses of pathological complete response and distant relapse-free survival using the available four-gene subset.
- DepMap CRISPR dependency analyses.
- GSE176078 single-cell analysis.

## Repository Structure

```text
TNBC-SV-DNA-Repair/
├── config/
│   └── data_sources.yaml
├── src/
│   ├── cohort.py
│   ├── config.py
│   ├── covariates.py
│   ├── data_loading.py
│   ├── expression.py
│   ├── figures_stage2.py
│   ├── figures_stage3.py
│   ├── null_models.py
│   ├── outcomes.py
│   ├── signatures.py
│   ├── stage1.py
│   ├── stage2.py
│   ├── stage3.py
│   ├── stage3_functional.py
│   ├── stage3_patient.py
│   ├── stage3_single_cell.py
│   ├── stage4.py
│   ├── stage4_scanb.py
│   ├── star.py
│   ├── statistics.py
│   └── utils.py
├── scripts/
│   ├── build_final_package.py
│   ├── download_public_data.py
│   ├── validate_final_package.py
│   └── verify_downloads.py
├── tests/
├── outputs/
│   └── v2/
│       └── final/
│           ├── figures/
│           └── tables/
├── environment.yml
├── requirements.txt
└── README.md
```

## Candidate and Comparison Gene Sets

### Candidate expression score

The fixed candidate set consists of nine meiosis-related genes:

```text
HORMAD1
HORMAD2
STAG3
REC8
SMC1B
SYCP2
SYCP3
MSH4
MSH5
```

These genes were selected from the biological framework of meiotic chromosome cohesion, chromosome pairing, and recombination. The set was fixed before testing the genomic outcomes.

### Comparison genes

Established homologous recombination and somatic cohesin genes are analysed separately:

```text
BRCA1
BRCA2
PALB2
RAD51
RAD51B
RAD51C
RAD51D
BRIP1
ATM
CHEK2
STAG2
RAD21
SMC1A
```

MLH3 is retained only in the locus-level structural variation analysis.

## Cohorts and Data Sources

### TCGA-BRCA

TCGA-BRCA is used for the primary genomic HRD analysis.

TNBC is defined using clinical ER, PR, and HER2 receptor status. ESR1, PGR, and ERBB2 expression thresholds are not used to define the final TNBC cohort.

Expression data are based on GDC augmented STAR counts. Candidate-gene expression is transformed and standardized before calculation of the equal-weight candidate score.

### SCAN-B

SCAN-B provides an independent patient cohort with genomic HRD measurements derived from whole-genome sequencing.

The principal endpoints are:

- HRDetect probability
- scarHRD score

The same nine-gene equal-weight expression score is used.

### METABRIC

METABRIC is used to examine the relationship between the candidate expression score and copy-number alteration burden and for additional patient-level analyses.

### GSE25066

GSE25066 is used for secondary clinical analyses of:

- pathological complete response
- distant relapse-free survival

Only four of the nine candidate genes are available on the required expression platform. The GSE25066 analysis therefore uses a reduced four-gene score and is not treated as validation of the complete nine-gene score.

### DepMap

DepMap CRISPR dependency data are used to test whether the candidate expression program is associated with selective genetic dependencies.

### GSE176078

GSE176078 single-cell RNA-seq data are used for secondary analysis of candidate-program expression in malignant and nonmalignant cellular compartments.

### Population structural variation

Population structural variation at the analysed loci is evaluated using:

- gnomAD v4.1 structural variants
- HGSVC2 structural variants

## Analysis Pipeline

The current study is implemented as a Python pipeline rather than the notebook workflow used in earlier development versions of the project.

The main analysis modules are contained in:

```text
src/
```

Supporting execution, data acquisition, validation, and final-output scripts are contained in:

```text
scripts/
```

Automated checks are contained in:

```text
tests/
```

## Reproducibility

Install the required Python dependencies using:

```bash
pip install -r requirements.txt
```

or create the specified environment using:

```bash
conda env create -f environment.yml
```

Public data-source definitions are provided in:

```text
config/data_sources.yaml
```

The repository does not redistribute controlled or externally hosted source datasets. Source datasets should be obtained from their respective repositories according to their access and licensing requirements.

## Main Outputs

The final analysis package contains manuscript-facing figures and tables under:

```text
outputs/v2/final/
```

The analysis pipeline also retains intermediate statistical outputs required to trace reported results to their source analyses.

The principal analyses include the TCGA primary HRD model, SCAN-B replication analyses, proliferation-adjusted models, leave-one-gene-out analyses, null-signature controls, METABRIC analyses, GSE25066 clinical analyses, DepMap dependency analyses, and GSE176078 single-cell analyses.

## Data Availability

The analyses use publicly available data from TCGA, METABRIC, SCAN-B, gnomAD, HGSVC2, GSE25066, GSE176078, DepMap, GTEx, and COSMIC where applicable.

Data remain subject to the terms and access requirements of their original providers.

## License

See the repository license for terms governing reuse of the code.
