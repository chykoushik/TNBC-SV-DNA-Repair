# Computational Analysis of Meiotic and Cohesin Gene Expression and Homologous Recombination Deficiency in Triple-Negative Breast Cancer

This repository contains the code and reproducibility materials for the study:

**Computational Analysis of Meiotic and Cohesin Gene Expression and Homologous Recombination Deficiency in Triple-Negative Breast Cancer**

The study evaluates whether expression of a predefined nine-gene program composed of meiosis-related and meiotic-cohesin genes is associated with genomic homologous recombination deficiency (HRD) in triple-negative breast cancer (TNBC).

The primary analysis uses TCGA-BRCA, with independent genomic evaluation in SCAN-B. Additional analyses use METABRIC, GSE25066, DepMap, and GSE176078. Population structural-variation evidence at DNA repair, cohesin, and meiosis-related loci is evaluated using gnomAD v4.1 and HGSVC2.

## Repository organization

The repository contains two analysis generations:

```text
TNBC-SV-DNA-Repair/
│
├── README.md
│
├── analysis/
│   ├── config/
│   ├── src/
│   ├── scripts/
│   ├── tests/
│   ├── outputs/
│   ├── results/
│   ├── environment.yml
│   └── requirements.txt
│
└── legacy/
    ├── notebook1_setup_data.ipynb
    ├── notebook2_sv_gips.ipynb
    ├── notebook3_immune_eqtl.ipynb
    ├── notebook4_survival_validation.ipynb
    ├── notebook5_updated_brca.ipynb
    ├── notebook6_new_cohorts.ipynb
    ├── notebook7_clinical_cohort_validation.ipynb
    └── results/
```

The `analysis/` directory contains the reproducible analysis pipeline corresponding to the study described above.

The `legacy/` directory preserves the earlier notebook-based exploratory workflow for provenance. These notebooks include analyses developed before the final study design and are not the primary workflow for reproducing the reported results.

## Study design

The analysis centers on a fixed nine-gene expression program tested against genomic measurements of HRD.

The principal analysis is performed in clinically defined TCGA-BRCA TNBC. Independent genomic evaluation is performed in SCAN-B using HRDetect and scarHRD.

Additional datasets address complementary questions rather than serving as interchangeable replications of the primary endpoint:

- METABRIC is used for copy-number burden and secondary patient-level analyses.
- GSE25066 is used for pathological complete response and distant relapse-free survival.
- DepMap is used for CRISPR dependency analysis.
- GSE176078 is used for single-cell analysis.
- gnomAD v4.1 and HGSVC2 are used for population structural-variation analysis.

## Candidate gene set

The predefined candidate expression set contains nine genes:

| Gene | Biological context |
|---|---|
| `HORMAD1` | Meiosis-related chromosome surveillance |
| `HORMAD2` | Meiosis-related chromosome surveillance |
| `STAG3` | Meiotic cohesin |
| `REC8` | Meiotic cohesin |
| `SMC1B` | Meiotic cohesin |
| `SYCP2` | Synaptonemal complex |
| `SYCP3` | Synaptonemal complex |
| `MSH4` | Meiotic recombination |
| `MSH5` | Meiotic recombination |

The candidate set was fixed before testing the genomic outcomes and was not selected according to associations observed in the analysed cohorts.

For the principal patient-level analyses, expression of each candidate gene is standardized within the specified cohort reference and the candidate score is calculated as the equal-weight mean of the standardized expression values.

The complete nine-gene score is used for the principal TCGA and SCAN-B analyses and for the relevant METABRIC and DepMap analyses.

GSE25066 is an exception because only four candidate genes are available under the required platform mapping. Its results therefore represent a reduced four-gene analysis rather than validation of the complete nine-gene score.

## Comparison gene set

Established homologous-recombination and somatic-cohesin genes are maintained as a separate comparison set:

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

These genes are not components of the nine-gene candidate score.

`MLH3` is retained only in the locus-level structural-variation analysis.

## Analyses

### TCGA-BRCA

TCGA-BRCA provides the primary analysis.

TNBC is defined using clinical ER, PR, and HER2 receptor information. Expression thresholds for `ESR1`, `PGR`, and `ERBB2` are not used to define TNBC.

The primary outcome is continuous genomic HRD score.

The principal model tests the association between the fixed nine-gene candidate expression score and genomic HRD.

Robustness and sensitivity analyses include:

- HC3 heteroscedasticity-consistent inference
- Spearman rank correlation
- proliferation-adjusted regression
- individual-gene analyses
- leave-one-gene-out analyses
- alternative candidate-score summaries
- permutation testing
- matched random-gene null analyses
- comparison with established HR and cohesin genes

### SCAN-B

SCAN-B provides independent genomic evaluation of the candidate expression program.

The principal outcomes are:

- HRDetect probability
- scarHRD score

SCAN-B gene-expression measurements are obtained from GEO accession `GSE96058` and matched to published whole-genome genomic measurements using SCAN-B sample identifiers.

The same fixed nine-gene equal-weight expression score is evaluated against both genomic HRD measurements.

Additional analyses include proliferation adjustment, secondary genomic outcomes, leave-one-gene-out analysis, permutation testing, matched random-gene signatures, and bootstrap analyses.

### METABRIC

METABRIC provides an additional patient-level analysis using expression, copy-number, clinical, and survival information.

The principal genomic analysis evaluates the relationship between the candidate expression score and copy-number alteration burden.

Proliferation-adjusted models are used to assess whether this association is independent of proliferative state.

Survival analyses are secondary.

### GSE25066

GSE25066 provides expression and clinical outcome data from patients receiving neoadjuvant chemotherapy.

The analysed outcomes are:

- pathological complete response
- distant relapse-free survival

Only four members of the nine-gene candidate set are available under the required platform mapping. No missing candidate genes are imputed.

The resulting four-gene score is therefore treated as a reduced secondary analysis.

### DepMap

DepMap model metadata, gene-expression data, and CRISPR gene-effect measurements are used to evaluate whether the candidate expression program is associated with selective genetic dependencies in TNBC models.

Genome-wide dependency analyses use multiple-testing correction.

### GSE176078

GSE176078 single-cell RNA-seq data are used to evaluate expression of the candidate program across malignant epithelial and nonmalignant cellular compartments.

Donor-level pseudobulk measurements are used for inference so that individual cells are not treated as independent biological replicates.

### Population structural variation

Population structural-variation evidence at DNA repair, cohesin, and meiosis-related loci is evaluated using:

- gnomAD v4.1 structural variants
- HGSVC2 structural variants

This locus-level analysis is separate from the patient-level expression-score analyses.

## Analysis directory

The reproducible Python pipeline is located in:

```text
analysis/
```

Its main structure is:

```text
analysis/
├── config/
│   └── data_sources.yaml
│
├── src/
│   ├── __init__.py
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
│
├── scripts/
│   ├── build_final_package.py
│   ├── download_public_data.py
│   ├── validate_final_package.py
│   └── verify_downloads.py
│
├── tests/
├── outputs/
├── results/
├── environment.yml
└── requirements.txt
```

## Data sources

Raw source datasets are not redistributed in this repository. Data should be obtained directly from the original providers and remain subject to their respective access requirements, licenses, and terms of use.

| Resource | Use in analysis | Access |
|---|---|---|
| TCGA-BRCA | Primary TNBC cohort, clinical information, expression and genomic HRD | [UCSC Xena](https://xenabrowser.net/) |
| SCAN-B genomic data | HRDetect, scarHRD and related whole-genome measurements | Published SCAN-B genomic study and associated data |
| GSE96058 | SCAN-B expression | [NCBI GEO GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) |
| METABRIC | Expression, copy-number, clinical and survival analyses | [cBioPortal METABRIC](https://www.cbioportal.org/study/summary?id=brca_metabric) |
| GSE25066 | Neoadjuvant response and DRFS analyses | [NCBI GEO GSE25066](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066) |
| DepMap | Expression, model metadata and CRISPR dependency | [DepMap Portal](https://depmap.org/portal/) |
| GSE176078 | Single-cell breast cancer transcriptomics | [NCBI GEO GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) |
| gnomAD v4.1 SV | Population structural variation | [gnomAD Downloads](https://gnomad.broadinstitute.org/downloads) |
| HGSVC2 | Population structural variation | [HGSVC2 Data](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/) |
| COSMIC Cancer Gene Census | Cancer-gene annotation | [COSMIC Cancer Gene Census](https://cancer.sanger.ac.uk/census) |

### TCGA-BRCA

TCGA-BRCA data are obtained through UCSC Xena:

```text
https://xenabrowser.net/
```

The analysis uses clinical receptor annotations for TNBC definition and GDC augmented STAR expression data for the expression analyses.

### SCAN-B

SCAN-B expression data are available through GEO accession:

```text
GSE96058
```

Source:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
```

The corresponding HRDetect, scarHRD, and related genomic measurements are obtained from the published SCAN-B whole-genome analysis and its associated supplementary data.

### METABRIC

METABRIC data are obtained through cBioPortal:

```text
https://www.cbioportal.org/study/summary?id=brca_metabric
```

Study identifier:

```text
brca_metabric
```

### GSE25066

GSE25066 is available through NCBI GEO:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066
```

Accession:

```text
GSE25066
```

### DepMap

DepMap data are available from:

```text
https://depmap.org/portal/
```

The analysis uses model metadata, gene expression, and CRISPR gene-effect measurements.

DepMap releases are versioned. Reproduction should therefore use the release specified by the analysis configuration rather than automatically substituting the latest release.

### GSE176078

The breast cancer single-cell dataset is available from NCBI GEO:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
```

Accession:

```text
GSE176078
```

### gnomAD v4.1

gnomAD structural-variant resources are available from:

```text
https://gnomad.broadinstitute.org/downloads
```

### HGSVC2

HGSVC2 data are available from:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/
```

### COSMIC Cancer Gene Census

COSMIC Cancer Gene Census data are available from:

```text
https://cancer.sanger.ac.uk/census
```

COSMIC data remain subject to COSMIC access and licensing conditions.

## Installation

The Python environment is specified by:

```text
analysis/requirements.txt
analysis/environment.yml
```

Clone the repository:

```bash
git clone https://github.com/chykoushik/TNBC-SV-DNA-Repair.git
cd TNBC-SV-DNA-Repair/analysis
```

Using `pip`:

```bash
python -m venv .venv
```

Linux/macOS:

```bash
source .venv/bin/activate
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Alternatively, with Conda:

```bash
conda env create -f environment.yml
```

## Configuration

Public data-source definitions and dataset configuration are stored in:

```text
analysis/config/data_sources.yaml
```

Large raw datasets are not committed to GitHub.

Dataset locations should be configured for the local computing environment before executing analyses that require source data.

## Running the analysis

Run the analysis from:

```text
analysis/
```

The principal implementation is contained in:

```text
src/
```

Supporting scripts for data acquisition, verification, output construction, and package validation are located in:

```text
scripts/
```

The analysis code is separated into modules for:

- cohort construction
- expression processing
- covariate construction
- genomic outcomes
- candidate signatures
- statistical models
- null models
- patient-level analyses
- functional analyses
- single-cell analyses
- SCAN-B analyses
- figure generation

## Tests

Automated tests are provided under:

```text
analysis/tests/
```

Run the test suite from the `analysis/` directory:

```bash
pytest tests/
```

The tests cover major components of data acquisition, cohort construction, analysis stages, and final-output validation.

## Outputs

Analysis outputs are stored under:

```text
analysis/outputs/
```

Manuscript-facing figures and tables are retained in the final-output directory generated by the pipeline.

The outputs include results for:

- TCGA genomic HRD
- SCAN-B HRDetect
- SCAN-B scarHRD
- proliferation-adjusted models
- leave-one-gene-out analyses
- permutation analyses
- matched random-gene controls
- METABRIC
- GSE25066
- DepMap
- GSE176078
- population structural variation

Intermediate outputs required for traceability and reproducibility are retained where appropriate.

## Legacy analyses

The directory:

```text
legacy/
```

preserves the earlier notebook-based development of the project.

It contains:

```text
notebook1_setup_data.ipynb
notebook2_sv_gips.ipynb
notebook3_immune_eqtl.ipynb
notebook4_survival_validation.ipynb
notebook5_updated_brca.ipynb
notebook6_new_cohorts.ipynb
notebook7_clinical_cohort_validation.ipynb
results/
```

These notebooks document an earlier exploratory framework that included GIPS, immune analyses, methylation analyses, survival analyses, TCGA-OV analyses, and other analyses that are not part of the principal workflow described in this README.

They are retained for project provenance rather than as the recommended reproduction route.

The population structural-variation analysis developed during the earlier workflow contributed locus-level evidence retained in the study. Relevant SV outputs are therefore preserved with the reproducibility materials.

For reproduction of the reported analyses, use:

```text
analysis/
```

## Reproducibility considerations

The nine candidate genes are fixed rather than selected according to their observed associations with the study outcomes.

The principal score uses equal weighting rather than coefficients fitted to an HRD outcome.

TCGA TNBC status is defined from clinical receptor information rather than expression-derived receptor thresholds.

The same complete nine-gene score is used for the principal TCGA and SCAN-B analyses.

GSE25066 uses a reduced four-gene score because the complete candidate set is not available under the required platform mapping.

HRD, HRDetect, scarHRD, and copy-number burden are measured on different scales. Raw regression coefficients across these endpoints should therefore not be interpreted as directly comparable effect sizes.

## Code availability

Repository:

```text
https://github.com/chykoushik/TNBC-SV-DNA-Repair
```

The reproducible study pipeline is located in:

```text
analysis/
```

Earlier exploratory work is archived in:

```text
legacy/
```

## How to cite

If you use the code, analysis outputs, or methodology from this repository, please cite the associated article.

### Article

The bibliographic information below is a placeholder and should be updated after publication.

```text
Chowdhury K. Computational Analysis of Meiotic and Cohesin Gene Expression
and Homologous Recombination Deficiency in Triple-Negative Breast Cancer.
Manuscript submitted for publication.
```

### BibTeX

```bibtex
@article{Chowdhury_TNBC_Meiotic_Cohesin_HRD,
  author  = {Chowdhury, Koushik},
  title   = {Computational Analysis of Meiotic and Cohesin Gene Expression and Homologous Recombination Deficiency in Triple-Negative Breast Cancer},
  journal = {To be updated},
  year    = {2026},
  volume  = {To be updated},
  number  = {To be updated},
  pages   = {To be updated},
  doi     = {To be updated}
}
```


## License

Code in this repository is governed by the repository license.

The source datasets are not distributed under the repository license. Each dataset remains subject to the terms, licenses, and access requirements of its original provider.
