# Epigenetic Age Prediction from Saliva DNA Methylation

A 19-CpG epigenetic clock for chronological age prediction in saliva-derived samples, developed as a final-year BSc Forensic Science dissertation at Northumbria University.

**Author:** Connor Wright
**Supervisor:** Dr Ed Schwalbe
**Module:** AP0600 — Research Project (2025/26)

---

## Overview

This project develops a saliva-optimised epigenetic age prediction model from publicly available DNA methylation data. A Lasso GLMnet regression was trained on 199 integrated saliva samples and benchmarked against three established epigenetic clocks. The final 19-CpG signature is intended as the first step toward a forensically deployable assay.

## Key Results

| Model | R² (test) | RMSE (test, years) | N CpGs |
|---|---|---|---|
| **Final Lasso (this work)** | **0.951** | **6.96** | **19** |
| Zhang Elastic Net (Zhang et al., 2019) | 0.933 | 6.25 | 514 |
| Zhang BLUP (Zhang et al., 2019) | 0.915 | 7.16 | ~320k |
| Horvath multi-tissue (Horvath, 2013) | 0.876 | 8.93 | 353 |

The final model uses lasso regression with λ ≈ 4.033 and tolerates up to ~20–25% missing data before performance deteriorates markedly.

## Dataset

199 saliva samples integrated from five GEO datasets, spanning ages 5–91 years:

| GEO Accession | n | Notes |
|---|---|---|
| GSE138279 | 65 | Paediatric cohort |
| GSE99091 | 57 | One duplicate removed (originally 58) |
| GSE92767 | 54 | |
| GSE59509 | 12 | |
| GSE111165 | 12 | |

All samples profiled on Illumina Infinium HumanMethylation450K or MethylationEPIC BeadChip arrays. A 70/30 age-stratified train/test split was used (140 training, 59 test).

## Repository Structure

```
.
├── README.md
├── .gitignore
├── R/
│   ├── 01_pca_analysis.R              # PCA: Figures 1–3
│   ├── 02_age_associated_cpgs.R       # 8 CpG correlations: Figure 4
│   ├── 03_published_clocks.R          # Horvath/EN/BLUP benchmarking: Figure 5
│   ├── 04_supervisor_todo_FINAL.R     # GLMnet + final model + simulation: Figures 6–8
│   └── 05_primer_design.R             # PrimerSuite preparation for 19-CpG signature
├── outputs/
│   ├── final_signature_CpGs_CORRECTED.csv
│   ├── glmnet_model_comparison_CORRECTED.csv
│   ├── missing_data_summary_CORRECTED.csv
│   └── cpg_primers_hg38.fasta
└── docs/
    └── dissertation_AP0600.docx       # Final dissertation document
```

## Data Availability

Beta matrices and phenotype RDS files are **not stored in this repository** due to file size (>100 MB). They are:

- Available locally in OneDrive: `Research Project/Methylation Data from Drive/Connor2/Main Files/16th Feb Correct files/`
- Reproducible from the original GEO accessions listed above

Required input files:
- `train_beta_CORRECTED.RDS`
- `test_beta_CORRECTED.RDS`
- `train_pheno_CORRECTED.RDS`
- `test_pheno_CORRECTED.RDS`

Drop these into the working directory before running any script.

## Requirements

- **R** ≥ 4.4.1
- **Core packages:** `glmnet`, `methylclock`, `impute`, `ggplot2`, `gridExtra`, `cowplot`, `dplyr`, `tidyr`
- **Annotation:** `IlluminaHumanMethylation450kanno.ilmn12.hg19`, `BSgenome.Hsapiens.UCSC.hg38`, `rtracklayer`

Install everything with:

```r
install.packages(c("glmnet", "ggplot2", "gridExtra", "cowplot", "dplyr", "tidyr"))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("methylclock", "impute",
                       "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       "BSgenome.Hsapiens.UCSC.hg38", "rtracklayer"))
```

> **Note:** `methylclock`'s `DNAmAge()` requires `dplyr` ≥ 1.2.0 and a beta matrix with CpG IDs in the first column rather than as row names. The clock benchmarking script handles this restructuring automatically.

## Usage

Run scripts in numerical order from the `R/` directory:

```r
setwd("path/to/this/repo")
source("R/01_pca_analysis.R")
source("R/02_age_associated_cpgs.R")
source("R/03_published_clocks.R")
source("R/04_supervisor_todo_FINAL.R")   # ~30 min with n_iterations = 1000
```

The simulation in script 04 defaults to 100 iterations for testing; set `n_iterations <- 1000` for the final run.

## Methodology Summary

1. **Data integration** — five GEO datasets harmonised; one GSE99091 duplicate removed
2. **Quality control** — k-NN imputation for missing values; sample ID corrections cross-referenced against GEO2R
3. **Feature selection** — top 10,000 CpGs by absolute correlation with age in the training set
4. **Benchmarking** — Horvath, Zhang EN, and Zhang BLUP clocks applied via `methylclock::DNAmAge()`
5. **Model development** — penalised regression (ridge α=0, elastic net α=0.5, lasso α=1) with 10-fold CV; lasso restricted to λ values yielding 5–20 non-zero CpGs
6. **Robustness simulation** — 1,000 iterations at each of 5–50% missing data (5% increments), using Ed's custom k-NN imputation with the distance matrix recomputed from degraded data per iteration
7. **Assay translation** — bisulfite-conversion PCR primers designed via PrimerSuite for each of the 19 CpGs, with sequences extracted from hg38 via `BSgenome.Hsapiens.UCSC.hg38`

## Citation

If referencing this work prior to publication:

> Wright, C. (2026). *Development of a saliva-derived 19-CpG epigenetic clock for forensic age estimation.* BSc Dissertation, Northumbria University.

## License

This work is shared for academic review only. Please do not redistribute without permission until the dissertation has been marked.

## Acknowledgements

- Dr Ed Schwalbe (project supervisor) — for the guidance across the entire process
- The original investigators behind GSE138279, GSE99091, GSE92767, GSE59509, and GSE111165 for making their data publicly available
- Northumbria University Department of Applied Sciences
