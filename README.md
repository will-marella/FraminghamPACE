# FraminghamPACE
The FraminghamPACE R package allows users to calculate FraminghamPACE in their own DNA methylation dataset. Given a set of methylation beta values, where rows are probes and columns are subjects, this tool will produce FraminghamPACE. The algorithm works on data produced by the Illumina 450k, EPIC, and EPIC v2 arrays.

There are two FraminghamPACE estimates, FraminghamPACEa and FraminghamPACEb. FraminghamPACEa was trained with age and sex included as additional covariates, whereas FraminghamPACEb was trained solely on CpGs. The methods perform very similarly across domains. Both are predicted with this package.

## Installation (via devtools):
To install and load FraminghamPACE, use the following code:
```r
devtools::install_github("will-marella/FraminghamPACE")

library(FraminghamPACE)
```

## Use
The primary function in the FraminghamPACE package is `getFraminghamPACE()`, which calculates FraminghamPACE for your dataset. Helper windows are available via `?`, e.g. `?getFraminghamPACE`.

If you prefer to load only the probes used for model fitting or background normalization, the CpGs are available via `getFraminghamPACEProbes()`. However, we highly recommend loading a beta matrix containing *all available* CpGs.

## Inputs

###### betas
A numeric matrix containing the percent-methylation for each probe. Missing data should be 'NA's. The rows should be probes, with the probe IDs as the row names, and the columns should be samples, with sample names as the column names.


###### model_type
Choose which model to use. Options are "FraminghamPACE" (primary), "Alternate1", "Alternate2", or "Alternate3".

###### sample_coverage_min
(default: 0.75) Minimum per-sample predictive coverage required to compute a score. Samples below this threshold are retained in the output but their score is set to `NA`. See Predictive coverage below.

###### return_coverage
(default: FALSE) If `TRUE`, return a list containing the per-sample coverage percentages and the threshold used.

###### return_normalized_betas
(default: FALSE) If `TRUE`, return the normalized beta matrix used for scoring.

#### Output
By default, a dataframe of FraminghamPACE values with sample names as row names and the model name as the column name.

## Predictive coverage

Predictive coverage quantifies how much of the model's signal is supported by observed (non-missing) probes in your data.

- **Per-sample coverage**: For each sample, we sum the model probe importances for probes that are present (non-`NA`) in that sample and multiply by 100. Because importances sum to 1 across scoring probes, per-sample coverage ranges from 0 to 100.
- **Dataset-level coverage (summary)**: Before imputation and after removing low-coverage samples, we report a dataset-level predictive coverage. This is the importance-weighted average fraction of non-missing values per scoring probe across the kept samples, expressed as a percentage, and labeled as **Good** (>90%), **Moderate** (>75%), or **Poor** (≤75%). We also report the count of available background and scoring probes detected in your input.

### Per-sample coverage and filtering

- **Thresholding**: Samples with per-sample coverage below `sample_coverage_min` (default 0.75 → 75%) have their scores set to `NA` and are excluded from normalization and model fitting. A message indicates how many samples remain above the threshold.
- **Tuning**: You can relax or tighten this filter by changing `sample_coverage_min`, e.g., `sample_coverage_min = 0.70`.

### Returning coverage and normalized betas

Set `return_coverage = TRUE` to receive a list with `coverage$sample` (named numeric vector of per-sample coverage percentages) and `coverage$threshold`. Set `return_normalized_betas = TRUE` to also receive the quantile-normalized beta matrix used for scoring.

### Replicates

If probe row names contain EPIC v2 replicate suffixes (e.g., "TC" or "BC"), replicates are automatically averaged prior to scoring.

## Example Usage
```r
# Load in beta matrix
my_betas <- read.csv("preprocessed_beta_matrix_or_dataframe.csv")

# Compute FraminghamPACE (primary model)
FraminghamPACE <- getFraminghamPACE(betas = my_betas, model_type = "FraminghamPACE")

```

### Example with coverage outputs
```r
res <- getFraminghamPACE(
  betas = my_betas,
  model_type = "FraminghamPACE",
  sample_coverage_min = 0.70,
  return_coverage = TRUE,
  return_normalized_betas = FALSE
)

# Scores (data.frame with one column)
res$pace_score

# Per-sample predictive coverage (%)
head(res$coverage$sample)

# Threshold used
res$coverage$threshold
```
## getFraminghamPACEProbes()

#### Inputs

###### backgroundList
(default: FALSE) If set to TRUE, it will return the full set of ~15,000 probes in the background set used for normalization.

###### model_type
Choose which model to use. Options are "FraminghamPACE" (primary), "Alternate1", "Alternate2", or "Alternate3".

#### Output
A vector of the CpGs used for FraminghamPACE predictions. Or, if backgroundList is set to TRUE, a vector of the CpGs in the background set used for normalization.
