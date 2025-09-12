## FraminghamPACE
FraminghamPACE is an R package for computing the Framingham Pace of Aging Calculated from the Epigenome (FraminghamPACE) from DNA methylation arrays. Provide a matrix of methylation beta values, with rows as CpG probes and columns as samples, and the package returns PACE estimates. The algorithm supports Illumina 450k, EPIC, EPIC v2, and MSA arrays.

The package can estimate four closely related model variants of FraminghamPACE, which were trained with slightly different feature sets and modelling strategies. The primary variance, FraminghamPACE, will be the selection for the majority of use cases, whereas the 'Alternate' versions 1, 2, and 3, ... See () for more information.

## Installation (via devtools)
To install and load FraminghamPACE:
```r
devtools::install_github("will-marella/FraminghamPACE")
library(FraminghamPACE)
```

## Use
The primary function is `getFraminghamPACE()`, which calculates FraminghamPACE for your dataset. Help is available via `?getFraminghamPACE`.

To retrieve the CpGs used by the models (either scoring probes or the larger background set for normalization), use `getFraminghamPACEProbes()`. For best estimates when calculating FraminghamPACE, we recommend providing a beta matrix containing all available CpGs.

---

## getFraminghamPACE()

### Arguments
- **betas (input matrix)**: Numeric matrix of methylation beta values representing the proportion methylated (0–1) at each CpG. Rows are CpG probe IDs (`cg...`) and columns are samples. Missing values should be `NA`.
  - Columns (sample names) should be unique. Duplicate sample names are permitted in R but can lead to ambiguous outputs and downstream merges; we strongly recommend making them unique.
  - Data frames are accepted and converted to a numeric matrix internally.
  - Replicate probes are handled automatically (see Replicates below).

- **model_type**: Which variant to compute. One of `"FraminghamPACE"` (default), `"Alternate1"`, `"Alternate2"`, or `"Alternate3"`.

- **sample_coverage_min**: Minimum per-sample predictive coverage required to compute a score. Samples below this threshold are kept in the result but their score is set to `NA`. This filter is applied before normalization and model scoring; see Predictive coverage below for details. Default: `0.75` (i.e., 75%).

- **return_coverage**: If `TRUE`, also return per-sample predictive coverage percentages and the threshold used. Default: `FALSE`.

- **return_normalized_betas**: If `TRUE`, also return the quantile-normalized beta matrix used for scoring. Default: `FALSE`.

### Output
- **Default**: A data frame of PACE scores with sample names as row names and the model name as the single column.
- **If `return_coverage = TRUE`**: The function returns a list with:
  - `pace_score`: as above
  - `coverage$sample`: named numeric vector of per-sample coverage percentages (0–100)
  - `coverage$threshold`: numeric value of `sample_coverage_min`
- **If `return_normalized_betas = TRUE`**: The list includes:
  - `normalized_betas`: the quantile-normalized beta matrix used in scoring (or `NULL` if `FALSE`)
- **If both are TRUE**: You receive both coverage and normalized betas in the list alongside `pace_score`.

### Predictive coverage
Model coefficients (weights) determine each probe’s contribution to the overall PACE estimate. We summarize each probe’s contribution as an importance score based on the coefficient estimate and variance of the probe in the training sample. Predictive coverage then quantifies, for the data you supplied, how much of the model’s weighted signal is actually supported by observed (non-missing) probes.

Per-sample predictive coverage is the sum of probe importance for probes that are present in that sample, multiplied by 100. Because probe importance is normalized to sum to 1 across scoring probes, per-sample coverage ranges from 0 to 100.

At the dataset level (before any imputation and after removing samples below `sample_coverage_min`), we report an importance-weighted percentage summarizing how well the scoring probes are covered across the kept samples. For orientation, we label this as Good (>90%), Moderate (>75%), or Poor (≤75%). We also report how many background and scoring probes were detected in your input.

### Per-sample filtering
Samples with coverage below `sample_coverage_min` have their scores set to `NA` and are excluded from normalization and model fitting. You can relax or tighten this threshold (e.g., `sample_coverage_min = 0.70`). Messages report how many samples remain and a brief predictive coverage summary.

### Returning coverage and normalized betas
If you set `return_coverage = TRUE`, the output includes `coverage$sample` (per-sample coverage percentages) and `coverage$threshold`. If you set `return_normalized_betas = TRUE`, the output includes `normalized_betas` (the quantile-normalized beta matrix used for scoring). If both are `TRUE`, both are returned.

### Replicates
In newer platforms (e.g., EPICv2), replicate probes can appear, with a suffix that may include strand (T/B), original or bisulfite-converted strand (C/O), probe type (I/II), and replicate index (1–10). For the calculation of FraminghamPACE these replicates are automatically averaged prior to scoring. Replicates are identified by a shared core probe ID with replicate-specific suffixes; these are collapsed to the core ID before normalization and scoring.

### Example
```r
# Load beta matrix (proportion methylated, rows = probes, cols = samples)
my_betas <- read.csv("preprocessed_beta_matrix_or_dataframe.csv")

# Compute FraminghamPACE (default model)
scores <- getFraminghamPACE(betas = my_betas)

# Example with an alternate model, relaxed coverage, and coverage output
res <- getFraminghamPACE(
  betas = my_betas,
  model_type = "Alternate1",
  sample_coverage_min = 0.70,
  return_coverage = TRUE,
  return_normalized_betas = FALSE
)
res$pace_score            # data.frame of scores
head(res$coverage$sample) # per-sample predictive coverage (%)
res$coverage$threshold    # threshold used
```

---

## getFraminghamPACEProbes()

### Arguments
- **backgroundList**: Logical. If `TRUE`, return the full set of background probes used for normalization (~15k CpGs). If `FALSE` (default), return only the scoring probes used for PACE predictions.
- **model_type**: One of `"FraminghamPACE"` (primary), `"Alternate1"`, `"Alternate2"`, or `"Alternate3"`.

### Output
- A character vector of CpG probe IDs. When `backgroundList = TRUE`, returns the background set; otherwise, returns the scoring set for the selected model.
