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

###### proportionOfProbesRequired
(default: 0.8) The proportion of probes required to have a non-missing value for FraminghamPACE for a sample, as well as to determine if we can impute the mean from the current cohort. Note: if the function detects EPIC v2 data based on row names that contain the replicate suffix, the proportion of probes required will automatically be lowered to 0.7.

###### model_type
Choose which model to use. Options are "FraminghamPACE" (primary), "Alternate1", "Alternate2", or "Alternate3".

#### Output
By default, a dataframe of FraminghamPACE values with sample names as row names and the model name as the column name.

## Example Usage
```r
# Load in beta matrix
my_betas <- read.csv("preprocessed_beta_matrix_or_dataframe.csv")

# Compute FraminghamPACE (primary model)
FraminghamPACE <- getFraminghamPACE(betas = my_betas, model_type = "FraminghamPACE")

```
## getFraminghamPACEProbes()

#### Inputs

###### backgroundList
(default: FALSE) If set to TRUE, it will return the full set of ~15,000 probes in the background set used for normalization.

###### model_type
Choose which model to use. Options are "FraminghamPACE" (primary), "Alternate1", "Alternate2", or "Alternate3".

#### Output
A vector of the CpGs used for FraminghamPACE predictions. Or, if backgroundList is set to TRUE, a vector of the CpGs in the background set used for normalization.
