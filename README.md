# FraminghamPACE
The FraminghamPACE R package allows users to calculate FraminghamPACE in their own DNA methylation dataset. Given a set of methylation beta values, where rows are probes and columns are subjects, this tool will produce FraminghamPACE. The algorithm works on data produced by the Illumina 450k, EPIC, EPIC V2, and MSA arrays.

There are two FraminghamPACE estimates, FraminghamPACEa and FraminghamPACEb. FraminghamPACEa was trained with age and sex included as additional covariates, whereas FraminghamPACEb was trained solely on CpGs. The methods perform very similarly across domains. Both are predicted with this package.

## Installation (via devtools):
To install and load FraminghamPACE, use the following code:
```r
devtools::install_github("will-marella/FraminghamPACE")

library(FraminghamPACE)
```

## Use
The primary functions in the FraminghamPACE package are the PACEProjector functions. There are two PACEProjectors: `FraminghamPACEaProjector()` and `FraminghamPACEbProjector()` for FraminghamPACEa and b, respectively. Helper windows are available for all functions via `?`, i.e. `?FraminghamPACEaProjector`

The PACEProjector functions calculate FraminghamPACE for the provided dataset. We recommend loading a beta matrix containing *all available* CpGs. For MSA array data, the "MSA" option must be selected via the `model_type` argument.

If it is preferred to load only the probes used for model fitting or background normalization, the CpGs are accessible via the `getFraminghamPACEaProbes()` and `getFraminghamPACEbProbes()` functions.

## FraminghamPACEaProjector():

#### Inputs:

###### betas
A numeric matrix containing the percent-methylation for each probe. Missing data should be 'NA's. The rows should be probes, with the probe ID as the row name, and the columns should be samples, with sample names as the column name.

###### age_sex_df
A numeric dataframe containing the age and sex values for each sample. The name of the age column should be 'Age' and the name of the sex column should be 'Female', where a 1 indicates the sample is female and 0 indicates the sample is male. The row names should indicate the sample, matching the column names of the 'betas' matrix.

###### proportionOfProbesRequired:
(default: 0.8). This value specifies the threshold for missing data (see description for more details on how missing data is handled). Note: if the function detects EPICv2 data based on row names that contain the replicate suffix, the proportion of probes required will automatically be lowered to 0.7.

###### model_type:
(default: "default"). Choose which model to use. Options are "default" for the original model or "MSA" for the alternative MSA model.

#### Output:
By default, a dataframe of FraminghamPACEa values with sample names as row names and the model name as column name. If return_normalized_betas=TRUE, a list containing both the pace value dataframe and the normalized_betas matrix.

## FraminghamPACEbProjector():

#### Inputs:

###### betas
A numeric matrix containing the percent-methylation for each probe. Missing data should be 'NA's. The rows should be probes, with the probe ID as the row name, and the columns should be samples, with sample names as the column name.

###### proportionOfProbesRequired:
(default: 0.8). This value specifies the threshold for missing data (see description for more details on how missing data is handled). Note: if the function detects EPICv2 data based on row names that contain the replicate suffix, the proportion of probes required will automatically be lowered to 0.7.

###### model_type:
(default: "default"). Choose which model to use. Options are "default" for the original model or "MSA" for the alternative MSA model.

#### Output:
By default, a dataframe of FraminghamPACEb values with sample names as row names and the model name as column name. If return_normalized_betas=TRUE, a list containing both the pace value dataframe and the normalized_betas matrix.

## Example Usage
```r
# Load in beta matrix
my_betas <- read.csv("preprocessed_beta_matrix_or_dataframe.csv")

# Load in age and sex dataframe (for FraminghamPACEa)
my_age_sex_df <- read.csv("demographics.csv")

# Run FraminghamPACEa
FraminghamPACEa <- FraminghamPACEaProjector(betas = my_betas, 
                                            age_sex_df = my_age_sex_df, 
                                            model_type = "default")

# Run FraminghamPACEb
FraminghamPACEb <- FraminghamPACEbProjector(betas = my_betas, 
                                            model_type = "default")
                                            
# NOTE: If your dataset is MSA, use model_type = "MSA"

```

## getFraminghamPACEaProbes():

#### Inputs:

###### backgroundList
(default: FALSE). If set to TRUE, it will return the full set of ~20,000 probes (~10,000 if MSA is selected) in the background set used for normalization.

###### model_type:
(default: "default"). Choose which model to use. Options are "default" for the original model or "MSA" for the alternative MSA model.

#### Output:
A vector of the CpGs used for FraminghamPACEa predictions. Or, if backgroundList is set to true, a vector of the CpGs in the background set used for normalization.

## getFraminghamPACEbProbes():

#### Inputs:

###### backgroundList
(default: FALSE). If set to TRUE, it will return the full set of ~20,000 probes (~10,000 if MSA is selected) in the background set used for normalization.

###### model_type:
(default: "default"). Choose which model to use. Options are "default" for the original model or "MSA" for the alternative MSA model.

#### Output:
A vector of the CpGs used for FraminghamPACEb predictions. Or, if backgroundList is set to true, a vector of the CpGs in the background set used for normalization.
