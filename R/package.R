#' FraminghamPACE: Methylation Pace of Aging Score Calculator
#'
#' This package provides functions to calculate FraminghamPACE estimates 
#' for methylation data generated from either the Illumina 450K array or 
#' the Illumina EPIC array.
#'
#' @keywords internal
#' @importFrom utils data
#' @importFrom stats aggregate
#' @importFrom preprocessCore normalize.quantiles.use.target
"_PACKAGE"

#' Framingham PACE Models Data
#'
#' A dataset containing the models and parameters for FraminghamPACEa analysis
#' with age and sex covariates.
#'
#' @format A list containing:
#' \describe{
#'   \item{model_names}{Names of available models}
#'   \item{model_probes}{Probes used in each model}
#'   \item{gold_standard_probes}{Background probes used for normalization}
#'   \item{gold_standard_probes}{Background probe means used for normalization}
#'   \item{model_weights}{Model coefficients}
#'   \item{model_means}{Mean values from training data}
#'   \item{model_intercept}{Model intercepts}
#'   \item{age_weight}{Age coefficient}
#'   \item{female_weight}{Sex coefficient}
#' }
#' @source Framingham Heart Study
"FraminghamPACEa_model"

#' No Covariates PACE Models Data
#'
#' A dataset containing the models and parameters for FraminghamPACEb analysis 
#' without age and sex covariates.
#'
#' @format A list containing:
#' \describe{
#'   \item{model_names}{Names of available models}
#'   \item{model_probes}{Probes used in each model}
#'   \item{gold_standard_probes}{Background probes used for normalization}
#'   \item{gold_standard_probes}{Background probe means used for normalization}
#'   \item{model_weights}{Model coefficients}
#'   \item{model_means}{Mean values from training data}
#'   \item{model_intercept}{Model intercepts}
#' }
#' @source Framingham Heart Study
"FraminghamPACEb_model"