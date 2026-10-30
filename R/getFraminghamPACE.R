# Preserve global variable declarations (not roxygen)
globalVariables(c("FraminghamPACEaProjector"))

#' Compute Framingham Pace of Aging Calculated from the Epigenome (FraminghamPACE) scores
#'
#' Computes FraminghamPACE scores for methylation arrays using one of two
#' supported WTM model variants. Provide your beta matrix to calculate FraminghamPACE.
#'
#' @param betas A numeric matrix of percent-methylation values with probe IDs as
#'   row names and samples as column names. Missing values should be `NA`.
#' @param model_type Character string selecting the model variant. One of
#'   "Ridge" or "ElasticNet" (default).
#' @param ... Optional arguments you can pass to control filtering and outputs:
#'   - `sample_coverage_min` (numeric; default 0.75): Minimum per-sample predictive
#'     coverage required to compute a score. Samples below the threshold are
#'     excluded from normalization/scoring and returned as `NA` in the score output.
#'   - `return_coverage` (logical; default FALSE): If `TRUE`, also return a list
#'     containing per-sample coverage percentages and the threshold used.
#'   - `return_normalized_betas` (logical; default FALSE): If `TRUE`, also return
#'     the quantile-normalized beta matrix used for scoring.
#'   - `timing` (logical; default FALSE): If `TRUE`, print timing for major
#'     processing stages (validation, replicate averaging, background subsetting,
#'     imputation, normalization, scoring). Has no effect on results.
#'
#' @details
#' Predictive coverage reflects the fraction of model signal supported by observed
#' (non-missing) probes. Coverage is computed prior to any imputation. Samples with
#' coverage below `sample_coverage_min` have scores set to `NA` and are excluded
#' from normalization and model fitting. Probes entirely absent in the input are
#' added at gold-standard means. For probes that are present but have missing
#' values, missing entries are imputed with cohort means when enough values are
#' observed; otherwise, the probe is replaced with model means.
#'
#' Replicate probes are automatically averaged prior to scoring.
#'
#' @return
#' - If `return_coverage = FALSE` and `return_normalized_betas = FALSE`:
#'   a data frame of FraminghamPACE scores with sample names as row names and the
#'   selected model type (`"ElasticNet"` or `"Ridge"`) as the single column name.
#' - Otherwise: a list with elements
#'   - `pace_score`: data frame of scores (as above)
#'   - `normalized_betas`: matrix or `NULL` depending on `return_normalized_betas`
#'   - `coverage`: list with `sample` (named numeric vector of per-sample coverage
#'      percentages) and `threshold` (numeric `sample_coverage_min`)
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' scores <- getFraminghamPACE(betas)
#'
#' # With a different model, coverage outputs, and a relaxed per-sample coverage threshold
#' res <- getFraminghamPACE(
#'   betas,
#'   model_type = "Ridge",
#'   sample_coverage_min = 0.70,
#'   return_coverage = TRUE,
#'   return_normalized_betas = FALSE
#' )
#' res$pace_score            # data.frame of scores
#' head(res$coverage$sample) # per-sample predictive coverage (%)
#' res$coverage$threshold    # threshold used
#' }
#'
#' @export
getFraminghamPACE <- function(betas, model_type = "ElasticNet", ...) {
  
  # Run the projector for the given model type
  if (model_type == "Ridge") {
    FraminghamPACE <- FraminghamPACEaProjector(betas, model_type = "WTM_RIDGE", ...)
  } else if (model_type == "ElasticNet") {
    FraminghamPACE <- FraminghamPACEaProjector(betas, model_type = "WTM_ENET", ...)
  } else {
    stop("Invalid model_type. Options are 'Ridge' or 'ElasticNet'.")
  }
  
  return(FraminghamPACE)
}
