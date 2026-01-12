#' List probes used by FraminghamPACE models
#'
#' Returns the set of probe IDs required by the selected model variant. Use
#' `backgroundList = TRUE` to return the full background set used during
#' normalization; otherwise, returns the scoring probes only.
#'
#' @param backgroundList Logical (default: FALSE). If TRUE, return background
#'   probes; if FALSE, return scoring probes.
#' @param model_type Character string selecting the model variant.
#'   Options are 'Ridge' or 'ElasticNet'.
#'
#' @return A character vector of probe IDs.
#' @examples
#' \dontrun{
#' probes <- getFraminghamPACEProbes(backgroundList = FALSE, model_type = "Ridge")
#' }
#' @importFrom utils data
#' @export
getFraminghamPACEProbes = function(backgroundList = FALSE, model_type = "ElasticNet") {
  if (model_type == "Ridge") {
    data("WTM_FraminghamPACEa_RIDGE", package = "FraminghamPACE")
    model_data <- WTM_FraminghamPACEa_RIDGE
  } else if (model_type == "ElasticNet") {
    data("WTM_FraminghamPACEa_ENET", package = "FraminghamPACE")
    model_data <- WTM_FraminghamPACEa_ENET
  } else {
    stop("Invalid model_type. Options are 'Ridge' or 'ElasticNet'.")
  }
  model_name <- model_data$model_names[1]
  if (isFALSE(backgroundList)) {
    return(model_data$model_probes[[model_name]])
  } else {
    return(model_data$gold_standard_probes[[model_name]])
  }
}
