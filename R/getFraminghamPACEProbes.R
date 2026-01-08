#' List probes used by FraminghamPACE models
#'
#' Returns the set of probe IDs required by the selected model variant. Use
#' `backgroundList = TRUE` to return the full background set used during
#' normalization; otherwise, returns the scoring probes only.
#'
#' @param backgroundList Logical (default: FALSE). If TRUE, return background
#'   probes; if FALSE, return scoring probes.
#' @param model_type Character string selecting the model variant.
#'   Options are 'FraminghamPACE', 'Alternate1', 'Alternate2', 'Alternate3', 'WTM_RIDGE', or 'WTM_ENET'.
#'
#' @return A character vector of probe IDs.
#' @examples
#' \dontrun{
#' probes <- getFraminghamPACEProbes(backgroundList = FALSE, model_type = "FraminghamPACE")
#' }
#' @importFrom utils data
#' @export
getFraminghamPACEProbes = function(backgroundList = FALSE, model_type = "FraminghamPACE") {
  if (model_type == "FraminghamPACE") {
    data("FraminghamPACEa_RIDGE", package = "FraminghamPACE")
    model_data <- FraminghamPACEa_RIDGE
  } else if (model_type == "Alternate1") {
    data("FraminghamPACEb_RIDGE", package = "FraminghamPACE")
    model_data <- FraminghamPACEb_RIDGE
  } else if (model_type == "Alternate2") {
    data("FraminghamPACEa_ENET", package = "FraminghamPACE")
    model_data <- FraminghamPACEa_ENET
  } else if (model_type == "WTM_RIDGE") {
    data("WTM_FraminghamPACEa_RIDGE", package = "FraminghamPACE")
    model_data <- WTM_FraminghamPACEa_RIDGE
  } else if (model_type == "WTM_ENET") {
    data("WTM_FraminghamPACEa_ENET", package = "FraminghamPACE")
    model_data <- WTM_FraminghamPACEa_ENET
  } else if (model_type == "Alternate3") {
    data("FraminghamPACEb_ENET", package = "FraminghamPACE")
    model_data <- FraminghamPACEb_ENET
  } else {
    stop("Invalid model_type. Options are 'FraminghamPACE', 'Alternate1', 'Alternate2', 'Alternate3', 'WTM_RIDGE', or 'WTM_ENET'.")
  }
  model_name <- model_data$model_names[1]
  if (isFALSE(backgroundList)) {
    return(model_data$model_probes[[model_name]])
  } else {
    return(model_data$gold_standard_probes[[model_name]])
  }
}
