# Preserve global variable declarations (not roxygen)
globalVariables(c("FraminghamPACEaProjector", "FraminghamPACEbProjector"))

#' Compute Framingham Pace of Aging Calculated from the Epigenome (FraminghamPACE) scores
#'
#' Computes FraminghamPACE scores for methylation arrays using one of four
#' supported model variants. Provide a probe-by-sample beta matrix and select a
#' variant via `model_type`.
#'
#' @param betas A numeric matrix of percent-methylation values with probe IDs as
#'   row names and samples as column names. Missing values should be `NA`.
#' @param model_type Character string selecting the model variant. One of
#'   "FraminghamPACE" (primary), "Alternate1", "Alternate2", or "Alternate3".
#'
#' @return A data frame of FraminghamPACE scores with sample names as row names
#'   and the model name as the single column.
#'
#' @examples
#' \dontrun{
#' FraminghamPACE <- getFraminghamPACE(betas, model_type = "FraminghamPACE")
#' }
#'
#' @export
getFraminghamPACE <- function(betas, model_type = "FraminghamPACE") {
  
  # Run the projector for the given model type
  if (model_type == "FraminghamPACE") {
    FraminghamPACE <- FraminghamPACEaProjector(betas, model_type = "RIDGE")
  } else if (model_type == "Alternate1") {
    FraminghamPACE <- FraminghamPACEbProjector(betas, model_type = "RIDGE")
  } else if (model_type == "Alternate2") {
    FraminghamPACE <- FraminghamPACEaProjector(betas, model_type = "ENET")
  } else if (model_type == "Alternate3") {
    FraminghamPACE <- FraminghamPACEbProjector(betas, model_type = "ENET")
  } else {
    stop("Invalid model_type. Options are 'FraminghamPACE', 'Alternate1', 'Alternate2', or 'Alternate3'.")
  }
  
  return(FraminghamPACE)
}
