# Preserve global variable declarations (not roxygen)
globalVariables(c("FraminghamPACEaProjector", "FraminghamPACEbProjector"))

#' Orchestrate FraminghamPACE projectors
#'
#' Runs the appropriate projector based on `model_type` and returns a
#' FraminghamPACE score data frame.
#'
#' @param betas A numeric matrix of percent-methylation values with probe IDs as
#'   row names and samples as column names. Missing values should be NA.
#' @param model_type Character string; one of "default", "Alternate1",
#'   "Alternate2", or "Alternate3".
#'
#' @return A data frame of FraminghamPACE scores with samples as row names and
#'   the model name as the single column.
#'
#' @examples
#' \dontrun{
#' FraminghamPACE <- getFraminghamPACE(betas, model_type = "default")
#' }
#'
#' @export
getFraminghamPACE <- function(betas, model_type = "default") {
  
  # Run the projector for the given model type
  if (model_type == "default") {
    FraminghamPACE <- FraminghamPACEaProjector(betas, model_type = "RIDGE")
  } else if (model_type == "Alternate1") {
    FraminghamPACE <- FraminghamPACEbProjector(betas, model_type = "RIDGE")
  } else if (model_type == "Alternate2") {
    FraminghamPACE <- FraminghamPACEaProjector(betas, model_type = "ENET")
  } else if (model_type == "Alternate3") {
    FraminghamPACE <- FraminghamPACEbProjector(betas, model_type = "ENET")
  } else {
    stop("Invalid model_type. Options are 'default', 'Alternate1', 'Alternate2', or 'Alternate3'.")
  }
  
  return(FraminghamPACE)
}






