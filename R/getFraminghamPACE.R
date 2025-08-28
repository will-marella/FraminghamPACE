#' @export

globalVariables(c("FraminghamPACEaProjector", "FraminghamPACEbProjector"))

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






