#' See the list of probes needed for the FraminghamPACEa Model
#'
#' \code{getFraminghamPACEaProbes} returns a list containing the IDs of the probes needed for the FraminghamPACEa Model.
#' @name getFraminghamPACEaProbes
#' @rdname getFraminghamPACEaProbes
#' @param backgroundList (default: FALSE). If set to TRUE, it will return the full set of ~20,000 probes (~10,000 if MSA is selected) in the background set used for normalization.
#' @param model_type (default: "default"). Choose which model to use. Options are "default" for the original model or "MSA" for the alternative MSA model.
#' @return A character vector containing the requested CpGs.
#' @importFrom utils data
#' @export
#'
globalVariables("FraminghamPACEa_model")
globalVariables("FraminghamPACEa_model_MSA")

getFraminghamPACEaProbes = function(backgroundList = FALSE, model_type = "default") {
  
  # Load the appropriate model based on model_type
  if (model_type == "default") {
    data("FraminghamPACEa_model", package = "FraminghamPACE")
    model_data <- FraminghamPACEa_model
  } else if (model_type == "MSA") {
    data("FraminghamPACEa_model_MSA", package = "FraminghamPACE")
    model_data <- FraminghamPACEa_model_MSA
  } else {
    stop("Invalid model_type. Options are 'default' or 'MSA'.")
  }
  
  # Get the model name
  model_name <- model_data$model_names[1]
  
  if(backgroundList == FALSE) {
    return(model_data$model_probes[[model_name]])
  } else {
    return(model_data$gold_standard_probes[[model_name]])
  }
}






