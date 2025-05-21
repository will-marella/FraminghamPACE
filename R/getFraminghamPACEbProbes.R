#' See the list of probes needed for the FraminghamPACEb Model
#'
#' \code{getFraminghamPACEbProbes} returns a list containing the IDs of the probes needed for the FraminghamPACEb Model.
#' @name getFraminghamPACEbProbes
#' @rdname getFraminghamPACEbProbes
#' @param backgroundList (default: FALSE). If set to TRUE, it will return the full set of ~20,000 probes (~10,000 if MSA is selected) in the background set used for normalization.
#' @param model_type (default: "default"). Choose which model to use. Options are "default" for the original model or "MSA" for the alternative MSA model.
#' @return A list of string vectors.  Each element in the list is a different model, and the character vector contains the IDs needed for that particular model.
#' @importFrom utils data
#' @export
globalVariables("FraminghamPACEb_model")
globalVariables("FraminghamPACEb_model_MSA")

getFraminghamPACEbProbes = function(backgroundList = FALSE, model_type = "default") {
  
  # Load the appropriate model based on model_type
  if (model_type == "default") {
    data("FraminghamPACEb_model", package = "FraminghamPACE")
    model_data <- FraminghamPACEb_model
  } else if (model_type == "MSA") {
    data("FraminghamPACEb_model_MSA", package = "FraminghamPACE")
    model_data <- FraminghamPACEb_model_MSA
  } else {
    stop("Invalid model_type. Options are 'default' or 'MSA'.")
  }
  
  # Get the model name
  model_name <- model_data$model_names[1]
  
  if( backgroundList == FALSE ) {
    return(model_data$model_probes[[model_name]])
  } else {
    return(model_data$gold_standard_probes[[model_name]])
  }
}






