#' See the list of probes needed for the FraminghamPACEa Model
#'
#' \code{getFraminghamPACEaProbes} returns a list containing the IDs of the probes needed for the FraminghamPACEa Model.
#' @name getFraminghamPACEaProbes
#' @rdname getFraminghamPACEaProbes
#' @param backgroundList (default: FALSE).  If set to true, it will return the full set of 20,000 probes in the background set used for normalization
#' @return A list of string vectors.  Each element in the list is a different model, and the character vector contains the IDs needed for that particular model.
#' @importFrom utils data
#'
globalVariables("FraminghamPACEa_model")

getFraminghamPACEaProbes = function(backgroundList=F) {
  data("FraminghamPACEa_model", package = "FraminghamPACE")
  if( backgroundList == FALSE ) {
    FraminghamPACEa_model$model_probes
  } else {
    FraminghamPACEa_model$gold_standard_probes
  }
}






