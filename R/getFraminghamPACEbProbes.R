#' See the list of probes needed for the FraminghamPACEa Model
#'
#' \code{getFraminghamPACEbProbes} returns a list containing the IDs of the probes needed for the FraminghamPACEb Model.
#' @name getFraminghamPACEaProbes
#' @rdname getFraminghamPACEaProbes
#' @param backgroundList (default: FALSE).  If set to true, it will return the full set of 20,000 probes in the background set used for normalization
#' @return A list of string vectors.  Each element in the list is a different model, and the character vector contains the IDs needed for that particular model.
#' @importFrom utils data
globalVariables("NoCov_mPACE_models")

getFraminghamPACEbProbes = function(backgroundList=F) {
  data("FraminghamPACEb", package = "FraminghamPACE")
  if( backgroundList == FALSE ) {
    NoCov_mPACE_models$model_probes
  } else {
    NoCov_mPACE_models$gold_standard_probes
  }
}






