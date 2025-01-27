#' Generate the Framingham Methylation Pace of Aging Scores!
#'
#' \code{FraminghamPACEbProjector} returns the FraminghamPACEb of Aging Methylation Scores
#'
#' @name FraminghamPACEbProjector
#' @rdname FraminghamPACEbProjector
#' @param betas A numeric matrix containing the percent-methylation for each probe. Missing data should be 'NA's. Rows are probes (rownames), and columns are samples (colnames).
#' @param proportionOfProbesRequired (default: 0.8). Threshold for missing data. If EPICv2 data is detected, this is automatically lowered to 0.7.
#' @return A list of mPACE values, one for each mPACE model. Each element is a named vector of scores, equivalent to the pace of aging (in years) over a 1-year period.
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @importFrom utils data
#' @details Missing data is handled by imputing means for probes with missing values and setting NA for samples missing too much data. Probes missing entirely are imputed with the mean from the training dataset.
#' @export
globalVariables("NoCov_mPACE_models")


FraminghamPACEbProjector = function( betas, proportionOfProbesRequired=0.8 ) {
  requireNamespace("preprocessCore")
  data("FraminghamPACEb", package = "FraminghamPACE")

  if( any(grepl("TC|BC", rownames(betas))) )
  {
    # Print message to user
    print("This looks like EPICv2 array data. If EPICv2, FraminghamPACE will lower the proportionOfProbesRequired to 0.7 and proceed with missing probes present on 450k and EPICv1. Averaging the opposite strand replicates may take a bit more time.")

    # Proceed with averaging: Code taken from ENmix package.
    cgid = sapply(strsplit(rownames(betas), split = "_"),  unlist)[1, ]
    dupcg = unique(cgid[duplicated(cgid)])
    betas2 = betas[cgid %in% dupcg, ]
    cid = sapply(strsplit(rownames(betas2), split = "_"),   unlist)[1, ]
    betas2 = aggregate(betas2, by = list(cid), FUN = function(x) mean(x, na.rm = TRUE))
    rownames(betas2) = betas2[, 1]
    betas2 = as.matrix(betas2[, -1])
    betas = betas[!(cgid %in% dupcg), ]
    rownames(betas) = sapply(strsplit(rownames(betas), split = "_"), unlist)[1, ]
    betas <-rbind(betas, betas2)

    # Set proportionOfProbesRequired to 0.7 if condition is TRUE
    proportionOfProbesRequired <- 0.7
  }

  # loop through models
  model_results <- lapply(NoCov_mPACE_models$model_names, function(model_name) {
    # make sure it has been converted to a matrix
    if( !is.numeric(as.matrix(betas)) ) { stop("betas matrix/data.frame is not numeric!") }
    probeOverlap <- length(which(rownames(betas) %in% NoCov_mPACE_models$model_probes[[model_name]])) / length(NoCov_mPACE_models$model_probes[[model_name]])
    probeOverlap_background <- length(which(rownames(betas) %in% NoCov_mPACE_models$gold_standard_probes[[model_name]])) / length(NoCov_mPACE_models$gold_standard_probes[[model_name]])
    # make sure enough of the probes are present in the data file
    if( probeOverlap < proportionOfProbesRequired | probeOverlap_background < proportionOfProbesRequired ) {
      result <- rep(NA, ncol(betas))
      names(result) <- colnames(betas)
      
      list(
        pace_score = result,
        normalized_betas = NULL
      )
      
    } else {
      # Work with a numeric matrix of betas
      betas.mat <- as.matrix(betas[which(rownames(betas) %in% NoCov_mPACE_models$gold_standard_probes[[model_name]]),])
      # If probes don't exist, we'll add them as rows of values based on their mean in the gold standard dataset
      probesNotInMatrix <- NoCov_mPACE_models$gold_standard_probes[[model_name]][which(NoCov_mPACE_models$gold_standard_probes[[model_name]] %in% rownames(betas.mat) == F)]
      if( length(probesNotInMatrix) > 0 ) {
        for( probe in probesNotInMatrix ) {
          tmp.mat <- matrix(0, nrow=1, ncol=ncol(betas.mat))
          rownames(tmp.mat) <- probe
          colnames(tmp.mat) <- colnames(betas.mat)
          tmp.mat[probe,] <- rep(NoCov_mPACE_models$gold_standard_means[[model_name]][probe], ncol(tmp.mat))
          betas.mat <- rbind(betas.mat, tmp.mat)
        }
      }

      # Identify samples with too many missing probes and remove them from the matrix
      samplesToRemove <- colnames(betas.mat)[which(apply(betas.mat, 2, function(x) { 1 - ( length(which(is.na(x))) / length(x) ) < proportionOfProbesRequired}))]
      if( length(samplesToRemove) > 0 ) {
        betas.mat <- betas.mat[,-which(colnames(betas.mat) %in% samplesToRemove)]
      }
      if(ncol(betas.mat) > 0) {
        # Identify missingness on a probe level
        pctValuesPresent <- apply( betas.mat, 1, function(x) { 1 - (length(which(is.na(x))) / length(x)) } )
        # If they're missing values, but less than the proportion required, we impute to the cohort mean
        probesToAdjust <- which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired)
        if( length(probesToAdjust) > 0 ) {
          if( length(probesToAdjust) > 1 ) {
            betas.mat[probesToAdjust,] <- t(apply( betas.mat[probesToAdjust,], 1 , function(x) {
              x[is.na(x)] = mean( x, na.rm = TRUE )
              x
            }))
          } else {
            betas.mat[probesToAdjust,which(is.na(betas.mat[probesToAdjust,]))] <- mean(betas.mat[probesToAdjust,], na.rm=T)
          }
        }
        # If they're missing too many values, everyones value gets replaced with the mean from the Framingham cohort
        if( length(which(pctValuesPresent < proportionOfProbesRequired)) > 0 ) {
          probesToReplaceWithMean <- rownames(betas.mat)[which(pctValuesPresent < proportionOfProbesRequired)]
          for( probe in probesToReplaceWithMean ) {
            betas.mat[probe,] <- rep(NoCov_mPACE_models$model_means[[model_name]][probe], ncol(betas.mat))
          }
        }

        # Normalize the matrix to the gold standard dataset
        betas.norm <- preprocessCore::normalize.quantiles.use.target(betas.mat, target=NoCov_mPACE_models$gold_standard_means[[model_name]])
        rownames(betas.norm) <- rownames(betas.mat)
        colnames(betas.norm) <- colnames(betas.mat)
        # Calculate score:
        score = NoCov_mPACE_models$model_intercept[[model_name]] + rowSums(t(betas.norm[NoCov_mPACE_models$model_probes[[model_name]],]) %*% diag(NoCov_mPACE_models$model_weights[[model_name]]))
        names(score) <- colnames(betas.norm)
        if( length(samplesToRemove) > 0 ) {
          score.tmp <- rep(NA, length(samplesToRemove))
          names(score.tmp) <- samplesToRemove
          score <- c(score, score.tmp)
        }
        score <- score[colnames(betas)]
        
        list(
          pace_score = score,
          normalized_betas = betas.norm
        )
        
      } else {
        result <- rep(NA, ncol(betas.mat))
        names(result) <- colnames(betas.mat)
        
        list(
          pace_score = result,
          normalized_betas = NULL  # or matrix(nrow=0, ncol=0) if you prefer
        )
      }
    }
  })
  names(model_results) <- NoCov_mPACE_models$model_names
  return(model_results)
}




