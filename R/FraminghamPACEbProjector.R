#' @noRd
globalVariables("FraminghamPACEb_RIDGE")
globalVariables("FraminghamPACEb_ENET")
globalVariables(c(".validate_betas"))


FraminghamPACEbProjector = function(betas, 
                                    proportionOfProbesRequired=0.8, 
                                    model_type="RIDGE",
                                    return_normalized_betas = FALSE) {
  requireNamespace("preprocessCore")
  
  # Load the appropriate model based on model_type
  if (model_type == "RIDGE") {
    data("FraminghamPACEb_RIDGE", package = "FraminghamPACE")
    model_data <- FraminghamPACEb_RIDGE
  } else if (model_type == "ENET") {
    data("FraminghamPACEb_ENET", package = "FraminghamPACE")
    model_data <- FraminghamPACEb_ENET
  } else {
    stop("Invalid model_type. Options are 'RIDGE' or 'ENET'.")
  }
  
  betas <- .validate_betas(betas)
  
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

  # Get the model name
  model_name <- model_data$model_names[1]

  # Assess predictive coverage before any imputation or threshold gating
  try({
    .assess_predictive_coverage(betas, model_data, model_name, verbose = TRUE)
  }, silent = TRUE)
  
  if( !is.numeric(as.matrix(betas)) ) { stop("betas matrix/data.frame is not numeric!") }
  probeOverlap <- length(which(rownames(betas) %in% model_data$model_probes[[model_name]])) / 
    length(model_data$model_probes[[model_name]])
  probeOverlap_background <- length(which(rownames(betas) %in% model_data$gold_standard_probes[[model_name]])) / 
    length(model_data$gold_standard_probes[[model_name]])
  
  if( probeOverlap < proportionOfProbesRequired | probeOverlap_background < proportionOfProbesRequired ) {
    result <- rep(NA, ncol(betas))
    names(result) <- colnames(betas)
    
    pace_result <- data.frame(result)
    colnames(pace_result) <- model_name
    rownames(pace_result) <- colnames(betas)
    
    cat("WARNING: Insufficient number of probes present in the data file\n")
    
    if(return_normalized_betas) {
      return(list(
        pace_score = pace_result,
        normalized_betas = NULL
      ))
    } else {
      return(pace_result)
    }
    
  }
  
  # Work with a numeric matrix of betas
  betas.mat <- as.matrix(betas[which(rownames(betas) %in% model_data$gold_standard_probes[[model_name]]),])
  
  # If probes don't exist, we'll add them as rows of values based on their mean in the gold standard dataset
  probesNotInMatrix <- model_data$gold_standard_probes[[model_name]][which(model_data$gold_standard_probes[[model_name]] %in% rownames(betas.mat) == F)]
  if( length(probesNotInMatrix) > 0 ) {
    for( probe in probesNotInMatrix ) {
      tmp.mat <- matrix(0, nrow=1, ncol=ncol(betas.mat))
      rownames(tmp.mat) <- probe
      colnames(tmp.mat) <- colnames(betas.mat)
      tmp.mat[probe,] <- rep(model_data$gold_standard_means[[model_name]][probe], ncol(tmp.mat))
      betas.mat <- rbind(betas.mat, tmp.mat)
    }
  }
  
  # Identify samples with too many missing probes and remove them from the matrix
  samplesToRemove <- colnames(betas.mat)[which(apply(betas.mat, 2, function(x) { 1 - ( length(which(is.na(x))) / length(x) ) < proportionOfProbesRequired}))]
  
  if( length(samplesToRemove) > 0 ) {
    betas.mat <- betas.mat[,-which(colnames(betas.mat) %in% samplesToRemove)]
  }
  
  if(ncol(betas.mat) == 0) {
    result <- rep(NA, ncol(betas))
    names(result) <- colnames(betas)
    
    pace_result <- data.frame(result)
    colnames(pace_result) <- model_name
    rownames(pace_result) <- colnames(betas)
    
    cat("WARNING: Empty beta matrix after removing samples with excessive missingness\n")
    
    if(return_normalized_betas) {
      return(list(
        pace_score = pace_result,
        normalized_betas = NULL
      ))
    } else {
      return(pace_result)
    }
  }
  
  # Identify missingness on a probe level
  pctValuesPresent <- apply( betas.mat, 1, function(x) { 1 - (length(which(is.na(x))) / length(x)) } )
  
  # If they're missing values, but less than the proportion required, we impute to the cohort mean
  probesToAdjust <- which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired)
  
  # Case 1: Multiple probes to adjust
  if( length(probesToAdjust) > 1 ) {
    betas.mat[probesToAdjust,] <- t(apply( betas.mat[probesToAdjust,], 1 , function(x) {
      x[is.na(x)] = mean( x, na.rm = TRUE )
      x
    }))
  }
  
  # Case 2: Exactly one probe to adjust
  if(length(probesToAdjust) == 1) {
    betas.mat[probesToAdjust,which(is.na(betas.mat[probesToAdjust,]))] <- mean(betas.mat[probesToAdjust,], na.rm=T)
  }
  
  # If they're missing too many values, everyones value gets replaced with the mean from the Framingham cohort
  if( length(which(pctValuesPresent < proportionOfProbesRequired)) > 0 ) {
    probesToReplaceWithMean <- rownames(betas.mat)[which(pctValuesPresent < proportionOfProbesRequired)]
    for( probe in probesToReplaceWithMean ) {
      betas.mat[probe,] <- rep(model_data$model_means[[model_name]][probe], ncol(betas.mat))
    }
  }
  
  # Normalize the matrix to the gold standard dataset
  betas.norm <- preprocessCore::normalize.quantiles.use.target(betas.mat, target=model_data$gold_standard_means[[model_name]])
  rownames(betas.norm) <- rownames(betas.mat)
  colnames(betas.norm) <- colnames(betas.mat)
  
  # Calculate score:
  score = model_data$model_intercept[[model_name]] + 
    rowSums(t(betas.norm[model_data$model_probes[[model_name]],]) %*% 
              diag(model_data$model_weights[[model_name]]))
  
  names(score) <- colnames(betas.norm)
  
  if( length(samplesToRemove) > 0 ) {
    score.tmp <- rep(NA, length(samplesToRemove))
    names(score.tmp) <- samplesToRemove
    score <- c(score, score.tmp)
  }
  
  score <- score[colnames(betas)]
  
  pace_result <- data.frame(score)
  colnames(pace_result) <- model_name
  
  if(return_normalized_betas) {
    return(list(
      pace_score = pace_result,
      normalized_betas = betas.norm
    ))
  } else {
    return(pace_result)
  }
  
}