#' @noRd
globalVariables("FraminghamPACEb_RIDGE")
globalVariables("FraminghamPACEb_ENET")
globalVariables(c(".validate_betas"))


FraminghamPACEbProjector = function(betas, 
                                    proportionOfProbesRequired=0.8, 
                                    model_type="RIDGE",
                                    return_normalized_betas = FALSE,
                                    sample_coverage_min = 0.75,
                                    return_coverage = FALSE) {
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
    print("Averaging replicates...")
    
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
  }

  # Get the model name
  model_name <- model_data$model_names[1]


  if( !is.numeric(as.matrix(betas)) ) { stop("betas matrix/data.frame is not numeric!") }
  
  # Compute per-sample coverage prior to imputation; mark low-coverage samples
  sample_cov <- .compute_sample_coverage(betas, model_data, model_name)
  low_cov_samples <- names(sample_cov)[which(sample_cov < (sample_coverage_min * 100))]
  if (length(low_cov_samples) > 0) {
    cli::cli_alert_warning("{length(low_cov_samples)} samples below coverage threshold ({sample_coverage_min*100}%). Scores will be set to NA for these samples.")
  }

  # Report how many samples remain after filtering; short-circuit if none remain
  remaining_count <- length(setdiff(colnames(betas), low_cov_samples))
  if (remaining_count == 0) {
    cli::cli_alert_danger("No samples remain after filtering at {sample_coverage_min*100}% coverage. Consider lowering 'sample_coverage_min'.")
    return(invisible(NULL))
  }
  cli::cli_alert_success("{remaining_count} samples remain after filtering at {sample_coverage_min*100}% coverage (from {ncol(betas)}).")

  # Report dataset-level predictive coverage on the kept samples (pre-imputation)
  try({
    kept <- setdiff(colnames(betas), low_cov_samples)
    betas_for_coverage <- if (length(kept) > 0) betas[, kept, drop = FALSE] else betas[, 0, drop = FALSE]
    .assess_predictive_coverage(betas_for_coverage, model_data, model_name, verbose = TRUE)
  }, silent = TRUE)
  
  # Work with a numeric matrix of betas
  betas.mat <- as.matrix(betas[which(rownames(betas) %in% model_data$gold_standard_probes[[model_name]]),])
  
  # If probes don't exist, we'll add them as rows of values based on their mean in the gold standard dataset
  probesNotInMatrix <- model_data$gold_standard_probes[[model_name]][which(model_data$gold_standard_probes[[model_name]] %in% rownames(betas.mat) == F)]
  if( length(probesNotInMatrix) > 0 ) {
    means_vec <- model_data$gold_standard_means[[model_name]][probesNotInMatrix]
    add_mat <- tcrossprod(means_vec, rep(1, ncol(betas.mat)))
    rownames(add_mat) <- probesNotInMatrix
    colnames(add_mat) <- colnames(betas.mat)
    betas.mat <- rbind(betas.mat, add_mat)
  }
  
  # Identify samples to remove based on coverage threshold
  samplesToRemove <- intersect(colnames(betas.mat), low_cov_samples)
  
  if( length(samplesToRemove) > 0 ) {
    betas.mat <- betas.mat[,-which(colnames(betas.mat) %in% samplesToRemove)]
  }
  
  if(ncol(betas.mat) == 0) {
    result <- rep(NA, ncol(betas))
    names(result) <- colnames(betas)
    
    pace_result <- data.frame(result)
    colnames(pace_result) <- model_name
    rownames(pace_result) <- colnames(betas)
    
    cat("WARNING: Empty beta matrix after removing low-coverage samples\n")
    
    if(return_normalized_betas || return_coverage) {
      return(list(
        pace_score = pace_result,
        normalized_betas = NULL,
        coverage = list(
          sample = sample_cov,
          threshold = sample_coverage_min
        )
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
  weighted_sum <- colSums(betas.norm[model_data$model_probes[[model_name]], , drop = FALSE] * 
                            model_data$model_weights[[model_name]])
  score = model_data$model_intercept[[model_name]] + weighted_sum
  
  names(score) <- colnames(betas.norm)
  
  if( length(samplesToRemove) > 0 ) {
    score.tmp <- rep(NA, length(samplesToRemove))
    names(score.tmp) <- samplesToRemove
    score <- c(score, score.tmp)
  }
  
  score <- score[colnames(betas)]
  
  pace_result <- data.frame(score)
  colnames(pace_result) <- model_name
  
  if(return_normalized_betas || return_coverage) {
    return(list(
      pace_score = pace_result,
      normalized_betas = if (isTRUE(return_normalized_betas)) betas.norm else NULL,
      coverage = list(
        sample = sample_cov,
        threshold = sample_coverage_min
      )
    ))
  }
  return(pace_result)
  
}