#' Generate the Framingham Methylation Pace of Aging Scores!
#'
#' \code{FraminghamPACEaProjector} calculates and returns FraminghamPACEa scores from an individual’s methylation array, age, and sex.
#'
#' @name FraminghamPACEaProjector
#' @rdname FraminghamPACEaProjector
#' @param betas A numeric matrix containing the percent-methylation for each probe.  Missing data should be 'NA's.  The rows should be probes, with the probe ID as the row name, and the columns should be samples, with sample names as the column name.
#' @param age_sex_df a numeric dataframe containing the age and sex values for each sample. The name of the age column should be 'Age' and the name of the sex column should be 'Female', where a 1 indicates the sample is female and 0 indicates the sample is male. The row names should indicate the sample, matching the column names of the `betas` matrix.
#' @param proportionOfProbesRequired (default: 0.8).  This value specifies the threshold for missing data (see description for more details on how missing data is handled). Note: if the function detects EPICv2 data based on row names that contain the replicate suffix, the proportion of probes required will automatically be lowered to 0.7.
#' @param model_type (default: "default"). Choose which model to use. Options are "default" for the original model or "MSA" for the alternative MSA model.
#' @param return_normalized_betas (default: FALSE). If FALSE, returns a dataframe of FraminghamPACEa values. If TRUE, returns a list containing both the dataframe of FraminghamPACEa values and the normalized beta matrix.
#' @return By default, a dataframe of FraminghamPACEa values with sample names as row names and the model name as column name. If return_normalized_betas=TRUE, a list containing both the pace_score dataframe and the normalized_betas matrix.
#' @details This function returns the Framingham Methylation Pace of Aging scores for methylation data generated from either the Illumina 450K array or the Illumina EPIC array.  First, if a sample is missing data for more probes than the threshold, the sample will get an `NA` back for a score.  If a particular probe is missing fewer samples than the threshold, then missing data is set to the mean in the provided 'betas' matrix.  If a probe is missing more samples than the threshold, then all samples in the 'betas' matrix have their value replaced with the mean of the training data for that particular model.
#' Because of how we handle missing data, it is recommended that entire cohorts be run at once as a large 'betas' matrix.
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @importFrom utils data
#' @export
globalVariables("FraminghamPACEa_model")
globalVariables("FraminghamPACEa_model_MSA")
globalVariables("FraminghamPACEas_model_TruDx")
globalVariables(c(".validate_betas", ".validate_age_sex", ".validate_and_align_samples"))


FraminghamPACEaProjector = function(betas, age_sex_df, 
                                    proportionOfProbesRequired=0.8, model_type="default",
                                    return_normalized_betas = FALSE) {
  requireNamespace("preprocessCore")
  
  # Load the appropriate model based on model_type
  if (model_type == "default") {
    data("FraminghamPACEa_model", package = "FraminghamPACE")
    model_data <- FraminghamPACEa_model
  } else if (model_type == "MSA") {
    data("FraminghamPACEa_model_MSA", package = "FraminghamPACE")
    model_data <- FraminghamPACEa_model_MSA
  } else if (model_type == "TruDx") {
    data("FraminghamPACEas_model_TruDx", package = "FraminghamPACE")
    model_data <- FraminghamPACEas_model_TruDx
  } else {
    stop("Invalid model_type. Options are 'default' or 'MSA'.")
  }
  
  # Individual validation of inputs
  betas <- .validate_betas(betas)
  age_sex_df <- .validate_age_sex(age_sex_df, colnames(betas))
  
  # Align samples between betas and age_sex_df
  aligned_data <- .validate_and_align_samples(betas, age_sex_df)
  betas <- aligned_data$betas
  age_sex_df <- aligned_data$age_sex_df
  
  ################################################################
  ### TEMPORARY ADDITION TO ATTEMPT AN AGE/SEX-INDEPENDENT SCORE
  age_sex_df$Age <- 0
  age_sex_df$Female <- 0.5
  ################################################################
  
  
  if( any(grepl("TC|BC", rownames(betas))) )
  {
    # Print message to user
    print("This looks like EPICv2 array data. If EPICv2, FraminghamPACE will lower the proportionOfProbesRequired to 0.7 and proceed with missing probes present on 450k and EPICv1. Averaging the opposite strand replicates may take a bit more time.")
    
    # Average replicates: Code taken from ENmix package.
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
  
  if( !is.numeric(as.matrix(betas)) ) { stop("betas matrix/data.frame is not numeric!") }
  probeOverlap <- length(which(rownames(betas) %in% model_data$model_probes[[model_name]])) / 
    length(model_data$model_probes[[model_name]])
  probeOverlap_background <- length(which(rownames(betas) %in% model_data$gold_standard_probes[[model_name]])) / 
    length(model_data$gold_standard_probes[[model_name]])
  
  # make sure enough of the probes are present in the data file
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
    rowSums(t(betas.norm[model_data$model_probes[[model_name]],]) %*% diag(model_data$model_weights[[model_name]])) +
    age_sex_df$Age * model_data$age_weight[[model_name]] +
    age_sex_df$Female * model_data$female_weight[[model_name]]
  
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
