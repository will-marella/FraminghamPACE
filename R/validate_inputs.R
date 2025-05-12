#' Validate and align samples between betas matrix and age_sex dataframe
#'
#' @param betas A numeric matrix of methylation values
#' @param age_sex_df A data frame containing age and sex information
#' @return A list containing aligned betas matrix and age_sex dataframe
#' @noRd
.validate_and_align_samples <- function(betas, age_sex_df) {
  # Get sample names from both inputs
  beta_samples <- colnames(betas)
  age_sex_samples <- rownames(age_sex_df)
  
  # Check for sample presence
  common_samples <- intersect(beta_samples, age_sex_samples)
  
  if (length(common_samples) == 0) {
    stop("No matching sample names found between betas and age_sex_df")
  }
  
  # Calculate removed samples
  removed_beta_samples <- setdiff(beta_samples, common_samples)
  removed_age_sex_samples <- setdiff(age_sex_samples, common_samples)
  
  # Print information about sample matching
  if (length(removed_beta_samples) > 0 || length(removed_age_sex_samples) > 0) {
    message(sprintf(
      "Sample alignment: %d samples matched. Removed %d samples from betas and %d from age_sex_df",
      length(common_samples),
      length(removed_beta_samples),
      length(removed_age_sex_samples)
    ))
  } else {
    message("All samples matched successfully")
  }
  
  # Align the data
  betas_aligned <- betas[, common_samples, drop = FALSE]
  age_sex_aligned <- age_sex_df[common_samples, , drop = FALSE]
  
  # Verify alignment
  stopifnot(
    "Sample alignment failed" = identical(colnames(betas_aligned), rownames(age_sex_aligned)),
    "No samples remained after alignment" = ncol(betas_aligned) > 0
  )
  
  return(list(
    betas = betas_aligned,
    age_sex_df = age_sex_aligned,
    removed_samples = list(
      betas = removed_beta_samples,
      age_sex = removed_age_sex_samples
    )
  ))
}

#' Validate beta values matrix/data frame
#'
#' @param betas A matrix or data frame containing methylation values
#' @return A validated numeric matrix of methylation values
#' @noRd
.validate_betas <- function(betas) {
  
  tryCatch({
    stopifnot(
      "'betas' must be provided" = !is.null(betas),
      "'betas' must be a matrix or data frame" = is.matrix(betas) || is.data.frame(betas),
      "'betas' cannot be empty" = nrow(betas) > 0 && ncol(betas) > 0
    )
    
    # Convert to matrix if needed
    if (is.data.frame(betas)) {
      message("Converting 'betas' from data frame to matrix")
      betas <- as.matrix(betas)
    }
    
    stopifnot(
      "'betas' must contain numeric values" = is.numeric(betas)
    )
    
    # Check probe names format
    if (!all(grepl("^cg", rownames(betas)))) {
      warning("Not all row names start with 'cg'. Please verify probe identifiers are correct.")
    }
    
    # Check matrix orientation
    if (nrow(betas) < ncol(betas)) {
      warning(
        "The 'betas' matrix has fewer rows (", nrow(betas), ") than columns (", ncol(betas), "). ",
        "Typically we expect more CpG probes than samples. Please verify matrix orientation."
      )
    }
    
    # Check NA patterns
    na_counts <- colSums(is.na(betas))
    if (any(na_counts > 0)) {
      message(sprintf(
        "Found samples with missing values. Range: %d-%d NAs per sample", 
        min(na_counts[na_counts > 0]), 
        max(na_counts)
      ))
    }
    
    return(betas)
  }, error = function(e) {
    stop(e$message, call. = FALSE)
  })
}

#' Validate age and sex data frame
#'
#' @param age_sex_df A data frame containing age and sex information
#' @param betas_samples Character vector of sample names from betas matrix
#' @noRd
.validate_age_sex <- function(age_sex_df, betas_samples) {
  
  # Basic validation
  stopifnot(
    "'age_sex_df' must be provided" = !is.null(age_sex_df),
    "'age_sex_df' must be a data frame" = is.data.frame(age_sex_df),
    "'age_sex_df' must contain 'Age' and 'Female' columns" = all(c("Age", "Female") %in% colnames(age_sex_df))
  )
  
  tryCatch({
    # Check for NA values first - this should error and explain options
    if(any(is.na(age_sex_df$Age) | is.na(age_sex_df$Female))) {
      na_samples <- rownames(age_sex_df)[is.na(age_sex_df$Age) | is.na(age_sex_df$Female)]
      stop(sprintf(
        "Found %d samples with missing age or sex values:\n%s\n\nOptions:\n%s",
        length(na_samples),
        paste(na_samples, collapse=", "),
        paste(
          "1. Remove these samples from age_sex_df to use FraminghamPACEa with complete cases only",
          "2. Use FraminghamPACEb instead for these samples (does not require age/sex)",
          "3. Impute values before passing to this function",
          sep="\n"
        )
      ))
    }
    
    # Basic numeric validation
    stopifnot(
      "'Age' column must be numeric" = is.numeric(age_sex_df$Age),
      "'Female' column must be numeric" = is.numeric(age_sex_df$Female),
      "'Age' column cannot contain negative values" = all(age_sex_df$Age >= 0),
      "'Female' column must contain values between 0 and 1" = all(age_sex_df$Female >= 0 & age_sex_df$Female <= 1)    )
    
    # 4. Warning for non-binary sex values
    if(any(age_sex_df$Female > 0 & age_sex_df$Female < 1)) {
      decimal_samples <- rownames(age_sex_df)[age_sex_df$Female > 0 & age_sex_df$Female < 1]
      warning(sprintf(
        paste(
          "Found %d samples with non-binary sex values (between 0 and 1):",
          "%s",
          "These values will be used as weights for the female coefficient in the model.",
          "If this is not intended, please update to use only 0 (male) or 1 (female).",
          sep="\n"
        ),
        length(decimal_samples),
        paste(decimal_samples, collapse=", ")
      ))
    }
    
    # Check age centering and range
    mean_age <- mean(age_sex_df$Age)
    if (abs(mean_age - 65) > 30) {  # Allow some deviation but catch major issues
      warning(sprintf(
        "Mean age (%.1f) is far from expected mean of 65. Please verify age values",
        mean_age
      ))
    }
    
    # Check for suspiciously high or low ages
    if (any(age_sex_df$Age > 120)) {
      warning("Some ages are over 120 years - please verify these values are correct")
    }
    if (any(age_sex_df$Age < 10)) {
      warning("Some ages are under 10 years - please verify these values are correct")
    }
    
    # Center age at 65
    age_sex_df$Age <- age_sex_df$Age - 65
    message("Ages have been centered by subtracting 65 years")
    
    return(age_sex_df)
  }, error = function(e) {
    stop(e$message, call. = FALSE)
  })
}