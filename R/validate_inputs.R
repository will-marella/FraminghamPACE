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
    
    # Check if dimensions suggest possible transposition
    # We expect many probes (rows) and fewer samples (columns)
    if (nrow(betas) < ncol(betas)) {
      warning(
        "The 'betas' matrix has fewer rows (", nrow(betas), ") than columns (", ncol(betas), "). ",
        "Rows should represent probes and columns should represent samples. ",
        "Please verify the matrix orientation is correct."
      )
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
  tryCatch({
    stopifnot(
      "'age_sex_df' must be provided" = !is.null(age_sex_df),
      "'age_sex_df' must be a data frame" = is.data.frame(age_sex_df),
      "'age_sex_df' must contain 'Age' and 'Female' columns" = all(c("Age", "Female") %in% colnames(age_sex_df)),
      "'Age' column must be numeric" = is.numeric(age_sex_df$Age),
      "'Female' column must be numeric" = is.numeric(age_sex_df$Female),
      "'Female' column must contain only 0 (male) or 1 (female)" = all(age_sex_df$Female %in% c(0, 1)),
      "'Age' column cannot contain negative values" = all(age_sex_df$Age >= 0),
      "All samples in betas must have corresponding entries in age_sex_df" = all(betas_samples %in% rownames(age_sex_df))
    )
    
    # Check age centering and range
    mean_age <- mean(age_sex_df$Age)
    if (abs(mean_age - 65) > 30) {  # Allow some deviation but catch major issues
      warning(sprintf(
        "Mean age (%.1f) is far from expected mean of 65. Please verify age values",
        mean_age
      ))
    }
    
    # Check for suspiciously high or low ages relative to centering
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