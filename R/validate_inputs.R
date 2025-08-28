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
