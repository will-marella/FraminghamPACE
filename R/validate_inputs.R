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

#' Assess dataset-level predictive coverage and probe availability
#'
#' Computes background/scoring probe availability and a weighted predictive
#' coverage percentage using model probe importance and probe missingness.
#' Prints a brief summary when verbose.
#'
#' @param betas Numeric matrix/data frame of beta values (probes x samples)
#' @param model_data Loaded model data list containing probe sets and importance
#' @param model_name Character model name within model_data
#' @param verbose Logical; print messages (default TRUE)
#' @return A list with background/scoring availability and coverage metrics
#' @noRd
.assess_predictive_coverage <- function(betas, model_data, model_name, verbose = TRUE) {
  # Ensure matrix
  if (is.data.frame(betas)) betas <- as.matrix(betas)

  # Background and scoring probe sets
  background_probes <- model_data$gold_standard_probes[[model_name]]
  scoring_probes <- model_data$model_probes[[model_name]]

  # Availability: probe present and has at least one non-NA value
  probe_has_data <- function(probe) {
    if (!probe %in% rownames(betas)) return(FALSE)
    any(!is.na(betas[probe, , drop = FALSE]))
  }

  background_total <- length(background_probes)
  background_available <- sum(vapply(background_probes, probe_has_data, logical(1)))

  scoring_total <- length(scoring_probes)
  scoring_available <- sum(vapply(scoring_probes, probe_has_data, logical(1)))

  # Predictive coverage: importance weighted by non-missing fraction per probe
  importance <- model_data$model_importance_pct[[model_name]]
  probes_in_imp <- intersect(names(importance), rownames(betas))

  # Fraction present per probe (0 if fully missing)
  non_missing_fraction <- function(probe) {
    vals <- betas[probe, , drop = FALSE]
    1 - (sum(is.na(vals)) / length(vals))
  }

  frac_vec <- vapply(probes_in_imp, non_missing_fraction, numeric(1))
  # Weight by importance and sum to percentage (importance sums to 1)
  coverage_pct <- round(100 * sum(importance[probes_in_imp] * frac_vec, na.rm = TRUE), 1)

  # Coverage label
  if (is.na(coverage_pct)) {
    coverage_label <- NA_character_
  } else if (coverage_pct > 90) {
    coverage_label <- "Good"
  } else if (coverage_pct > 75) {
    coverage_label <- "Moderate"
  } else {
    coverage_label <- "Poor"
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_info("Probes: {background_available}/{background_total} background | {scoring_available}/{scoring_total} scoring")
    if (!is.na(coverage_pct)) {
      cli::cli_alert_info("Predictive coverage: {coverage_pct}% [{coverage_label}]")
    }
    cli::cli_text("More: https://github.com/will-marella/FraminghamPACE")
  }

  return(list(
    background = list(available = background_available, total = background_total),
    scoring = list(available = scoring_available, total = scoring_total),
    coverage_pct = coverage_pct,
    coverage_label = coverage_label
  ))
}

#' Compute per-sample predictive coverage percentages
#'
#' Coverage per sample is the sum of model probe importances for probes that are
#' non-missing in that sample, multiplied by 100. Importance sums to 1 across
#' scoring probes, so coverage ranges from 0 to 100.
#'
#' @param betas Numeric matrix/data frame of beta values (probes x samples)
#' @param model_data Loaded model data list containing probe sets and importance
#' @param model_name Character model name within model_data
#' @return Named numeric vector of per-sample coverage percentages
#' @noRd
.compute_sample_coverage <- function(betas, model_data, model_name) {
  if (is.data.frame(betas)) betas <- as.matrix(betas)
  scoring_probes <- model_data$model_probes[[model_name]]
  importance <- model_data$model_importance_pct[[model_name]]
  probes <- intersect(scoring_probes, rownames(betas))
  if (length(probes) == 0) {
    return(setNames(rep(0, ncol(betas)), colnames(betas)))
  }
  imp_vec <- importance[probes]
  present <- !is.na(betas[probes, , drop = FALSE])
  # colSums of importance weights for present probes; importance sums to 1
  coverage <- as.numeric(crossprod(imp_vec, present)) * 100
  names(coverage) <- colnames(betas)
  return(coverage)
}
