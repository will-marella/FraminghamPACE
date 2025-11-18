#' Validate beta values matrix/data frame
#'
#' @param betas A matrix or data frame containing methylation values
#' @return A validated numeric matrix of methylation values
#' @noRd
.validate_betas <- function(betas, timing = FALSE) {
  
  tryCatch({
    time_stage <- function(label, expr) {
      if (!isTRUE(timing)) {
        return(eval.parent(substitute(expr)))
      }
      t <- system.time(res <- eval.parent(substitute(expr)))
      cli::cli_alert_info("{label}: {round(t['elapsed'], 3)}s")
      res
    }
    
    stopifnot(
      "'betas' must be provided" = !is.null(betas),
      "'betas' must be a matrix or data frame" = is.matrix(betas) || is.data.frame(betas),
      "'betas' cannot be empty" = nrow(betas) > 0 && ncol(betas) > 0
    )
    
    # Convert to matrix if needed
    if (is.data.frame(betas)) {
      message("Converting 'betas' from data frame to matrix")
      betas <- time_stage("Convert data.frame to matrix", as.matrix(betas))
    }
    
    stopifnot(
      "'betas' must contain numeric values" = is.numeric(betas)
    )
    
    # Check probe names format; filter to cg-only rows for robustness
    cg_mask <- time_stage("Compute cg probe mask", grepl("^cg", rownames(betas)))
    if (!all(cg_mask)) {
      removed <- sum(!cg_mask)
      message(sprintf("Filtering out %d non-'cg' rows prior to scoring", removed))
      betas <- time_stage("Filter to cg probes", betas[cg_mask, , drop = FALSE])
      if (nrow(betas) == 0) {
        stop("After filtering, no 'cg' probes remain. Please provide Illumina CpG IDs as row names.")
      }
    }
    
    # Check matrix orientation
    if (nrow(betas) < ncol(betas)) {
      warning(
        "The 'betas' matrix has fewer rows (", nrow(betas), ") than columns (", ncol(betas), "). ",
        "Typically we expect more CpG probes than samples. Please verify matrix orientation."
      )
    }
    
    # Check NA patterns
    na_counts <- time_stage("NA scan across matrix", colSums(is.na(betas)))
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
  rn <- rownames(betas)
  n_samples <- ncol(betas)

  # Vectorized precomputations
  idx_bg <- match(background_probes, rn, nomatch = 0L)
  idx_sc <- match(scoring_probes, rn, nomatch = 0L)
  nonmiss_counts <- rowSums(!is.na(betas))

  # Availability counts (present in betas and at least one non-missing value)
  background_total <- length(background_probes)
  scoring_total <- length(scoring_probes)
  background_available <- sum(idx_bg > 0L & nonmiss_counts[idx_bg] > 0L)
  scoring_available <- sum(idx_sc > 0L & nonmiss_counts[idx_sc] > 0L)

  # Predictive coverage: importance weighted by non-missing fraction per probe
  importance <- model_data$model_importance_pct[[model_name]]
  imp_names <- names(importance)
  idx_imp <- match(imp_names, rn, nomatch = 0L)
  keep <- idx_imp > 0L
  if (any(keep)) {
    frac_vec <- nonmiss_counts[idx_imp[keep]] / n_samples
    coverage_pct <- round(100 * sum(importance[keep] * frac_vec, na.rm = TRUE), 1)
  } else {
    coverage_pct <- NA_real_
  }

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
  
  # Collapse presence by base ID so suffixes (e.g., _TC/_BC/_orig) don't break matching
  base_ids <- sub("_.*$", "", rownames(betas))
  present_counts_by_base <- rowsum((!is.na(betas)) * 1L, group = base_ids, reorder = FALSE)
  present_by_base <- present_counts_by_base > 0L
  
  # Build present matrix aligned to scoring_probes
  n_samples <- ncol(betas)
  present <- matrix(FALSE, nrow = length(scoring_probes), ncol = n_samples,
                    dimnames = list(scoring_probes, colnames(betas)))
  idx <- match(scoring_probes, rownames(present_by_base), nomatch = 0L)
  if (any(idx > 0L)) {
    present[idx > 0L, ] <- present_by_base[idx[idx > 0L], , drop = FALSE]
  }
  
  # Importance vector aligned to scoring_probes; importance sums to 1
  imp_vec <- importance[scoring_probes]
  imp_vec[is.na(imp_vec)] <- 0
  
  # colSums of importance weights for present probes
  coverage <- as.numeric(crossprod(imp_vec, present)) * 100
  names(coverage) <- colnames(betas)
  return(coverage)
}
