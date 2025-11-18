#' Internal helpers (not exported)
#' @noRd
.with_timing <- function(label, expr, timing) {
  if (!isTRUE(timing)) {
    return(eval.parent(substitute(expr)))
  }
  t <- system.time(res <- eval.parent(substitute(expr)))
  cli::cli_alert_info("{label}: {round(t['elapsed'], 3)}s")
  res
}

#' Extract base CpG IDs (strip suffix after first underscore if present)
#' @noRd
.extract_base_ids <- function(rownames_vec) {
  sub("_.*$", "", rownames_vec)
}

#' Compute per-sample presence matrix by base CpG ID
#' Returns a logical matrix: bases x samples
#' @noRd
.present_by_base <- function(betas) {
  base_ids <- .extract_base_ids(rownames(betas))
  present_counts <- rowsum((!is.na(betas)) * 1L, group = base_ids, reorder = FALSE)
  present_counts > 0L
}

#' Replicate-aware prefilter to background by base ID (no renaming)
#' @noRd
.prefilter_to_background <- function(betas, background_cg_ids, timing = FALSE) {
  .with_timing("Prefilter to background by base ID", {
    rn <- rownames(betas)
    base_ids <- .extract_base_ids(rn)
    in_bg <- base_ids %in% background_cg_ids
    if (!all(in_bg)) {
      betas <- betas[in_bg, , drop = FALSE]
    }
    betas
  }, timing = timing)
}

#' Collapse replicate rows to base IDs using NA-aware means via rowsum
#' Mirrors previous semantics: if duplicates exist, collapse; otherwise rename.
#' @noRd
.collapse_replicates_rowsum <- function(betas, timing = FALSE) {
  .with_timing("Replicate averaging", {
    rn <- rownames(betas)
    base_ids <- .extract_base_ids(rn)
    tab <- table(base_ids)
    dup_ids <- names(tab)[tab > 1L]
    if (length(dup_ids) > 0L) {
      sel_dup <- base_ids %in% dup_ids
      betas_dup <- betas[sel_dup, , drop = FALSE]
      groups <- base_ids[sel_dup]
      sums <- rowsum(betas_dup, groups, na.rm = TRUE)
      counts <- rowsum((!is.na(betas_dup)) * 1L, groups)
      counts_num <- pmax(counts, 1L)
      means <- sums / counts_num
      means[counts == 0] <- NA_real_
      betas_nondup <- betas[!sel_dup, , drop = FALSE]
      rownames(betas_nondup) <- base_ids[!sel_dup]
      betas <- rbind(betas_nondup, as.matrix(means))
    } else {
      rownames(betas) <- base_ids
    }
    betas
  }, timing = timing)
}


