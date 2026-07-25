###############################################################################
# Exported Functions
###############################################################################

.validate_intensity_matrix <- function(data,
                                       context = "data",
                                       require_nonnegative = FALSE,
                                       allow_nonfinite = FALSE) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop(context, " must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(data) < 1L || ncol(data) < 1L) {
    stop(context, " must contain at least one row and one column.", call. = FALSE)
  }
  if (!allow_nonfinite && any(!is.finite(data))) {
    stop(
      context,
      " must contain only finite values. Convert declared non-detections to ",
      "missing values and impute them before running WiNN.",
      call. = FALSE
    )
  }
  if (require_nonnegative && any(data < 0, na.rm = TRUE)) {
    stop(
      context,
      " must contain non-negative intensities because WiNN uses log1p. ",
      "Values at or below -1 are undefined and other negative values are not ",
      "valid abundance inputs.",
      call. = FALSE
    )
  }
  invisible(data)
}

.validate_batch <- function(batch, n_samples, required = TRUE) {
  if (is.null(batch)) {
    if (isTRUE(required)) {
      stop("batch must be supplied.", call. = FALSE)
    }
    return(invisible(NULL))
  }
  if (is.matrix(batch) || is.data.frame(batch) || length(batch) != n_samples) {
    stop("batch must be a vector with one value per data column.", call. = FALSE)
  }
  if (anyNA(batch)) {
    stop("batch must not contain missing values.", call. = FALSE)
  }
  if (is.numeric(batch) && any(!is.finite(batch))) {
    stop("Numeric batch labels must be finite.", call. = FALSE)
  }
  batch_text <- trimws(as.character(batch))
  if (any(!nzchar(batch_text))) {
    stop("batch labels must not be empty.", call. = FALSE)
  }
  invisible(batch)
}

.validate_run_order <- function(run_order, batch, n_samples) {
  if (is.null(run_order)) {
    return(invisible(NULL))
  }
  if (!is.numeric(run_order) || is.matrix(run_order) ||
      length(run_order) != n_samples || any(!is.finite(run_order))) {
    stop(
      "run_order must be a finite numeric vector with one value per data column.",
      call. = FALSE
    )
  }
  groups <- if (is.null(batch)) {
    list(all_samples = seq_len(n_samples))
  } else {
    split(seq_len(n_samples), as.character(batch), drop = TRUE)
  }
  duplicate_groups <- names(groups)[vapply(groups, function(idx) {
    anyDuplicated(run_order[idx]) > 0L
  }, logical(1))]
  if (length(duplicate_groups)) {
    stop(
      "run_order must be unique within each batch; duplicates found in: ",
      paste(duplicate_groups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(run_order)
}

.validate_control_samples <- function(control_samples, n_samples) {
  if (is.null(control_samples)) {
    return(invisible(NULL))
  }
  if (!is.numeric(control_samples) || any(!is.finite(control_samples)) ||
      any(control_samples != round(control_samples)) ||
      any(control_samples < 1L | control_samples > n_samples) ||
      anyDuplicated(control_samples)) {
    stop(
      "control_samples must contain unique, finite column indices in data.",
      call. = FALSE
    )
  }
  invisible(as.integer(control_samples))
}

.validate_probability <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0 || value >= 1) {
    stop(name, " must be one finite number strictly between 0 and 1.", call. = FALSE)
  }
  invisible(value)
}

.validate_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(value)
}

.validate_ljung_box_fitdf <- function(ljung_box_fitdf) {
  if (!is.numeric(ljung_box_fitdf) || length(ljung_box_fitdf) != 1L ||
      !is.finite(ljung_box_fitdf) || ljung_box_fitdf < 0 ||
      ljung_box_fitdf != round(ljung_box_fitdf)) {
    stop("ljung_box_fitdf must be a single non-negative integer.", call. = FALSE)
  }
  as.integer(ljung_box_fitdf)
}

.ordered_batch_indices <- function(batch, batch_value, run_order = NULL) {
  idx <- which(batch == batch_value)
  if (!is.null(run_order)) {
    idx <- idx[order(run_order[idx])]
  }
  idx
}

#' Normalize Data Matrix by Dilution Factor
#'
#' This function normalizes a data matrix by dilution factors or alternatively shrinks sample measurements
#' if any samples have dilution factors that are more than one standard deviation away from the mean dilution factor.
#' Numeric zeros remain observed values. Features whose reference median is zero
#' are omitted from quotient fitting, but retained when the estimated sample
#' factors are applied. A sample with a non-positive dilution factor is rejected.
#'
#' @param data A numeric matrix or data frame where rows represent metabolites and columns represent samples.
#' @param processing A character string specifying the processing method to use.
#' Options are `"shrink"` (default) or `"normalize"`.
#' @param control_samples An optional numeric vector specifying which columns correspond to control samples.
#' If provided, the reference spectrum is calculated using only these samples.
#' @param return_diagnostics Logical; if `TRUE`, return the corrected matrix and
#' sample-level dilution diagnostics in a list.
#' @return A numeric matrix of normalized intensities, or a list containing the
#' matrix and diagnostics when `return_diagnostics = TRUE`.
#' @examples
#' # Example usage:
#' your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
#' normalized_data <- normalize_by_dilution_factor(your_data_matrix, control_samples = 1:4)
#' @export
normalize_by_dilution_factor <- function(data,
                                         processing = "shrink",
                                         control_samples = NULL,
                                         return_diagnostics = FALSE) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame.")
  }
  if (is.data.frame(data))
    data <- as.matrix(data)
  .validate_intensity_matrix(data, require_nonnegative = TRUE)
  processing <- match.arg(processing, c("shrink", "normalize"))
  .validate_control_samples(control_samples, ncol(data))
  .validate_flag(return_diagnostics, "return_diagnostics")

  if (!is.null(control_samples)) {
    # Calculate reference spectrum using control samples only
    reference_spectrum <- apply(data[, control_samples, drop = FALSE], 1, median)
  } else {
    reference_spectrum <- apply(data, 1, median)
  }

  # Features with a zero reference cannot define a quotient. Zeros remain valid
  # measurements; only zero-reference features are omitted from factor fitting.
  reference_eligible <- is.finite(reference_spectrum) & reference_spectrum > 0
  if (!any(reference_eligible)) {
    stop(
      "No feature has a positive reference median; dilution factors cannot be estimated.",
      call. = FALSE
    )
  }
  quotients <- data[reference_eligible, , drop = FALSE] /
    reference_spectrum[reference_eligible]
  dilution_factors <- apply(quotients, 2, median)
  if (any(!is.finite(dilution_factors) | dilution_factors <= 0)) {
    stop(
      "Dilution factors must be finite and positive. Check zero-heavy samples ",
      "and declared non-detection values before normalization.",
      call. = FALSE
    )
  }

  mean_dilution <- mean(dilution_factors)
  stdev_dilution <- sd(dilution_factors)
  if (!is.finite(stdev_dilution)) {
    if (length(dilution_factors) == 1L) {
      stdev_dilution <- 0
    } else {
      stop("Dilution-factor standard deviation is not finite.", call. = FALSE)
    }
  }
  lower_threshold <- mean_dilution - stdev_dilution
  upper_threshold <- mean_dilution + stdev_dilution
  
  normalized_data <- data
  scaling_factors <- rep(1, length(dilution_factors))
  if (processing == "shrink") {
    for (i in seq_along(dilution_factors)) {
      if (dilution_factors[i] < lower_threshold) {
        scaling_factors[i] <- lower_threshold / dilution_factors[i]
        normalized_data[, i] <- (normalized_data[, i] / dilution_factors[i]) *
          lower_threshold
      }
      if (dilution_factors[i] > upper_threshold) {
        scaling_factors[i] <- upper_threshold / dilution_factors[i]
        normalized_data[, i] <- (normalized_data[, i] / dilution_factors[i]) *
          upper_threshold
      }
    }
  } else {
    for (i in seq_along(dilution_factors)) {
      scaling_factors[i] <- 1 / dilution_factors[i]
      normalized_data[, i] <- normalized_data[, i] / dilution_factors[i]
    }
  }
  if (!isTRUE(return_diagnostics)) {
    return(normalized_data)
  }

  sample_ids <- colnames(data)
  if (is.null(sample_ids)) {
    sample_ids <- paste0("sample_", seq_len(ncol(data)))
  }
  changed <- vapply(seq_len(ncol(data)), function(i) {
    any(abs(normalized_data[, i] - data[, i]) > 1e-12, na.rm = TRUE)
  }, logical(1))
  diagnostics <- data.frame(
    sample_id = sample_ids,
    dilution_factor = as.numeric(dilution_factors),
    lower_threshold = lower_threshold,
    upper_threshold = upper_threshold,
    altered = changed,
    scaling_factor = as.numeric(scaling_factors),
    reference_spectrum_type = if (is.null(control_samples)) "all_samples_median" else "control_samples_median",
    reference_features_used = sum(reference_eligible),
    zero_reference_features_omitted = sum(!reference_eligible),
    processing = processing,
    stringsAsFactors = FALSE
  )
  list(data = normalized_data, diagnostics = diagnostics)
}

#' Adjust Outliers in Data Matrix Using MAD
#'
#' This function adjusts outliers in a data matrix on a per-metabolite basis
#' using the median absolute deviation (MAD). For each metabolite, the median,
#' MAD, thresholds, and upper/lower outlier masks are computed once from the
#' original finite values. Upper- and lower-tail values are then shrunk
#' independently in one non-iterative pass; adjusting one tail cannot change
#' the threshold or membership of the other tail.
#' Non-finite values are not eligible for threshold estimation and are returned
#' unchanged by this standalone helper. The full [winn()] pipeline instead
#' requires finite, non-negative input. Numeric zeros are treated as observed
#' values and are never automatically converted to missing values.
#'
#' @param data A numeric matrix or data frame where rows represent metabolites and columns represent samples.
#' @return A numeric matrix with adjusted outliers.
#' @examples
#' your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
#' adjusted_data <- adjust_outliers_mad(your_data_matrix)
#' @export
adjust_outliers_mad <- function(data) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame.")
  }
  if (is.data.frame(data))
    data <- as.matrix(data)
  .validate_intensity_matrix(data, allow_nonfinite = TRUE)

  adjusted_data <- data
  for (i in seq_len(nrow(data))) {
    original <- data[i, ]
    eligible <- is.finite(original)
    if (!any(eligible)) {
      next
    }

    med <- median(original[eligible])
    mad_val <- mad(original[eligible])
    if (!is.finite(mad_val) || mad_val <= 0) {
      next
    }

    lower_threshold <- med - 4 * mad_val
    upper_threshold <- med + 4 * mad_val
    upper_outlier <- eligible & original >= upper_threshold
    lower_outlier <- eligible & original <= lower_threshold
    central <- eligible & !upper_outlier & !lower_outlier

    if (any(upper_outlier) && any(central)) {
      upper_reference <- max(original[central])
      upper_extreme <- max(original[upper_outlier])
      adjusted_data[i, upper_outlier] <- approx(
        x = c(upper_reference, upper_extreme),
        y = c(med + 3 * mad_val, upper_threshold),
        xout = original[upper_outlier],
        rule = 2
      )$y
    }

    if (any(lower_outlier) && any(central)) {
      lower_extreme <- min(original[lower_outlier])
      lower_reference <- min(original[central])
      adjusted_data[i, lower_outlier] <- approx(
        x = c(lower_extreme, lower_reference),
        y = c(lower_threshold, med - 3 * mad_val),
        xout = original[lower_outlier],
        rule = 2
      )$y
    }
  }
  adjusted_data
}

# ---- internal: correct inverse for log1p ----
.inv_log1p <- function(z, clamp_nonneg = TRUE) {
  x <- expm1(z)
  if (clamp_nonneg) {
    x <- pmax(x, 0)
  }
  x
}

.compute_drift_offsets <- function(segment, seg_run, spline_method) {
  apply(segment, 1, function(x) {
    log_seg <- log1p(x)
    if (spline_method == "conservative") {
      fit <- .fit_conservative_spline(log_seg, seg_run)
    } else {
      fit <- tryCatch(
        mgcv::gam(log_seg ~ mgcv::s(seg_run, bs = "cr")),
        error = function(e)
          NULL
      )
    }
    if (is.null(fit)) {
      linear_fit <- lm(log_seg ~ seg_run)
      linear_fit$fitted.values - mean(linear_fit$fitted.values)
    } else {
      fit$fitted.values - mean(fit$fitted.values)
    }
  })
}

.detrend_segment_rows <- function(segment, seg_run, spline_method) {
  if (!nrow(segment)) {
    return(segment)
  }
  spline_vals <- .compute_drift_offsets(segment, seg_run, spline_method)
  z <- log1p(segment) - t(spline_vals)
  .inv_log1p(z, clamp_nonneg = TRUE)
}

.compute_ljung_box_pvalues <- function(segment, segment_lag, model_df) {
  if (is.na(segment_lag)) {
    return(rep(NA_real_, nrow(segment)))
  }
  log_segment <- log1p(segment)
  if (anyNA(log_segment)) {
    return(apply(log_segment, 1, function(x_log) {
      tryCatch(
        Box.test(
          x_log,
          lag = segment_lag,
          type = "Ljung-Box",
          fitdf = model_df
        )$p.value,
        error = function(e)
          NA_real_
      )
    }))
  }
  centered <- sweep(log_segment, 1, rowMeans(log_segment), "-")
  denom <- rowSums(centered^2)
  valid <- denom > .Machine$double.eps
  if (!any(valid)) {
    return(rep(NA_real_, nrow(segment)))
  }
  n_obs <- ncol(centered)
  q_stat <- numeric(nrow(centered))
  for (lag_idx in seq_len(segment_lag)) {
    numer <- rowSums(
      centered[, (lag_idx + 1):n_obs, drop = FALSE] *
        centered[, seq_len(n_obs - lag_idx), drop = FALSE]
    )
    rho_sq <- numeric(length(numer))
    rho_sq[valid] <- (numer[valid] / denom[valid])^2
    q_stat <- q_stat + rho_sq / (n_obs - lag_idx)
  }
  q_stat <- n_obs * (n_obs + 2) * q_stat
  p_vals <- rep(NA_real_, nrow(centered))
  p_vals[valid] <- stats::pchisq(
    q_stat[valid],
    df = segment_lag - model_df,
    lower.tail = FALSE
  )
  p_vals
}

.compute_segment_autocorrelation_pvalues <- function(segment,
                                                     seg_run,
                                                     test,
                                                     segment_lag,
                                                     model_df) {
  if (identical(test, "Ljung-Box") && is.na(segment_lag)) {
    return(rep(NA_real_, nrow(segment)))
  }
  if (identical(test, "DW") && !requireNamespace("lmtest", quietly = TRUE)) {
    stop("Package 'lmtest' is required for DW test. Please install it.")
  }
  if (identical(test, "Ljung-Box")) {
    return(p.adjust(
      .compute_ljung_box_pvalues(
        segment = segment,
        segment_lag = segment_lag,
        model_df = model_df
      ),
      method = "fdr"
    ))
  }
  p_vals <- apply(segment, 1, function(x) {
    x_log <- log1p(x)
    if (test == "DW") {
      tryCatch(
        lmtest::dwtest(x_log ~ seg_run)$p.value,
        error = function(e)
          NA_real_
      )
    } else {
      stop("Invalid test method. Use 'Ljung-Box' or 'DW'.")
    }
  })
  p.adjust(p_vals, method = "fdr")
}

.prepare_autocorrelation_correction <- function(data,
                                                run_order,
                                                batch,
                                                lag,
                                                test,
                                                ljung_box_fitdf,
                                                detrend,
                                                spline_method,
                                                max_fdr_threshold) {
  batch_vec <- batch
  unique_batch <- unique(batch_vec)
  prepared <- list(
    data = data,
    mode = detrend,
    batches = vector("list", length(unique_batch))
  )
  for (i in seq_along(unique_batch)) {
    b <- unique_batch[i]
    idx <- .ordered_batch_indices(batch_vec, b, run_order)
    segment <- data[, idx, drop = FALSE]
    seg_run <- if (!is.null(run_order)) run_order[idx] else seq_along(idx)
    batch_info <- list(idx = idx)
    if (detrend == "spline") {
      batch_info$detrended <- .detrend_segment_rows(segment, seg_run, spline_method)
    } else if (detrend == "mean") {
      segment_n <- length(idx)
      model_df <- .autocorrelation_model_df(
        test = test,
        ljung_box_fitdf = ljung_box_fitdf
      )
      segment_lag <- .resolve_autocorrelation_lag(
        lag = lag,
        n_obs = segment_n,
        model_df = model_df
      )
      adjusted_p <- .compute_segment_autocorrelation_pvalues(
        segment = segment,
        seg_run = seg_run,
        test = test,
        segment_lag = segment_lag,
        model_df = model_df
      )
      candidate_ids <- which(adjusted_p < max_fdr_threshold)
      batch_info$adjusted_p <- adjusted_p
      batch_info$candidate_ids <- candidate_ids
      batch_info$detrended <- if (length(candidate_ids) > 0) {
        .detrend_segment_rows(
          segment[candidate_ids, , drop = FALSE],
          seg_run = seg_run,
          spline_method = spline_method
        )
      } else {
        segment[candidate_ids, , drop = FALSE]
      }
    } else {
      stop("Invalid detrend method. Use 'mean' or 'spline'.")
    }
    prepared$batches[[i]] <- batch_info
  }
  prepared
}

.materialize_autocorrelation_correction <- function(prepared, fdr_threshold = NULL) {
  corrected <- prepared$data
  for (batch_info in prepared$batches) {
    idx <- batch_info$idx
    if (prepared$mode == "spline") {
      corrected[, idx] <- batch_info$detrended
      next
    }
    correct_ids <- which(batch_info$adjusted_p < fdr_threshold)
    if (!length(correct_ids)) {
      next
    }
    candidate_pos <- match(correct_ids, batch_info$candidate_ids, nomatch = 0L)
    if (any(candidate_pos == 0L)) {
      stop("Autocorrelation cache does not cover the requested threshold.")
    }
    corrected[correct_ids, idx] <- batch_info$detrended[candidate_pos, , drop = FALSE]
  }
  corrected
}

.materialize_autocorrelation_candidates <- function(prepared, thresholds) {
  if (prepared$mode != "mean") {
    return(list(list(
      threshold = NA_real_,
      thresholds = NA_real_,
      matrix = .materialize_autocorrelation_correction(prepared, fdr_threshold = 1)
    )))
  }
  ascending_thresholds <- sort(unique(thresholds))
  current <- prepared$data
  threshold_groups <- list()
  matrices <- list()
  sorted_group_ids <- integer(length(ascending_thresholds))
  prev_threshold <- 0
  last_group_id <- 0L
  
  for (i in seq_along(ascending_thresholds)) {
    threshold_value <- ascending_thresholds[i]
    added_any <- FALSE
    for (batch_info in prepared$batches) {
      if (!length(batch_info$candidate_ids)) {
        next
      }
      candidate_p <- batch_info$adjusted_p[batch_info$candidate_ids]
      new_positions <- which(candidate_p >= prev_threshold & candidate_p < threshold_value)
      if (!length(new_positions)) {
        next
      }
      current[
        batch_info$candidate_ids[new_positions],
        batch_info$idx
      ] <- batch_info$detrended[new_positions, , drop = FALSE]
      added_any <- TRUE
    }
    if (!added_any && last_group_id > 0L) {
      sorted_group_ids[i] <- last_group_id
      threshold_groups[[last_group_id]] <- c(threshold_groups[[last_group_id]], threshold_value)
    } else {
      last_group_id <- last_group_id + 1L
      sorted_group_ids[i] <- last_group_id
      threshold_groups[[last_group_id]] <- threshold_value
      matrices[[last_group_id]] <- current
    }
    prev_threshold <- threshold_value
  }
  
  candidates <- list()
  seen_groups <- integer(0)
  for (threshold_value in thresholds) {
    group_id <- sorted_group_ids[match(threshold_value, ascending_thresholds)]
    if (group_id %in% seen_groups) {
      next
    }
    seen_groups <- c(seen_groups, group_id)
    candidates[[length(candidates) + 1L]] <- list(
      threshold = threshold_value,
      thresholds = threshold_groups[[group_id]],
      matrix = matrices[[group_id]]
    )
  }
  candidates
}

#' Correct for Drift in Data Using Autocorrelation Correction
#'
#' This function corrects for drift effects in metabolomics data by detrending
#' based on run order within each batch segment. When `run_order` is supplied,
#' every segment is sorted by that vector for both testing and smoothing, and
#' the corrected values are then restored to the original matrix-column order.
#'
#' @param data A numeric matrix with rows representing metabolites and columns representing samples.
#' @param run_order An optional numeric vector representing the run order of the samples.
#' If `NULL`, the current matrix-column order is treated as acquisition order.
#' @param batch A vector indicating the batch (or segment) assignment for each sample.
#' @param lag An optional integer specifying the lag to be used in the autocorrelation test.
#' If `NULL`, the lag is selected adaptively for each batch segment using
#' `max(min(10, floor(n / 5)), df + 3)` and then capped at `n - 1`.
#' @param test A character string specifying the autocorrelation test to use ("Ljung-Box" or "DW").
#' @param ljung_box_fitdf Non-negative integer passed to [stats::Box.test()] as
#' `fitdf`. The default is `0` because WiNN does not fit an ARMA model before
#' testing. Set this to `1` only to reproduce the legacy WiNN gate or for a
#' documented sensitivity analysis.
#' @param detrend A character string indicating the method for detrending ("mean" or "spline").
#' @param fdr_threshold A numeric value specifying the FDR threshold for significance.
#' @param spline_method A character string specifying the spline method when detrend="spline" ("conservative" or "standard").
#' @param gate Optional gate mode: `"selective"`, `"all"`, or `"none"`.
#' When `NULL`, the historical behavior is retained: `detrend = "mean"` is
#' selective and `detrend = "spline"` corrects all testable profiles.
#' @param return_diagnostics Logical; if `TRUE`, return profile-level test and
#' correction diagnostics together with the corrected matrix.
#' @return A numeric matrix with drift corrected, or a list containing the
#' matrix and diagnostics when `return_diagnostics = TRUE`.
#' @examples
#' your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
#' batch <- rep(1:4, length.out = ncol(your_data_matrix))
#' run_order <- seq_len(ncol(your_data_matrix))
#' drift_corrected <- autocorrelation_correct(your_data_matrix, run_order, batch)
#' @export
autocorrelation_correct <- function(data,
                                    run_order = NULL,
                                    batch,
                                    lag = NULL,
                                    test = "Ljung-Box",
                                    ljung_box_fitdf = 0L,
                                    detrend = "mean",
                                    fdr_threshold = 0.05,
                                    spline_method = "conservative",
                                    gate = NULL,
                                    return_diagnostics = FALSE) {
  .validate_intensity_matrix(data, require_nonnegative = TRUE)
  .validate_batch(batch, ncol(data))
  .validate_run_order(run_order, batch, ncol(data))
  .validate_probability(fdr_threshold, "fdr_threshold")
  .validate_flag(return_diagnostics, "return_diagnostics")
  ljung_box_fitdf <- .validate_ljung_box_fitdf(ljung_box_fitdf)
  if (!identical(test, "Ljung-Box") && !identical(test, "DW")) {
    stop("test must be either 'Ljung-Box' or 'DW'.")
  }
  if (!is.null(lag) && (length(lag) != 1 || !is.numeric(lag) || is.na(lag) ||
    lag <= 0 || lag != round(lag))) {
    stop("lag must be NULL or a positive integer.")
  }
  if (!detrend %in% c("mean", "spline")) {
    stop("detrend must be either 'mean' or 'spline'.")
  }
  if (!spline_method %in% c("conservative", "standard")) {
    stop("spline_method must be either 'conservative' or 'standard'.")
  }
  legacy_all <- is.null(gate) && identical(detrend, "spline")
  if (is.null(gate)) {
    gate <- if (identical(detrend, "spline")) "all" else "selective"
  }
  gate <- match.arg(gate, c("selective", "all", "none"))

  detrended_data <- data
  batch_vec <- batch
  unique_batch <- unique(batch_vec)
  feature_ids <- rownames(data)
  if (is.null(feature_ids)) {
    feature_ids <- paste0("feature_", seq_len(nrow(data)))
  }
  diagnostic_rows <- vector("list", length(unique_batch))

  for (batch_index in seq_along(unique_batch)) {
    b <- unique_batch[[batch_index]]
    idx <- .ordered_batch_indices(batch_vec, b, run_order)
    segment <- data[, idx, drop = FALSE]
    segment_n <- length(idx)
    model_df <- .autocorrelation_model_df(
      test = test,
      ljung_box_fitdf = ljung_box_fitdf
    )
    segment_lag <- .resolve_autocorrelation_lag(
      lag = lag,
      n_obs = segment_n,
      model_df = model_df
    )
    seg_run <- if (!is.null(run_order)) run_order[idx] else seq_along(idx)
    p_values <- rep(NA_real_, nrow(segment))
    if (!identical(gate, "none")) {
      p_values <- apply(segment, 1, function(x) {
        x_log <- log1p(x)
        if (test == "Ljung-Box") {
          if (is.na(segment_lag)) {
            return(NA_real_)
          }
          tryCatch(
            Box.test(
              x_log,
              lag = segment_lag,
              type = "Ljung-Box",
              fitdf = model_df
            )$p.value,
            error = function(e) NA_real_
          )
        } else if (test == "DW") {
          if (!requireNamespace("lmtest", quietly = TRUE)) {
            stop("Package 'lmtest' is required for DW test. Please install it.")
          }
          tryCatch(
            lmtest::dwtest(x_log ~ seg_run)$p.value,
            error = function(e) NA_real_
          )
        } else {
          stop("Invalid test method. Use 'Ljung-Box' or 'DW'.")
        }
      })
    }
    adjusted_values <- p.adjust(p_values, method = "fdr")
    testable <- is.finite(p_values)
    selective_decision <- testable & is.finite(adjusted_values) &
      adjusted_values < fdr_threshold
    requested_correction <- switch(
      gate,
      selective = selective_decision,
      all = if (isTRUE(legacy_all)) rep(TRUE, nrow(segment)) else testable,
      none = rep(FALSE, nrow(segment))
    )

    actual_correction <- rep(FALSE, nrow(segment))
    decision_reason <- ifelse(
      !testable & !identical(gate, "none") & !isTRUE(legacy_all),
      "not_testable",
      ifelse(requested_correction,
        if (identical(gate, "all")) "forced_all" else "selected",
        "not_selected"
      )
    )
    smoother_type <- rep(NA_character_, nrow(segment))
    basis_dimension <- rep(NA_real_, nrow(segment))
    effective_df <- rep(NA_real_, nrow(segment))
    fallback_path <- rep(NA_character_, nrow(segment))
    warning_status <- rep(FALSE, nrow(segment))
    warning_message <- rep("", nrow(segment))
    correction_magnitude <- rep(0, nrow(segment))

    correct_ids <- which(requested_correction)
    for (feature_index in correct_ids) {
      log_seg <- log1p(segment[feature_index, ])
      fit <- if (identical(spline_method, "conservative")) {
        .fit_conservative_spline(log_seg, seg_run)
      } else {
        standard_fit <- tryCatch(
          mgcv::gam(log_seg ~ mgcv::s(seg_run, bs = "cr")),
          error = function(e) NULL
        )
        if (is.null(standard_fit)) {
          standard_fit <- tryCatch(lm(log_seg ~ seg_run), error = function(e) NULL)
          if (!is.null(standard_fit)) {
            standard_fit <- .tag_spline_fit(
              standard_fit,
              smoother_type = "linear",
              fallback_path = "standard_gam->linear",
              basis_dimension = 2
            )
          }
        } else {
          standard_basis_dimension <- tryCatch(
            standard_fit$smooth[[1]]$bs.dim,
            error = function(e) NA_real_
          )
          standard_fit <- .tag_spline_fit(
            standard_fit,
            smoother_type = "standard_gam_cr",
            fallback_path = "standard_gam",
            basis_dimension = standard_basis_dimension
          )
        }
        standard_fit
      }

      if (is.null(fit) || is.null(fit$fitted.values) ||
          any(!is.finite(fit$fitted.values))) {
        decision_reason[feature_index] <- "smoother_failed"
        next
      }
      centered_fit <- fit$fitted.values - mean(fit$fitted.values)
      corrected_log <- log_seg - centered_fit
      corrected_values <- .inv_log1p(corrected_log, clamp_nonneg = TRUE)
      detrended_data[feature_index, idx] <- corrected_values
      actual_correction[feature_index] <- TRUE
      fit_smoother_type <- attr(fit, "winn_smoother_type")
      if (length(fit_smoother_type) == 1L) {
        smoother_type[feature_index] <- fit_smoother_type
      }
      fit_basis_dimension <- attr(fit, "winn_basis_dimension")
      if (length(fit_basis_dimension) == 1L) {
        basis_dimension[feature_index] <- fit_basis_dimension
      }
      fit_fallback_path <- attr(fit, "winn_fallback_path")
      if (length(fit_fallback_path) == 1L) {
        fallback_path[feature_index] <- fit_fallback_path
      }
      fit_warnings <- attr(fit, "winn_warnings")
      if (!is.null(fit_warnings) && length(fit_warnings)) {
        warning_status[feature_index] <- TRUE
        warning_message[feature_index] <- paste(fit_warnings, collapse = " | ")
      }
      if (inherits(fit, "gam") && !is.null(fit$edf)) {
        effective_df[feature_index] <- sum(fit$edf)
      } else if (!is.null(fit$rank)) {
        effective_df[feature_index] <- fit$rank
      }
      correction_magnitude[feature_index] <- median(
        abs(log1p(corrected_values) - log_seg),
        na.rm = TRUE
      )
    }

    diagnostic_rows[[batch_index]] <- data.frame(
      feature_id = feature_ids,
      segment = as.character(b),
      segment_size = segment_n,
      testable = testable,
      lag = if (identical(test, "Ljung-Box")) segment_lag else NA_integer_,
      ljung_box_fitdf = if (identical(test, "Ljung-Box")) model_df else NA_integer_,
      raw_p_value = as.numeric(p_values),
      adjusted_p_value = as.numeric(adjusted_values),
      selective_significance = selective_decision,
      actual_correction = actual_correction,
      decision_reason = decision_reason,
      gate = gate,
      smoother_type = smoother_type,
      basis_dimension = basis_dimension,
      effective_df = effective_df,
      fallback_path = fallback_path,
      warning_status = warning_status,
      warning_message = warning_message,
      correction_magnitude = correction_magnitude,
      stringsAsFactors = FALSE
    )
  }
  if (!isTRUE(return_diagnostics)) {
    return(detrended_data)
  }
  list(data = detrended_data, diagnostics = do.call(rbind, diagnostic_rows))
}

#' Perform ANOVA-based Mean-Only Batch Correction
#'
#' This function runs an ANOVA test on each metabolite to detect batch effects,
#' and then corrects significant batch effects by subtracting the estimated
#' batch-specific shifts (while preserving the overall mean). With the default
#' selective gate, features that do not pass the FDR threshold are unchanged;
#' this differs from [combat_batch_correction()], which applies an empirical
#' Bayes location/scale model to every feature.
#'
#' @param data A numeric matrix (metabolites × samples).
#' @param batch A factor or numeric vector indicating batch for each sample.
#' @param fdr_threshold Significance threshold for FDR-adjusted p-values.
#' @param gate Batch gate mode: `"selective"` (default), `"all"`, or
#' `"none"`.
#' @param return_diagnostics Logical; if `TRUE`, return feature-level test and
#' correction diagnostics together with the corrected matrix.
#' @return A numeric matrix of corrected intensities, or a list containing the
#' matrix and diagnostics when `return_diagnostics = TRUE`.
#' @examples
#' mat <- matrix(rnorm(200, mean = 100, sd = 15), nrow=20)
#' batch <- rep(1:4, length.out=ncol(mat))
#' corrected <- anova_batch_correction(mat, batch, fdr_threshold=0.05)
#' @export
anova_batch_correction <- function(data,
                                   batch,
                                   fdr_threshold = 0.05,
                                   gate = c("selective", "all", "none"),
                                   return_diagnostics = FALSE) {
  .validate_intensity_matrix(data, require_nonnegative = TRUE)
  .validate_batch(batch, ncol(data))
  .validate_probability(fdr_threshold, "fdr_threshold")
  .validate_flag(return_diagnostics, "return_diagnostics")
  gate <- match.arg(gate)
  batch <- factor(batch)
  feature_ids <- rownames(data)
  if (is.null(feature_ids)) {
    feature_ids <- paste0("feature_", seq_len(nrow(data)))
  }
  if (identical(gate, "none")) {
    diagnostics <- data.frame(
      feature_id = feature_ids,
      raw_p_value = NA_real_,
      adjusted_p_value = NA_real_,
      selective_significance = FALSE,
      eligible = apply(data, 1, function(v) all(is.finite(v))),
      actual_correction = FALSE,
      decision_reason = "not_selected",
      n_batches = length(unique(batch)),
      gate = gate,
      correction_magnitude = 0,
      stringsAsFactors = FALSE
    )
    if (!isTRUE(return_diagnostics)) {
      return(data)
    }
    return(list(data = data, diagnostics = diagnostics))
  }
  if (length(unique(batch)) < 2) {
    message("Only one batch detected. Skipping ANOVA-based correction.")
    if (!isTRUE(return_diagnostics)) {
      return(data)
    }
    return(list(
      data = data,
      diagnostics = data.frame(
        feature_id = feature_ids,
        raw_p_value = NA_real_,
        adjusted_p_value = NA_real_,
        selective_significance = FALSE,
        eligible = FALSE,
        actual_correction = FALSE,
        decision_reason = "not_testable",
        n_batches = 1L,
        gate = gate,
        correction_magnitude = 0,
        stringsAsFactors = FALSE
      )
    ))
  }
  z <- log1p(data)
  n_met <- nrow(z)
  pvals <- rep(NA_real_, n_met)
  for (i in seq_len(n_met)) {
    df <- data.frame(value = z[i, ], batch = batch)
    pvals[i] <- tryCatch({
      a <- aov(value ~ batch, data = df)
      summary(a)[[1]]$`Pr(>F)`[1]
    }, error = function(e) NA_real_)
  }
  padj <- p.adjust(pvals, method = "fdr")

  corrected <- z
  eligible <- apply(z, 1, function(v) all(is.finite(v)))
  selective_decision <- is.finite(padj) & padj < fdr_threshold
  correct_mask <- switch(
    gate,
    selective = selective_decision,
    all = eligible,
    none = rep(FALSE, n_met)
  )
  sig <- which(correct_mask)
  if (length(sig) > 0L) {
    overall_means <- rowMeans(z, na.rm = TRUE)
    for (i in sig) {
      # compute batch-specific and overall means
      batch_means <- tapply(z[i, ], batch, mean, na.rm = TRUE)
      shift <- batch_means - overall_means[i]
      # subtract the shift for each sample
      corrected[i, ] <- z[i, ] - shift[as.character(batch)]
    }
  }
  corrected_intensity <- .inv_log1p(corrected, clamp_nonneg = TRUE)
  if (!isTRUE(return_diagnostics)) {
    return(corrected_intensity)
  }
  correction_magnitude <- apply(abs(corrected - z), 1, median, na.rm = TRUE)
  decision_reason <- ifelse(
    !eligible,
    "not_testable",
    ifelse(
      correct_mask,
      if (identical(gate, "all")) "forced_all" else "selected",
      "not_selected"
    )
  )
  list(
    data = corrected_intensity,
    diagnostics = data.frame(
      feature_id = feature_ids,
      raw_p_value = as.numeric(pvals),
      adjusted_p_value = as.numeric(padj),
      selective_significance = selective_decision,
      eligible = eligible,
      actual_correction = correct_mask,
      decision_reason = decision_reason,
      n_batches = nlevels(batch),
      gate = gate,
      correction_magnitude = as.numeric(correction_magnitude),
      stringsAsFactors = FALSE
    )
  )
}

#' Perform ComBat Batch Correction by batch
#'
#' This function applies the empirical Bayes ComBat method to correct batch effects
#' by batch, adjusting both location and scale parameters across batch. Unlike
#' the default selective ANOVA gate in [anova_batch_correction()], ComBat is a
#' global correction: every feature is passed to the empirical Bayes model.
#'
#' @param data A numeric matrix (metabolites × samples).
#' @param batch A factor or numeric vector indicating batch for each sample.
#' @param par_prior Logical indicating whether to use parametric prior (default TRUE).
#' @param mean_only Logical indicating mean-only adjustment (default FALSE).
#' @param ref_batch Optional reference batch level for anchoring (default NULL).
#' @return A numeric matrix of corrected intensities.
#' @examples
#' mat <- matrix(rnorm(200, mean = 100, sd = 15), nrow=20)
#' batch <- rep(1:4, length.out=ncol(mat))
#' corrected <- combat_batch_correction(mat, batch, par_prior = TRUE)
#' @export
combat_batch_correction <- function(data,
                                    batch,
                                    par_prior = TRUE,
                                    mean_only = FALSE,
                                    ref_batch = NULL) {
  .validate_intensity_matrix(data, require_nonnegative = TRUE)
  .validate_batch(batch, ncol(data))
  .validate_flag(par_prior, "par_prior")
  .validate_flag(mean_only, "mean_only")
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for ComBat correction. Please install it.")
  }
  batch <- factor(batch)
  if (!is.null(ref_batch) &&
      (length(ref_batch) != 1L || is.na(ref_batch) ||
       !as.character(ref_batch) %in% levels(batch))) {
    stop("ref_batch must identify exactly one observed batch level.", call. = FALSE)
  }
  if (length(unique(batch)) < 2) {
    message("Only one batch detected. Skipping ComBat-based correction.")
    return(data)
  }
  z <- log1p(data)
  # design matrix with intercept only to preserve global mean
  mod <- model.matrix( ~ 1, data = data.frame(batch = batch))
  
  # apply ComBat
  corrected <- sva::ComBat(
    dat = z,
    batch = batch,
    mod = mod,
    par.prior = par_prior,
    mean.only = mean_only,
    ref.batch = ref_batch
  )
  return(.inv_log1p(corrected, clamp_nonneg = TRUE))
}


#' Scale Data by batch
#'
#' This function scales the values for each metabolite within each batch by subtracting the batch mean and dividing by the batch standard deviation.
#'
#' @param data A numeric matrix with rows representing metabolites and columns representing samples.
#' @param batch A numeric vector indicating the batch (or segment) assignment for each sample.
#' @return A numeric matrix of scaled intensities.
#' @examples
#' your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
#' batch <- rep(1:4, length.out = ncol(your_data_matrix))
#' scaled_data <- scale_by_batch(your_data_matrix, batch)
#' @export
scale_by_batch <- function(data, batch) {
  .validate_intensity_matrix(data)
  .validate_batch(batch, ncol(data))
  
  scaled_data <- data
  batch_vec <- batch
  unique_batch <- unique(batch_vec)
  for (b in unique_batch) {
    idx <- which(batch_vec == b)
    if (length(idx) < 2L) {
      stop(
        "scale_by_batch() requires at least two samples in every batch; batch '",
        as.character(b),
        "' has ",
        length(idx),
        ".",
        call. = FALSE
      )
    }
    row_means <- rowMeans(data[, idx, drop = FALSE])
    row_sds <- apply(data[, idx, drop = FALSE], 1, sd)
    if (any(!is.finite(row_sds))) {
      stop("Within-batch standard deviations must be finite.", call. = FALSE)
    }
    # A constant feature is centered but not divided by zero.
    row_sds[row_sds == 0] <- 1
    scaled_data[, idx] <- (data[, idx, drop = FALSE] - row_means) / row_sds
  }
  return(scaled_data)
}

#' Winn Correction for Metabolomics Data
#'
#' This function performs a series of corrections on metabolomics data to adjust for dilution effects, outliers, drift, and batch effects.
#' If batch information is not supplied, segments are automatically detected using an fkPELT-based approach and labeled as batch.
#'
#' The correction pipeline is as follows:
#' \enumerate{
#'   \item Outlier adjustment using MAD.
#'   \item Drift correction using autocorrelation-based detrending.
#'   \item Batch effect correction using ANOVA-based mean-adjustment or ComBat.
#'   \item Per-sample median adjustment using dilution factor.
#'   \item Optional scaling by batch.
#' }
#'
#' When control samples are provided and parameters="auto", the function performs comprehensive parameter
#' optimization by testing multiple combinations of settings and selecting those that maximize control sample
#' correlation while minimizing coefficient of variation. This searches spline
#' methods, FDR thresholds, and normalization approaches while respecting the
#' supplied autocorrelation-test family and `scale_by_batch` setting.
#' Auto mode requires at least two controls. In auto mode, `fdr_threshold`,
#' `spline_method`, and `median_adjustment` are fixed-mode settings and are
#' replaced by the documented internal candidate grids. Other arguments still
#' constrain the candidate family. Set `return_details = TRUE` to retain the
#' selected settings and stage decisions in machine-readable form.
#'
#' The full pipeline expects finite, non-negative quantitative intensities.
#' Zeros are treated as observed values and are not converted to missingness.
#' Users should recode declared non-detection sentinels and impute missing
#' values using a documented preprocessing policy before calling `winn()`.
#'
#' @param data A numeric matrix or data frame where rows represent metabolites and columns represent samples.
#' @param batch An optional numeric vector indicating batch assignments for each sample. If NULL, segments will be auto-detected.
#' @param run_order An optional numeric vector representing the run order of samples.
#' @param control_samples An optional numeric vector representing the columns corresponding to control samples. If provided,
#' these will be used for normalization and parameter tuning.
#' @param parameters A character string specifying fixed (`"fixed"`, the
#' default) or QC-tuned (`"auto"`) parameters. Auto mode requires controls.
#' @param fdr_threshold Fixed-mode FDR threshold for drift and selective ANOVA
#' batch correction. Auto mode searches 0.10, 0.05, and 0.01.
#' @param median_adjustment Fixed-mode median adjustment (`"shrink"`,
#' `"normalize"`, or `"none"`). Auto mode searches `"shrink"` and
#' `"normalize"`.
#' @param detrend_non_autocorrelated A character string specifying the method for detrending non-autocorrelated metabolites ("mean" or "spline").
#' @param spline_method Fixed-mode spline method (`"conservative"` or
#' `"standard"`). Auto mode searches both.
#' @param remove_batch_effects Batch method. `"anova"` selectively corrects
#' features that pass its FDR gate; `"combat"` applies global empirical Bayes
#' location/scale correction to all features.
#' @param test A character vector specifying the autocorrelation test(s)
#' to use. Fixed mode requires a single value. Auto mode evaluates only the
#' supplied tests, so `DW` is included only when explicitly requested.
#' @param lag An optional integer specifying the lag for the autocorrelation test.
#' If `NULL`, `autocorrelation_correct()` selects the lag adaptively for each
#' batch segment.
#' @param ljung_box_fitdf Non-negative integer passed to [stats::Box.test()].
#' The default is `0` because no ARMA model is fitted. Use `1` only for legacy
#' reproduction or a documented sensitivity analysis.
#' @param scale_by_batch Logical indicating whether to scale data by batch after corrections.
#' @param pelt_penalty Optional fkPELT penalty specification. Use a positive
#' numeric value, `"bic"`, `"mbic"`, or `NULL`. When `NULL`, fkPELT uses
#' the conservative MBIC default.
#' @param return_details Logical. The default `FALSE` returns the corrected
#' matrix for backward compatibility. `TRUE` returns a `winn_result` list with
#' the matrix, selected parameters, batch assignments, and stage decisions.
#' @return A numeric matrix of corrected intensities, or a structured
#' `winn_result` when `return_details = TRUE`.
#' @examples
#' your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
#' batch <- rep(1:4, length.out = ncol(your_data_matrix))
#' run_order <- seq_len(ncol(your_data_matrix))
#' corrected_data <- winn(your_data_matrix, batch = batch, run_order = run_order)
#' @export
winn <- function(data,
                 batch = NULL,
                 run_order = NULL,
                 control_samples = NULL,
                 parameters = "fixed",
                 fdr_threshold = 0.05,
                 median_adjustment = "shrink",
                 detrend_non_autocorrelated = "mean",
                 spline_method = "conservative",
                 remove_batch_effects = "anova",
                 test = "Ljung-Box",
                 lag = NULL,
                 ljung_box_fitdf = 0L,
                 scale_by_batch = FALSE,
                 pelt_penalty = NULL,
                 return_details = FALSE) {
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }
  if (is.data.frame(data))
    data <- as.matrix(data)
  .validate_intensity_matrix(data, require_nonnegative = TRUE)
  
  n_samples <- ncol(data)
  n_metabolites <- nrow(data)
  
  if (n_samples < 3)
    stop("At least 3 samples required.")
  if (n_metabolites < 1)
    stop("At least 1 metabolite required.")
  
  # Parameter validation
  if (!is.character(parameters) || length(parameters) != 1L ||
      !parameters %in% c("fixed", "auto")) {
    stop("parameters must be either 'fixed' or 'auto'.")
  }
  .validate_batch(batch, n_samples, required = FALSE)
  .validate_run_order(run_order, batch, n_samples)
  .validate_control_samples(control_samples, n_samples)
  if (!is.null(control_samples)) {
    control_samples <- as.integer(control_samples)
  }
  if (identical(parameters, "auto") && is.null(control_samples)) {
    stop("parameters = 'auto' requires at least two control samples.")
  }
  if (identical(parameters, "auto") && length(control_samples) < 2L) {
    stop("At least two control samples are required for parameter optimization.")
  }
  if (!spline_method %in% c("conservative", "standard")) {
    stop("spline_method must be either 'conservative' or 'standard'.")
  }
  if (!detrend_non_autocorrelated %in% c("mean", "spline")) {
    stop("detrend_non_autocorrelated must be either 'mean' or 'spline'.")
  }
  if (!median_adjustment %in% c("shrink", "normalize", "none")) {
    stop("median_adjustment must be 'shrink', 'normalize', or 'none'.")
  }
  if (!remove_batch_effects %in% c("anova", "combat")) {
    stop("remove_batch_effects must be either 'anova' or 'combat'.")
  }
  if (!is.character(test) || length(test) < 1L || any(!test %in% c("Ljung-Box", "DW"))) {
    stop("test must contain only 'Ljung-Box' and/or 'DW'.")
  }
  if (parameters == "fixed" && length(test) != 1L) {
    stop("test must be length 1 when parameters = 'fixed'.")
  }
  
  .validate_probability(fdr_threshold, "fdr_threshold")
  ljung_box_fitdf <- .validate_ljung_box_fitdf(ljung_box_fitdf)
  .validate_flag(scale_by_batch, "scale_by_batch")
  .validate_flag(return_details, "return_details")
  if (!is.null(lag) && (length(lag) != 1 || !is.numeric(lag) || is.na(lag) ||
    lag <= 0 || lag != round(lag))) {
    stop("lag must be NULL or a positive integer.")
  }
  if (!is.null(pelt_penalty) &&
    !(.is_valid_pelt_penalty_value(pelt_penalty))) {
    stop(
      "pelt_penalty must be NULL, a single positive numeric value, ",
      "'bic', or 'mbic'."
    )
  }
  
  message("Starting Winn correction...")

  batch_source <- if (is.null(batch)) "auto_detected" else "supplied"
  stage_decisions <- list()
  selected_parameters <- list()

  # Outlier adjustment
  message("Adjusting outliers using MAD...")
  norm_data <- adjust_outliers_mad(data)
  if (isTRUE(return_details)) {
    feature_ids <- rownames(data)
    if (is.null(feature_ids)) {
      feature_ids <- paste0("feature_", seq_len(nrow(data)))
    }
    outlier_delta <- abs(log1p(norm_data) - log1p(data))
    stage_decisions$outlier <- data.frame(
      feature_id = feature_ids,
      entries_changed = rowSums(outlier_delta > 1e-12),
      median_absolute_log1p_change = apply(outlier_delta, 1, median),
      stringsAsFactors = FALSE
    )
  }
  
  # If control samples are provided and auto parameter detection is enabled, perform grid search
  if (!is.null(control_samples) && parameters == "auto") {
    message("Auto-detecting optimal parameters using control samples...")
    
    # Define grid for parameter search
    batch_options <- if (!is.null(batch))
      c("provided")
    else
      c("auto")
    tests <- unique(test)
    if (detrend_non_autocorrelated == "spline") {
      tests <- tests[1]
    }
    normalizations <- c("shrink", "normalize")
    acorr_fdr_options <- c(0.1, 0.05, 0.01)
    anova_fdr_options <- c(0.1, 0.05, 0.01)
    drift_fdr_grid <- if (detrend_non_autocorrelated == "mean") {
      acorr_fdr_options
    } else {
      NA_real_
    }
    batch_fdr_grid <- if (remove_batch_effects == "anova") {
      anova_fdr_options
    } else {
      NA_real_
    }
    scale_options <- if (scale_by_batch) {
      c(TRUE)
    } else {
      c(FALSE)
    }
    spline_methods <- c("conservative", "standard")  # Add spline method optimization
    
    best_score <- -Inf
    best_final_data <- NULL
    best_params <- list()
    best_batch <- NULL
    
    # Progress tracking for parameter optimization
    total_combinations <- length(batch_options) * length(tests) * length(spline_methods) *
      length(drift_fdr_grid) * length(batch_fdr_grid) *
      length(normalizations) * length(scale_options)
    current_combination <- 0
    
    for (batch_option in batch_options) {
      current_batch <- if (batch_option == "provided") {
        batch
      } else {
        .auto_detect_batch(
          norm_data,
          pelt_penalty = pelt_penalty,
          run_order = run_order
        )
      }
      
      for (current_test in tests) {
        for (current_spline_method in spline_methods) {
          prepared_drift <- tryCatch({
            .prepare_autocorrelation_correction(
              data = norm_data,
              run_order = run_order,
              batch = current_batch,
              lag = lag,
              test = current_test,
              ljung_box_fitdf = ljung_box_fitdf,
              detrend = detrend_non_autocorrelated,
              spline_method = current_spline_method,
              max_fdr_threshold = if (detrend_non_autocorrelated == "mean") {
                max(acorr_fdr_options)
              } else {
                1
              }
            )
          }, error = function(e) {
            NULL
          })
          
          if (is.null(prepared_drift))
            next
          
          drift_candidates <- .materialize_autocorrelation_candidates(
            prepared_drift,
            thresholds = if (prepared_drift$mode == "mean") {
              acorr_fdr_options
            } else {
              NA_real_
            }
          )
          
          for (drift_candidate in drift_candidates) {
            acorr_fdr <- drift_candidate$threshold
            drift_corrected <- drift_candidate$matrix
            
            if (remove_batch_effects == "anova") {
              prepared_batch <- tryCatch({
                .prepare_anova_batch_correction(drift_corrected, current_batch)
              }, error = function(e) {
                NULL
              })
              
              if (is.null(prepared_batch))
                next
              
              batch_candidates <- .materialize_anova_candidates(
                prepared_batch,
                thresholds = anova_fdr_options
              )
            } else {
              batch_candidates <- list(list(
                threshold = NA_real_,
                thresholds = NA_real_,
                matrix = tryCatch({
                  combat_batch_correction(drift_corrected, current_batch)
                }, error = function(e) {
                  NULL
                })
              ))
            }
            
            if (is.null(batch_candidates[[1]]$matrix))
              next
            
            for (batch_candidate in batch_candidates) {
              batch_corrected <- batch_candidate$matrix
              anova_fdr <- batch_candidate$threshold
              
              if (is.null(batch_corrected))
                next
              
              for (normalize_opts in normalizations) {
                normalized_data <- tryCatch({
                  normalize_by_dilution_factor(
                    batch_corrected,
                    processing = normalize_opts,
                    control_samples = control_samples
                  )
                }, error = function(e) {
                  NULL
                })
                
                if (is.null(normalized_data))
                  next
                
                for (scale_opt in scale_options) {
                  current_combination <- current_combination + 1
                  
                  if (current_combination %% max(1, floor(total_combinations / 10)) == 0) {
                    message(
                      "Parameter optimization progress: ",
                      round(
                        100 * current_combination / total_combinations,
                        1
                      ),
                      "%"
                    )
                  }
                  
                  final_data <- if (scale_opt) {
                    tryCatch({
                      scale_by_batch(normalized_data, current_batch)
                    }, error = function(e)
                      normalized_data)
                  } else {
                    normalized_data
                  }
                  
                  quality_metrics <- .calculate_quality_score(final_data, control_samples)
                  if (is.na(quality_metrics$score))
                    next
                  
                  if (quality_metrics$score > best_score) {
                    best_score <- quality_metrics$score
                    best_final_data <- final_data
                    best_batch <- current_batch
                    best_params <- list(
                      parameters = "auto",
                      batch_option = batch_option,
                      pelt_penalty = attr(current_batch, "pelt_penalty"),
                      test = current_test,
                      ljung_box_fitdf = if (identical(current_test, "Ljung-Box")) {
                        ljung_box_fitdf
                      } else {
                        "not used"
                      },
                      lag = if (is.null(lag) && identical(current_test, "Ljung-Box")) {
                        "adaptive"
                      } else if (identical(current_test, "Ljung-Box")) {
                        as.character(as.integer(lag))
                      } else {
                        "not used"
                      },
                      spline_method = current_spline_method,
                      acorr_fdr = if (prepared_drift$mode == "mean") {
                        acorr_fdr
                      } else {
                        "not used"
                      },
                      anova_fdr = if (remove_batch_effects == "anova") {
                        anova_fdr
                      } else {
                        "not used"
                      },
                      normalization = normalize_opts,
                      scale_by_batch = scale_opt,
                      quality_score = quality_metrics$score,
                      mean_cv = quality_metrics$mean_cv,
                      mean_correlation = quality_metrics$mean_correlation
                    )
                  }
                }
              }
            }
          }
        }
      }
    }
    
    if (is.null(best_final_data)) {
      stop("Parameter optimization failed. No valid parameter combination found.")
    }
    
    message("Optimal parameters selected based on control samples:")
    message("  Batch detection: ", best_params$batch_option)
    if (identical(best_params$batch_option, "auto") && !is.null(best_params$pelt_penalty)) {
      message("  fkPELT penalty: ", round(best_params$pelt_penalty, 4))
    }
    message("  Autocorrelation test: ",
            best_params$test,
            " (FDR: ",
            best_params$acorr_fdr,
            ")")
    if (identical(best_params$test, "Ljung-Box")) {
      message("  Ljung-Box lag: ", best_params$lag)
      message("  Ljung-Box fitdf: ", best_params$ljung_box_fitdf)
    }
    message("  Spline method: ", best_params$spline_method)
    message("  Batch correction FDR: ", best_params$anova_fdr)
    message("  Normalization: ", best_params$normalization)
    message("  Scale by batch: ", best_params$scale_by_batch)
    message(
      "  Quality metrics - CV: ",
      round(best_params$mean_cv, 4),
      ", Correlation: ",
      round(best_params$mean_correlation, 4)
    )
    
    final_data <- best_final_data
    batch <- best_batch
    selected_parameters <- best_params
    if (isTRUE(return_details)) {
      drift_fdr <- if (identical(detrend_non_autocorrelated, "mean")) {
        as.numeric(best_params$acorr_fdr)
      } else {
        fdr_threshold
      }
      drift_detail <- autocorrelation_correct(
        norm_data,
        run_order = run_order,
        batch = batch,
        lag = lag,
        test = best_params$test,
        ljung_box_fitdf = ljung_box_fitdf,
        detrend = detrend_non_autocorrelated,
        fdr_threshold = drift_fdr,
        spline_method = best_params$spline_method,
        return_diagnostics = TRUE
      )
      stage_decisions$drift <- drift_detail$diagnostics

      if (identical(remove_batch_effects, "anova")) {
        batch_detail <- anova_batch_correction(
          drift_detail$data,
          batch,
          fdr_threshold = as.numeric(best_params$anova_fdr),
          return_diagnostics = TRUE
        )
        batch_replayed <- batch_detail$data
        stage_decisions$batch <- batch_detail$diagnostics
      } else {
        batch_replayed <- combat_batch_correction(drift_detail$data, batch)
        stage_decisions$batch <- data.frame(
          feature_id = if (is.null(rownames(data))) {
            paste0("feature_", seq_len(nrow(data)))
          } else {
            rownames(data)
          },
          method = "ComBat",
          scope = "global_all_features",
          actual_correction = TRUE,
          stringsAsFactors = FALSE
        )
      }

      normalization_detail <- normalize_by_dilution_factor(
        batch_replayed,
        processing = best_params$normalization,
        control_samples = control_samples,
        return_diagnostics = TRUE
      )
      stage_decisions$dilution <- normalization_detail$diagnostics
      replayed_final <- if (isTRUE(best_params$scale_by_batch)) {
        scale_by_batch(normalization_detail$data, batch)
      } else {
        normalization_detail$data
      }
      if (!isTRUE(all.equal(
        replayed_final,
        best_final_data,
        tolerance = 1e-8,
        check.attributes = TRUE
      ))) {
        stop(
          "Internal auto-mode replay did not reproduce the selected result.",
          call. = FALSE
        )
      }
      final_data <- replayed_final
    }
  } else {
    # Use fixed parameters
    if (is.null(batch)) {
      message("No batch information provided. Auto-detecting segments using fkPELT...")
      batch <- .auto_detect_batch(
        norm_data,
        pelt_penalty = pelt_penalty,
        run_order = run_order
      )
      message("Auto-detected ", max(batch), " segments as batch.")
      message("fkPELT penalty used: ", round(attr(batch, "pelt_penalty"), 4))
    } else {
      .validate_batch(batch, n_samples)
    }
    
    # Drift correction
    message(
      "Correcting drift using autocorrelation test (",
      test,
      ") with FDR threshold: ",
      fdr_threshold,
      "..."
    )
    if (identical(test, "Ljung-Box") && is.null(lag)) {
      message("Using adaptive Ljung-Box lag selection per batch segment.")
    }
    drift_result <- autocorrelation_correct(
      norm_data,
      run_order = run_order,
      batch = batch,
      lag = lag,
      test = test,
      ljung_box_fitdf = ljung_box_fitdf,
      detrend = detrend_non_autocorrelated,
      fdr_threshold = fdr_threshold,
      spline_method = spline_method,
      return_diagnostics = return_details
    )
    drift_corrected <- if (isTRUE(return_details)) drift_result$data else drift_result
    if (isTRUE(return_details)) {
      stage_decisions$drift <- drift_result$diagnostics
    }
    
    
    # Batch effect correction via ANOVA
    
    if (remove_batch_effects == "anova") {
      message(
        "Correcting batch effects using ANOVA-based residualization with FDR threshold: ",
        fdr_threshold,
        "..."
      )
      batch_result <- anova_batch_correction(
        drift_corrected,
        batch,
        fdr_threshold = fdr_threshold,
        return_diagnostics = return_details
      )
      batch_corrected <- if (isTRUE(return_details)) batch_result$data else batch_result
      if (isTRUE(return_details)) {
        stage_decisions$batch <- batch_result$diagnostics
      }
    } else {
      message("Testing and removing batch (batch) effects using ComBat")
      batch_corrected <- combat_batch_correction(drift_corrected, batch)
      if (isTRUE(return_details)) {
        stage_decisions$batch <- data.frame(
          feature_id = if (is.null(rownames(data))) {
            paste0("feature_", seq_len(nrow(data)))
          } else {
            rownames(data)
          },
          method = "ComBat",
          scope = "global_all_features",
          actual_correction = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Median adjustment using control samples if provided
    if (median_adjustment == "none") {
      message("Skipping median adjustment.")
      if (isTRUE(return_details)) {
        sample_ids <- colnames(data)
        if (is.null(sample_ids)) {
          sample_ids <- paste0("sample_", seq_len(ncol(data)))
        }
        stage_decisions$dilution <- data.frame(
          sample_id = sample_ids,
          altered = FALSE,
          processing = "none",
          stringsAsFactors = FALSE
        )
      }
    } else {
      message("Performing median adjustment using method: ",
              median_adjustment)
      dilution_result <- normalize_by_dilution_factor(
        batch_corrected,
        processing = median_adjustment,
        control_samples = control_samples,
        return_diagnostics = return_details
      )
      batch_corrected <- if (isTRUE(return_details)) {
        dilution_result$data
      } else {
        dilution_result
      }
      if (isTRUE(return_details)) {
        stage_decisions$dilution <- dilution_result$diagnostics
      }
    }
    
    # Optional scaling by batch
    if (scale_by_batch) {
      message("Scaling data by batch...")
      final_data <- scale_by_batch(batch_corrected, batch)
    } else {
      final_data <- batch_corrected
    }
    selected_parameters <- list(
      parameters = "fixed",
      batch_source = batch_source,
      pelt_penalty = attr(batch, "pelt_penalty"),
      test = test,
      ljung_box_fitdf = if (identical(test, "Ljung-Box")) {
        ljung_box_fitdf
      } else {
        "not used"
      },
      lag = if (is.null(lag) && identical(test, "Ljung-Box")) {
        "adaptive"
      } else if (identical(test, "Ljung-Box")) {
        as.integer(lag)
      } else {
        "not used"
      },
      fdr_threshold = fdr_threshold,
      median_adjustment = median_adjustment,
      detrend = detrend_non_autocorrelated,
      spline_method = spline_method,
      remove_batch_effects = remove_batch_effects,
      scale_by_batch = scale_by_batch
    )
  }
  
  message("Winn correction completed.")
  if (!isTRUE(return_details)) {
    return(final_data)
  }
  result <- list(
    data = final_data,
    selected_parameters = selected_parameters,
    batch = batch,
    batch_source = batch_source,
    run_order = if (is.null(run_order)) seq_len(n_samples) else run_order,
    control_samples = control_samples,
    stage_decisions = stage_decisions,
    input_policy = list(
      finite_nonnegative_required = TRUE,
      zero_is_observed = TRUE,
      automatic_zero_to_missing = FALSE
    )
  )
  class(result) <- "winn_result"
  result
}

###############################################################################
# Internal Helper Functions (not exported)
###############################################################################

.autocorrelation_model_df <- function(test, ljung_box_fitdf = 0L) {
  if (identical(test, "Ljung-Box")) {
    return(.validate_ljung_box_fitdf(ljung_box_fitdf))
  }
  0L
}

.resolve_autocorrelation_lag <- function(lag, n_obs, model_df = 0L) {
  max_allowed_lag <- n_obs - 1L
  min_required_lag <- model_df + 3L
  if (max_allowed_lag < min_required_lag) {
    return(NA_integer_)
  }
  if (is.null(lag)) {
    lag <- min(10L, floor(n_obs / 5))
  }
  lag <- max(as.integer(lag), min_required_lag)
  lag <- min(lag, max_allowed_lag)
  as.integer(lag)
}

.compute_batch_anova_pvalues <- function(z, batch) {
  batch <- factor(batch)
  n_groups <- nlevels(batch)
  n_samples <- ncol(z)
  if (n_groups < 2L || n_samples <= n_groups) {
    return(rep(1, nrow(z)))
  }
  if (anyNA(z)) {
    pvals <- rep(1, nrow(z))
    for (i in seq_len(nrow(z))) {
      keep <- !is.na(z[i, ])
      y <- z[i, keep]
      batch_i <- droplevels(batch[keep])
      group_count <- nlevels(batch_i)
      if (group_count < 2L || length(y) <= group_count) {
        next
      }
      design <- model.matrix(~ batch_i)
      fit <- stats::lm.fit(design, y)
      rss_full <- sum(fit$residuals^2)
      rss_null <- sum((y - mean(y))^2)
      df1 <- ncol(design) - 1L
      df2 <- length(y) - ncol(design)
      if (df1 <= 0L || df2 <= 0L) {
        next
      }
      if (rss_full <= .Machine$double.eps) {
        pvals[i] <- if (rss_null <= .Machine$double.eps) 1 else 0
        next
      }
      f_stat <- ((rss_null - rss_full) / df1) / (rss_full / df2)
      pvals[i] <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)
    }
    return(pvals)
  }
  group_indices <- split(seq_len(n_samples), batch)
  group_sizes <- lengths(group_indices)
  group_means <- vapply(
    group_indices,
    function(idx) rowMeans(z[, idx, drop = FALSE], na.rm = TRUE),
    numeric(nrow(z))
  )
  if (!is.matrix(group_means)) {
    group_means <- matrix(group_means, ncol = 1L)
  }
  overall_means <- rowMeans(z, na.rm = TRUE)
  centered_means <- sweep(group_means, 1, overall_means, "-")
  ss_between <- rowSums(sweep(centered_means^2, 2, group_sizes, "*"))
  ss_within <- numeric(nrow(z))
  for (j in seq_along(group_indices)) {
    idx <- group_indices[[j]]
    residuals <- z[, idx, drop = FALSE] - group_means[, j]
    ss_within <- ss_within + rowSums(residuals^2)
  }
  df1 <- n_groups - 1L
  df2 <- n_samples - n_groups
  f_stat <- (ss_between / df1) / (ss_within / df2)
  stable <- ss_within <= .Machine$double.eps
  f_stat[stable & ss_between <= .Machine$double.eps] <- 0
  f_stat[stable & ss_between > .Machine$double.eps] <- Inf
  stats::pf(f_stat, df1, df2, lower.tail = FALSE)
}

.prepare_anova_batch_correction <- function(data, batch) {
  z <- log1p(data)
  batch <- factor(batch)
  group_indices <- split(seq_len(ncol(z)), batch)
  group_means <- vapply(
    group_indices,
    function(idx) rowMeans(z[, idx, drop = FALSE], na.rm = TRUE),
    numeric(nrow(z))
  )
  if (!is.matrix(group_means)) {
    group_means <- matrix(group_means, ncol = 1L)
  }
  overall_means <- rowMeans(z, na.rm = TRUE)
  corrected_all <- z
  for (j in seq_along(group_indices)) {
    idx <- group_indices[[j]]
    corrected_all[, idx] <- z[, idx, drop = FALSE] - group_means[, j] + overall_means
  }
  list(
    original = z,
    corrected_all = corrected_all,
    padj = p.adjust(.compute_batch_anova_pvalues(z, batch), method = "fdr")
  )
}

.materialize_anova_batch_correction <- function(prepared, fdr_threshold) {
  corrected <- prepared$original
  sig <- which(prepared$padj < fdr_threshold)
  if (length(sig) > 0) {
    corrected[sig, ] <- prepared$corrected_all[sig, , drop = FALSE]
  }
  .inv_log1p(corrected, clamp_nonneg = TRUE)
}

.materialize_anova_candidates <- function(prepared, thresholds) {
  ascending_thresholds <- sort(unique(thresholds))
  current <- prepared$original
  threshold_groups <- list()
  matrices <- list()
  sorted_group_ids <- integer(length(ascending_thresholds))
  prev_threshold <- 0
  last_group_id <- 0L
  
  for (i in seq_along(ascending_thresholds)) {
    threshold_value <- ascending_thresholds[i]
    new_sig <- which(prepared$padj >= prev_threshold & prepared$padj < threshold_value)
    if (!length(new_sig) && last_group_id > 0L) {
      sorted_group_ids[i] <- last_group_id
      threshold_groups[[last_group_id]] <- c(threshold_groups[[last_group_id]], threshold_value)
    } else {
      if (length(new_sig)) {
        current[new_sig, ] <- prepared$corrected_all[new_sig, , drop = FALSE]
      }
      last_group_id <- last_group_id + 1L
      sorted_group_ids[i] <- last_group_id
      threshold_groups[[last_group_id]] <- threshold_value
      matrices[[last_group_id]] <- .inv_log1p(current, clamp_nonneg = TRUE)
    }
    prev_threshold <- threshold_value
  }
  
  candidates <- list()
  seen_groups <- integer(0)
  for (threshold_value in thresholds) {
    group_id <- sorted_group_ids[match(threshold_value, ascending_thresholds)]
    if (group_id %in% seen_groups) {
      next
    }
    seen_groups <- c(seen_groups, group_id)
    candidates[[length(candidates) + 1L]] <- list(
      threshold = threshold_value,
      thresholds = threshold_groups[[group_id]],
      matrix = matrices[[group_id]]
    )
  }
  candidates
}

.is_valid_pelt_penalty_value <- function(pelt_penalty) {
  if (is.numeric(pelt_penalty) && length(pelt_penalty) == 1L &&
    !is.na(pelt_penalty) && pelt_penalty > 0) {
    return(TRUE)
  }
  if (is.character(pelt_penalty) && length(pelt_penalty) == 1L) {
    return(tolower(pelt_penalty) %in% c("bic", "mbic"))
  }
  FALSE
}

.make_batch_from_change_points <- function(change_points, n) {
  tau <- c(0, change_points, n)
  batch <- rep(NA_integer_, n)
  for (i in seq_len(length(tau) - 1L)) {
    batch[(tau[i] + 1L):tau[i + 1L]] <- i
  }
  batch
}

.make_fk_knots <- function(n) {
  num_knots <- floor(n / 60) + 2
  if (num_knots <= 2L) {
    return(numeric(0))
  }
  knots <- seq(1, n, length.out = num_knots)
  knots <- knots[-c(1, length(knots))]
  ifelse(floor(knots) == knots, knots + 0.5, knots)
}

.estimate_segmentation_variance <- function(signal) {
  diffs <- diff(signal)
  sigma <- if (length(diffs) > 1L) {
    mad(diffs, center = 0, constant = 1.4826, na.rm = TRUE)
  } else {
    sd(signal, na.rm = TRUE)
  }
  sigma <- max(sigma, .Machine$double.eps)
  sigma^2
}

.resolve_pelt_penalty <- function(agg_signal, pelt_penalty = NULL) {
  n <- length(agg_signal)
  sigma2 <- .estimate_segmentation_variance(agg_signal)
  bic_penalty <- sigma2 * log(max(n, 2L))
  mbic_penalty <- 3 * log(max(n, 2L))
  if (is.numeric(pelt_penalty)) {
    return(as.numeric(pelt_penalty))
  }
  penalty_key <- if (is.null(pelt_penalty)) {
    "mbic"
  } else {
    tolower(pelt_penalty)
  }
  if (identical(penalty_key, "bic")) {
    return(bic_penalty)
  }
  mbic_penalty
}

.auto_detect_batch <- function(data, pelt_penalty = NULL, run_order = NULL) {
  # Auto-detect segments as batch using fkPELT based on aggregated median signal
  .validate_run_order(run_order, batch = NULL, n_samples = ncol(data))
  acquisition_order <- if (is.null(run_order)) {
    seq_len(ncol(data))
  } else {
    order(run_order)
  }
  agg_signal <- apply(data[, acquisition_order, drop = FALSE], 2, median)
  n <- length(agg_signal)
  knots <- .make_fk_knots(n)
  resolved_penalty <- .resolve_pelt_penalty(
    agg_signal = agg_signal,
    pelt_penalty = pelt_penalty
  )
  change_points <- .fkPELT(agg_signal, knots, penalty = resolved_penalty)
  detected_in_order <- .make_batch_from_change_points(change_points, n)
  batch <- integer(n)
  batch[acquisition_order] <- detected_in_order
  attr(batch, "pelt_penalty") <- resolved_penalty
  return(batch)
}

.mean_control_correlation <- function(data, control_samples) {
  # Calculate mean pairwise correlation among control sample columns
  control_data <- data[, control_samples, drop = FALSE]
  if (ncol(control_data) < 2)
    return(NA)
  corr_matrix <- cor(control_data, use = "pairwise.complete.obs")
  lower_tri <- corr_matrix[lower.tri(corr_matrix)]
  return(mean(lower_tri, na.rm = TRUE))
}

.calculate_quality_score <- function(data, control_samples) {
  # Calculate comprehensive quality score for parameter optimization
  tryCatch({
    control_data <- data[, control_samples, drop = FALSE]
    
    # Calculate coefficient of variation (CV) for each metabolite in controls
    cvs <- apply(control_data, 1, function(x) {
      mu <- mean(x, na.rm = TRUE)
      if (mu == 0)
        return(NA)
      sd(x, na.rm = TRUE) / abs(mu)
    })
    mean_cv <- mean(cvs, na.rm = TRUE)
    
    # Calculate mean pairwise correlation within controls
    mean_correlation <- .mean_control_correlation(data, control_samples)
    
    # Combined score: higher correlation is better, lower CV is better
    # Use weighted combination with correlation being primary metric
    if (is.na(mean_correlation) || is.na(mean_cv)) {
      score <- NA
    } else {
      # Normalize and combine: correlation (0 to 1) - penalty for high CV
      score <- mean_correlation - (mean_cv * 0.5)  # CV penalty weighted at 50%
    }
    
    return(list(
      score = score,
      mean_cv = mean_cv,
      mean_correlation = mean_correlation
    ))
  }, error = function(e) {
    return(list(
      score = NA,
      mean_cv = NA,
      mean_correlation = NA
    ))
  })
}

.tag_spline_fit <- function(fit, smoother_type, fallback_path, basis_dimension = NA_real_, warnings = character()) {
  if (is.null(fit)) {
    return(NULL)
  }
  attr(fit, "winn_smoother_type") <- smoother_type
  attr(fit, "winn_fallback_path") <- fallback_path
  attr(fit, "winn_basis_dimension") <- as.numeric(basis_dimension)
  attr(fit, "winn_warnings") <- as.character(warnings)
  fit
}

.fit_conservative_spline <- function(y, x) {
  # Conservative spline fitting with optimal parameters for drift removal
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for spline detrending. Please install it.")
  }
  n <- length(y)
  
  # For very short segments, use linear regression
  if (n < 5) {
    fit <- tryCatch(
      lm(y ~ x),
      error = function(e)
        NULL
    )
    return(.tag_spline_fit(
      fit,
      smoother_type = "linear",
      fallback_path = "short_segment_linear",
      basis_dimension = 2
    ))
  }
  
  # For longer segments, use conservative GAM parameters
  # Calculate optimal k (number of basis functions) - conservative approach
  k_max <- max(4, min(10, floor(n / 5)))
  
  # Try multiple approaches in order of conservatism
  # 1. P-splines with REML and gamma penalty
  fit <- tryCatch({
    mgcv::gam(y ~ mgcv::s(x, bs = "ps", k = k_max),
              method = "REML",
              gamma = 1.4)  # Conservative gamma > 1
  }, error = function(e)
    NULL)
  
  if (!is.null(fit))
    return(.tag_spline_fit(
      fit,
      smoother_type = "conservative_ps_reml",
      fallback_path = "ps_reml_gamma_1.4",
      basis_dimension = k_max
    ))
  
  # 2. Thin plate splines with conservative settings
  fit <- tryCatch({
    mgcv::gam(y ~ mgcv::s(x, bs = "tp", k = k_max),
              method = "REML",
              gamma = 1.2)
  }, error = function(e)
    NULL)
  
  if (!is.null(fit))
    return(.tag_spline_fit(
      fit,
      smoother_type = "conservative_tp_reml",
      fallback_path = "ps_failed->tp_reml_gamma_1.2",
      basis_dimension = k_max
    ))
  
  # 3. Cubic regression splines with fixed df (very conservative)
  fit <- tryCatch({
    mgcv::gam(y ~ mgcv::s(
      x,
      bs = "cr",
      k = min(6, k_max),
      fx = TRUE
    ))
  }, error = function(e)
    NULL)
  
  if (!is.null(fit))
    return(.tag_spline_fit(
      fit,
      smoother_type = "conservative_cr_fixed",
      fallback_path = "ps_failed->tp_failed->cr_fixed",
      basis_dimension = min(6, k_max)
    ))
  
  # 4. Fallback to LOESS
  fit <- tryCatch({
    loess_fit <- loess(y ~ x, span = 0.3, degree = 1)
    list(fitted.values = loess_fit$fitted)
  }, error = function(e)
    NULL)
  
  if (!is.null(fit))
    return(.tag_spline_fit(
      fit,
      smoother_type = "loess",
      fallback_path = "ps_failed->tp_failed->cr_failed->loess",
      basis_dimension = NA_real_
    ))
  
  # 5. Final fallback to linear regression
  fit <- tryCatch(
    lm(y ~ x),
    error = function(e)
      NULL
  )
  .tag_spline_fit(
    fit,
    smoother_type = "linear",
    fallback_path = "ps_failed->tp_failed->cr_failed->loess_failed->linear",
    basis_dimension = 2
  )
}

.fkPELT <- function(data, knots, penalty = NULL) {
  if (is.null(data))
    stop("Data cannot be NULL in .fkPELT")
  n <- length(data)
  if (!is.numeric(penalty) || length(penalty) != 1 || penalty <= 0) {
    stop("penalty must be a single positive numeric value.")
  }
  pruning_constant <- penalty
  f <- numeric(n + 1)
  f[1] <- -penalty
  cp <- vector("list", n + 1)
  cp[[1]] <- numeric(0)
  R <- vector("list", n)
  R[[1]] <- 0
  
  for (t in 1:n) {
    m <- length(R[[t]])
    neglog <- numeric(m)
    for (r in seq_len(m)) {
      start_idx <- R[[t]][r] + 1
      seg_data <- data[start_idx:t]
      neglog[r] <- .fksplinecost(seg_data,
                                 knots,
                                 index1 = R[[t]][r] + 1,
                                 index2 = t)
    }
    stat <- numeric(m)
    for (r in seq_len(m)) {
      stat[r] <- f[R[[t]][r] + 1] + penalty + neglog[r]
    }
    f[t + 1] <- min(stat)
    t1 <- R[[t]][which.min(stat)]
    cp[[t + 1]] <- c(cp[[t1 + 1]], t1)
    R[[t + 1]] <- numeric(0)
    for (r in seq_len(m)) {
      tot <- f[R[[t]][r] + 1] + pruning_constant + neglog[r]
      if (tot <= f[t + 1] && t < n) {
        R[[t + 1]] <- c(R[[t + 1]], R[[t]][r])
      }
    }
    if (t < n && t > 79) {
      R[[t + 1]] <- c(R[[t + 1]], t - 39)
    }
  }
  cp_final <- cp[[n + 1]]
  cp_final <- cp_final[cp_final != 0]
  return(cp_final)
}

.fksplinecost <- function(data,
                          knots,
                          index1 = 1,
                          index2 = length(data)) {
  size <- length(data)
  if (size == 1)
    return(0)
  if (size < 5) {
    sd_val <- sd(data, na.rm = TRUE)
    mu_val <- mean(data, na.rm = TRUE)
    if (sd_val <= 1e-5)
      return(0)
    neglog <- 2 * (sum((data - mu_val)^2 / (2 * sd_val^2)) + size * log(sd_val * sqrt(2 * pi)))
    return(neglog)
  }
  cov <- index1:index2
  newknots <- knots[knots < index2 & knots > index1]
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required for .fksplinecost. Please install it.")
  }
  fit <- tryCatch(
    lm(data ~ splines::ns(
      cov, knots = newknots, intercept = TRUE
    )),
    error = function(e)
      NULL
  )
  if (is.null(fit))
    return(Inf)
  mu_val <- fit$fitted.values
  sd_val <- sd(data - mu_val, na.rm = TRUE)
  if (sd_val <= 1e-5)
    return(0)
  neglog <- 2 * (sum((data - mu_val)^2 / (2 * sd_val^2)) + size * log(sd_val * sqrt(2 * pi)))
  return(neglog)
}
