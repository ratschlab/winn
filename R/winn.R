###############################################################################
# Exported Functions
###############################################################################

#' Normalize Data Matrix by Dilution Factor
#'
#' This function normalizes a data matrix by dilution factors or alternatively shrinks sample measurements
#' if any samples have dilution factors that are more than one standard deviation away from the mean dilution factor.
#'
#' @param data A numeric matrix or data frame where rows represent metabolites and columns represent samples.
#' @param processing A character string specifying the processing method to use. Options are "shrink" (default) or "normalize".
#' @param control_samples An optional numeric vector specifying which columns correspond to control samples.
#' If provided, the reference spectrum is calculated using only these samples.
#' @return A numeric matrix of normalized intensities.
#' @examples
#' # Example usage:
#' your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
#' normalized_data <- normalize_by_dilution_factor(your_data_matrix, control_samples = 1:4)
#' @export
normalize_by_dilution_factor <- function(data,
                                         processing = "shrink",
                                         control_samples = NULL) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame.")
  }
  if (is.data.frame(data))
    data <- as.matrix(data)
  
  if (!is.null(control_samples)) {
    # Calculate reference spectrum using control samples only
    reference_spectrum <- apply(data[, control_samples, drop = FALSE], 1, median, na.rm = TRUE)
  } else {
    reference_spectrum <- apply(data, 1, median, na.rm = TRUE)
  }
  
  quotients <- data / reference_spectrum
  dilution_factors <- apply(quotients, 2, median, na.rm = TRUE)
  
  mean_dilution <- mean(dilution_factors, na.rm = TRUE)
  stdev_dilution <- sd(dilution_factors, na.rm = TRUE)
  
  normalized_data <- data
  if (processing == "shrink") {
    for (i in seq_along(dilution_factors)) {
      if (dilution_factors[i] < (mean_dilution - stdev_dilution)) {
        normalized_data[, i] <- (normalized_data[, i] / dilution_factors[i]) *
          (mean_dilution - stdev_dilution)
      }
      if (dilution_factors[i] > (mean_dilution + stdev_dilution)) {
        normalized_data[, i] <- (normalized_data[, i] / dilution_factors[i]) *
          (mean_dilution + stdev_dilution)
      }
    }
  } else if (processing == "normalize") {
    for (i in seq_along(dilution_factors)) {
      normalized_data[, i] <- normalized_data[, i] / dilution_factors[i]
    }
  }
  return(normalized_data)
}

#' Adjust Outliers in Data Matrix Using MAD
#'
#' This function adjusts outliers in a data matrix on a per-metabolite basis using the Median Absolute Deviation (MAD).
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
  
  adjusted_data <- data
  for (i in 1:nrow(data)) {
    med <- median(data[i, ], na.rm = TRUE)
    mad_val <- mad(data[i, ], na.rm = TRUE)
    lower_threshold <- med - 4 * mad_val
    upper_threshold <- med + 4 * mad_val
    is_outlier <- data[i, ] <= lower_threshold |
      data[i, ] >= upper_threshold
    if (any(is_outlier)) {
      # For upper outliers
      if (any(is_outlier & data[i, ] >= upper_threshold)) {
        ref_val <- max(data[i, !is_outlier &
                              data[i, ] < upper_threshold], na.rm = TRUE)
        interp_val <- med + 3 * mad_val
        adjusted_data[i, is_outlier &
                        data[i, ] >= upper_threshold] <-
          approx(c(ref_val, max(data[i, is_outlier &
                                       data[i, ] >= upper_threshold])),
                 c(interp_val, upper_threshold),
                 xout = data[i, is_outlier &
                               data[i, ] >= upper_threshold])$y
      } else {
        ref_val <- min(data[i, is_outlier & data[i, ] <= lower_threshold])
        interp_val <- med - 3 * mad_val
        adjusted_data[i, is_outlier &
                        data[i, ] <= lower_threshold] <-
          approx(c(ref_val, min(data[i, !is_outlier &
                                       data[i, ] > lower_threshold], na.rm = TRUE)),
                 c(lower_threshold, interp_val),
                 xout = data[i, is_outlier &
                               data[i, ] <= lower_threshold])$y
      }
    }
  }
  return(adjusted_data)
}

# ---- internal: correct inverse for log1p ----
.inv_log1p <- function(z, clamp_nonneg = TRUE) {
  x <- expm1(z)
  if (clamp_nonneg) {
    x <- pmax(x, 0)
  }
  x
}

#' Correct for Drift in Data Using Autocorrelation Correction
#'
#' This function corrects for drift effects in metabolomics data by detrending based on run order within each batch segment.
#'
#' @param data A numeric matrix with rows representing metabolites and columns representing samples.
#' @param run_order An optional numeric vector representing the run order of the samples.
#' @param batch A numeric vector indicating the batch (or segment) assignment for each sample.
#' @param lag An optional integer specifying the lag to be used in the autocorrelation test.
#' If `NULL`, the lag is selected adaptively for each batch segment using
#' `max(min(10, floor(n / 5)), df + 3)` and then capped at `n - 1`.
#' @param test A character string specifying the autocorrelation test to use ("Ljung-Box" or "DW").
#' @param detrend A character string indicating the method for detrending ("mean" or "spline").
#' @param fdr_threshold A numeric value specifying the FDR threshold for significance.
#' @param spline_method A character string specifying the spline method when detrend="spline" ("conservative" or "standard").
#' @return A numeric matrix with drift corrected.
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
                                    detrend = "mean",
                                    fdr_threshold = 0.05,
                                    spline_method = "conservative") {
  if (!is.matrix(data))
    stop("Data must be a numeric matrix.")
  if (!is.null(run_order) && length(run_order) != ncol(data)) {
    stop("Length of run_order must match number of columns in data.")
  }
  if (is.null(batch) || length(batch) != ncol(data)) {
    stop("batch vector must be provided and its length must equal number of columns in data.")
  }
  if (!is.null(lag) && (length(lag) != 1 || !is.numeric(lag) || is.na(lag) ||
    lag <= 0 || lag != round(lag))) {
    stop("lag must be NULL or a positive integer.")
  }
  
  detrended_data <- data
  batch_vec <- batch
  unique_batch <- unique(batch_vec)
  
  for (b in unique_batch) {
    idx <- which(batch_vec == b)
    segment <- data[, idx, drop = FALSE]
    segment_n <- length(idx)
    model_df <- .autocorrelation_model_df(test = test)
    segment_lag <- .resolve_autocorrelation_lag(
      lag = lag,
      n_obs = segment_n,
      model_df = model_df
    )
    seg_run <- if (!is.null(run_order)) run_order[idx] else seq_along(idx)
    if (detrend == "spline") {
      spline_vals <- apply(segment, 1, function(x) {
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
        if (is.null(fit))
          lm(log_seg ~ seg_run)$fitted.values - mean(lm(log_seg ~ seg_run)$fitted.values)
        else
          (fit$fitted.values - mean(fit$fitted.values))
      })
      z <- log1p(data[, idx, drop = FALSE]) - t(spline_vals)
      detrended_data[, idx] <- .inv_log1p(z, clamp_nonneg = TRUE)
    } else if (detrend == "mean") {
      # Use log1p-transformed data for autocorrelation testing (consistent with detrending)
      p_vals <- apply(segment, 1, function(x) {
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
            error = function(e)
              NA
          )
        } else if (test == "DW") {
          if (!requireNamespace("lmtest", quietly = TRUE)) {
            stop("Package 'lmtest' is required for DW test. Please install it.")
          }
          tryCatch(
            lmtest::dwtest(x_log ~ seg_run)$p.value,
            error = function(e)
              NA
          )
        } else {
          stop("Invalid test method. Use 'Ljung-Box' or 'DW'.")
        }
      })
      p_vals <- p.adjust(p_vals, method = "fdr")
      correct_ids <- which(p_vals < fdr_threshold)
      if (length(correct_ids) > 0) {
        spline_vals <- apply(segment[correct_ids, , drop = FALSE], 1, function(x) {
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
          if (is.null(fit))
            lm(log_seg ~ seg_run)$fitted.values - mean(lm(log_seg ~ seg_run)$fitted.values)
          else
            (fit$fitted.values - mean(fit$fitted.values))
        })
        z <- log1p(data[correct_ids, idx, drop = FALSE]) - t(spline_vals)
        detrended_data[correct_ids, idx] <- .inv_log1p(z, clamp_nonneg = TRUE)
      }
      non_correct_ids <- setdiff(seq_len(nrow(data)), correct_ids)
      if (length(non_correct_ids) > 0) {
        detrended_data[non_correct_ids, idx] <- data[non_correct_ids, idx]
      }
    }
  }
  return(detrended_data)
}

#' Perform ANOVA-based Mean-Only Batch Correction
#'
#' This function runs an ANOVA test on each metabolite to detect batch effects,
#' and then corrects significant batch effects by subtracting the estimated
#' batch-specific shifts (while preserving the overall mean).
#'
#' @param data A numeric matrix (metabolites × samples).
#' @param batch A factor or numeric vector indicating batch for each sample.
#' @param fdr_threshold Significance threshold for FDR-adjusted p-values.
#' @return A numeric matrix of corrected intensities.
#' @examples
#' mat <- matrix(rnorm(200, mean = 100, sd = 15), nrow=20)
#' batch <- rep(1:4, length.out=ncol(mat))
#' corrected <- anova_batch_correction(mat, batch, fdr_threshold=0.05)
#' @export
anova_batch_correction <- function(data, batch, fdr_threshold = 0.05) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("Data must be a numeric matrix.")
  }
  if (length(batch) != ncol(data)) {
    stop("Length of batch must match number of columns in data.")
  }
  batch <- factor(batch)
  if (length(unique(batch)) < 2) {
    message("Only one batch detected. Skipping ANOVA-based correction.")
    return(data)
  }
  z <- log1p(data)
  n_met <- nrow(z)
  pvals <- numeric(n_met)
  for (i in seq_len(n_met)) {
    df <- data.frame(value = z[i, ], batch = batch)
    a <- aov(value ~ batch, data = df)
    pvals[i] <- summary(a)[[1]]$`Pr(>F)`[1]
  }
  padj <- p.adjust(pvals, method = "fdr")
  
  corrected <- z
  sig <- which(padj < fdr_threshold)
  if (length(sig) > 0) {
    overall_means <- rowMeans(z, na.rm = TRUE)
    for (i in sig) {
      # compute batch-specific and overall means
      batch_means <- tapply(z[i, ], batch, mean, na.rm = TRUE)
      shift <- batch_means - overall_means[i]
      # subtract the shift for each sample
      corrected[i, ] <- z[i, ] - shift[as.character(batch)]
    }
  }
  return(.inv_log1p(corrected, clamp_nonneg = TRUE))
}

#' Perform ComBat Batch Correction by batch
#'
#' This function applies the empirical Bayes ComBat method to correct batch effects
#' by batch, adjusting both location and scale parameters across batch.
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
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("Data must be a numeric matrix.")
  }
  if (length(batch) != ncol(data)) {
    stop("Length of batch must match number of columns in data.")
  }
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for ComBat correction. Please install it.")
  }
  batch <- factor(batch)
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
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("Data must be a numeric matrix.")
  }
  if (length(batch) != ncol(data)) {
    stop("Length of batch must match number of columns in data.")
  }
  
  scaled_data <- data
  batch_vec <- batch
  unique_batch <- unique(batch_vec)
  for (b in unique_batch) {
    idx <- which(batch_vec == b)
    row_means <- rowMeans(data[, idx, drop = FALSE], na.rm = TRUE)
    row_sds <- apply(data[, idx, drop = FALSE], 1, sd, na.rm = TRUE)
    # Avoid division by zero by setting zero standard deviations to 1
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
#' correlation while minimizing coefficient of variation. This includes optimization of spline methods,
#' FDR thresholds, normalization approaches, and scaling options.
#'
#' @param data A numeric matrix or data frame where rows represent metabolites and columns represent samples.
#' @param batch An optional numeric vector indicating batch assignments for each sample. If NULL, segments will be auto-detected.
#' @param run_order An optional numeric vector representing the run order of samples.
#' @param control_samples An optional numeric vector representing the columns corresponding to control samples. If provided,
#' these will be used for normalization and parameter tuning.
#' @param parameters An optional character string specifying whether to use fixed ("fixed") or auto-detected ("auto") parameters in presence of control samples (default: "auto").
#' @param fdr_threshold A numeric value specifying the FDR threshold for drift and batch corrections (default: 0.05).
#' @param median_adjustment A character string specifying the method for median adjustment ("shrink", "normalize", or "none").
#' @param detrend_non_autocorrelated A character string specifying the method for detrending non-autocorrelated metabolites ("mean" or "spline").
#' @param spline_method A character string specifying the spline method when detrend="spline" ("conservative" for robust drift removal or "standard" for traditional approach).
#' @param remove_batch_effects A character string specifying the method for removing batch effects ("anova" or "combat").
#' @param test A character string specifying the autocorrelation test ("Ljung-Box" or "DW").
#' @param lag An optional integer specifying the lag for the autocorrelation test.
#' If `NULL`, `autocorrelation_correct()` selects the lag adaptively for each
#' batch segment.
#' @param scale_by_batch Logical indicating whether to scale data by batch after corrections.
#' @param pelt_penalty Optional fkPELT penalty specification. Use a positive
#' numeric value, `"bic"`, `"mbic"`, or `NULL`. When `NULL`, fkPELT uses
#' the conservative MBIC default.
#' @return A numeric matrix of corrected intensities.
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
                 scale_by_batch = FALSE,
                 pelt_penalty = NULL) {
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }
  if (is.data.frame(data))
    data <- as.matrix(data)
  if (!is.numeric(data))
    stop("Data must be numeric.")
  if (any(is.na(data) | is.infinite(data))) {
    warning("Data contains NA or infinite values. Results may be unreliable.")
  }
  
  n_samples <- ncol(data)
  n_metabolites <- nrow(data)
  
  if (n_samples < 3)
    stop("At least 3 samples required.")
  if (n_metabolites < 1)
    stop("At least 1 metabolite required.")
  
  if (!is.null(run_order) && length(run_order) != n_samples) {
    stop("Length of run_order must match number of columns in data.")
  }
  
  if (!is.null(control_samples)) {
    if (any(is.na(match(control_samples, seq_len(n_samples))))) {
      stop("Control samples should refer to valid column numbers in the data matrix.")
    }
    if (length(control_samples) < 2 && parameters == "auto") {
      stop("At least 2 control samples required for parameter optimization.")
    }
  }
  
  # Parameter validation
  if (!parameters %in% c("fixed", "auto")) {
    stop("parameters must be either 'fixed' or 'auto'.")
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
  if (!test %in% c("Ljung-Box", "DW")) {
    stop("test must be either 'Ljung-Box' or 'DW'.")
  }
  
  if (fdr_threshold <= 0 || fdr_threshold >= 1) {
    stop("fdr_threshold must be between 0 and 1.")
  }
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
  
  
  # Outlier adjustment
  message("Adjusting outliers using MAD...")
  norm_data <- adjust_outliers_mad(data)
  
  # If control samples are provided and auto parameter detection is enabled, perform grid search
  if (!is.null(control_samples) && parameters == "auto") {
    message("Auto-detecting optimal parameters using control samples...")
    
    # Define grid for parameter search
    batch_options <- if (!is.null(batch))
      c("provided")
    else
      c("auto")
    tests <- c("Ljung-Box", "DW")
    normalizations <- c("shrink", "normalize")
    acorr_fdr_options <- c(0.1, 0.05, 0.01)
    anova_fdr_options <- c(0.1, 0.05, 0.01)
    scale_options <- if (scale_by_batch) {
      c(TRUE)
    } else {
      c(FALSE)
    }
    spline_methods <- c("conservative", "standard")  # Add spline method optimization
    
    best_score <- -Inf
    best_final_data <- NULL
    best_params <- list()
    
    # Progress tracking for parameter optimization
    total_combinations <- length(batch_options) * length(tests) * length(spline_methods) *
      length(acorr_fdr_options) * length(anova_fdr_options) *
      length(normalizations) * length(scale_options)
    current_combination <- 0
    
    for (batch_option in batch_options) {
      current_batch <- if (batch_option == "provided") {
        batch
      } else {
        .auto_detect_batch(
          norm_data,
          pelt_penalty = pelt_penalty
        )
      }
      
      for (current_test in tests) {
        for (spline_method in spline_methods) {
          for (acorr_fdr in acorr_fdr_options) {
            # Drift correction with current autocorrelation parameters
            drift_corrected <- tryCatch({
              autocorrelation_correct(
                norm_data,
                run_order = run_order,
                batch = current_batch,
                lag = lag,
                test = current_test,
                detrend = detrend_non_autocorrelated,
                fdr_threshold = acorr_fdr,
                spline_method = spline_method
              )
            }, error = function(e) {
              # Silently skip failed parameter combinations during optimization
              return(NULL)
            })
            
            if (is.null(drift_corrected))
              next
            
            for (anova_fdr in anova_fdr_options) {
              # Batch effect correction with current ANOVA FDR
              batch_corrected <- tryCatch({
                if (remove_batch_effects == "anova") {
                  anova_batch_correction(drift_corrected,
                                         current_batch,
                                         fdr_threshold = anova_fdr)
                } else {
                  combat_batch_correction(drift_corrected, current_batch)
                }
              }, error = function(e) {
                return(NULL)
              })
              
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
                  return(NULL)
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
                  
                  # Calculate quality metrics
                  quality_metrics <- .calculate_quality_score(final_data, control_samples)
                  if (is.na(quality_metrics$score))
                    next
                  
                  if (quality_metrics$score > best_score) {
                    best_score <- quality_metrics$score
                    best_final_data <- final_data
                    best_params <- list(
                      batch_option = batch_option,
                      pelt_penalty = attr(current_batch, "pelt_penalty"),
                      test = current_test,
                      lag = if (is.null(lag) && identical(current_test, "Ljung-Box")) {
                        "adaptive"
                      } else if (identical(current_test, "Ljung-Box")) {
                        as.character(as.integer(lag))
                      } else {
                        "not used"
                      },
                      spline_method = spline_method,
                      acorr_fdr = acorr_fdr,
                      anova_fdr = anova_fdr,
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
  } else {
    # Use fixed parameters
    if (is.null(batch)) {
      message("No batch information provided. Auto-detecting segments using fkPELT...")
      batch <- .auto_detect_batch(
        norm_data,
        pelt_penalty = pelt_penalty
      )
      message("Auto-detected ", max(batch), " segments as batch.")
      message("fkPELT penalty used: ", round(attr(batch, "pelt_penalty"), 4))
    } else {
      if (length(batch) != n_samples) {
        stop("Length of batch must match number of columns in data.")
      }
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
    drift_corrected <- autocorrelation_correct(
      norm_data,
      run_order = run_order,
      batch = batch,
      lag = lag,
      test = test,
      detrend = detrend_non_autocorrelated,
      fdr_threshold = fdr_threshold,
      spline_method = spline_method
    )
    
    
    # Batch effect correction via ANOVA
    
    if (remove_batch_effects == "anova") {
      message(
        "Correcting batch effects using ANOVA-based residualization with FDR threshold: ",
        fdr_threshold,
        "..."
      )
      batch_corrected <- anova_batch_correction(drift_corrected, batch, fdr_threshold = fdr_threshold)
    } else {
      message("Testing and removing batch (batch) effects using ComBat")
      batch_corrected <- combat_batch_correction(drift_corrected, batch)
      
    }
    
    # Median adjustment using control samples if provided
    if (median_adjustment == "none") {
      message("Skipping median adjustment.")
      batch_corrected <- batch_corrected
    } else {
      message("Performing median adjustment using method: ",
              median_adjustment)
      batch_corrected <- normalize_by_dilution_factor(batch_corrected,
                                                      processing = median_adjustment,
                                                      control_samples = control_samples)
    }
    
    # Optional scaling by batch
    if (scale_by_batch) {
      message("Scaling data by batch...")
      final_data <- scale_by_batch(batch_corrected, batch)
    } else {
      final_data <- batch_corrected
    }
  }
  
  message("Winn correction completed.")
  return(final_data)
}

###############################################################################
# Internal Helper Functions (not exported)
###############################################################################

.autocorrelation_model_df <- function(test) {
  if (identical(test, "Ljung-Box")) {
    return(1L)
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

.auto_detect_batch <- function(data, pelt_penalty = NULL) {
  # Auto-detect segments as batch using fkPELT based on aggregated median signal
  agg_signal <- apply(data, 2, median, na.rm = TRUE)
  n <- length(agg_signal)
  knots <- .make_fk_knots(n)
  resolved_penalty <- .resolve_pelt_penalty(
    agg_signal = agg_signal,
    pelt_penalty = pelt_penalty
  )
  change_points <- .fkPELT(agg_signal, knots, penalty = resolved_penalty)
  batch <- .make_batch_from_change_points(change_points, n)
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
    return(fit)
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
    return(fit)
  
  # 2. Thin plate splines with conservative settings
  fit <- tryCatch({
    mgcv::gam(y ~ mgcv::s(x, bs = "tp", k = k_max),
              method = "REML",
              gamma = 1.2)
  }, error = function(e)
    NULL)
  
  if (!is.null(fit))
    return(fit)
  
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
    return(fit)
  
  # 4. Fallback to LOESS
  fit <- tryCatch({
    loess_fit <- loess(y ~ x, span = 0.3, degree = 1)
    list(fitted.values = loess_fit$fitted)
  }, error = function(e)
    NULL)
  
  if (!is.null(fit))
    return(fit)
  
  # 5. Final fallback to linear regression
  return(tryCatch(
    lm(y ~ x),
    error = function(e)
      NULL
  ))
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

###############################################################################
# Unit Tests (for interactive sessions)
###############################################################################
if (interactive()) {
  library(testthat)
  
  test_that("normalize_by_dilution_factor works (with and without control samples)",
            {
              set.seed(1)
              mat <- matrix(rnorm(100, mean = 100, sd = 15), nrow = 10)
              res1 <- normalize_by_dilution_factor(mat)
              res2 <- normalize_by_dilution_factor(mat, control_samples = 1:5)
              expect_equal(dim(res1), dim(mat))
              expect_equal(dim(res2), dim(mat))
            })
  
  test_that("winn works with auto-detected batch", {
    set.seed(1)
    mat <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
    res <- winn(mat)
    expect_equal(dim(res), dim(mat))
  })
  
  test_that("winn works with provided batch", {
    set.seed(1)
    mat <- matrix(rnorm(400, mean = 100, sd = 15), nrow = 20)
    batch <- rep(1:4, each = 5)
    run_order <- seq_len(ncol(mat))
    res <- winn(mat, batch = batch, run_order = run_order)
    expect_equal(dim(res), dim(mat))
  })
  
  test_that("winn auto detection with control samples works", {
    set.seed(1)
    mat <- matrix(rnorm(400, mean = 100, sd = 15), nrow = 20)
    batch <- rep(1:4, each = 5)
    run_order <- seq_len(ncol(mat))
    control_samples <- 1:4
    res <- winn(
      mat,
      batch = batch,
      run_order = run_order,
      control_samples = control_samples,
      parameters = "auto"
    )
    expect_equal(dim(res), dim(mat))
    # Optionally, one might test that the mean correlation among controls is above a minimal threshold.
  })
}
