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
#' your_data_matrix <- matrix(rnorm(200), nrow = 20)
#' normalized_data <- normalize_by_dilution_factor(your_data_matrix, control_samples = 1:4)
#' @export
normalize_by_dilution_factor <- function(data, processing = "shrink", control_samples = NULL) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame.")
  }
  if (is.data.frame(data)) data <- as.matrix(data)
  
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
#' your_data_matrix <- matrix(rnorm(200), nrow = 20)
#' adjusted_data <- adjust_outliers_mad(your_data_matrix)
#' @export
adjust_outliers_mad <- function(data) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame.")
  }
  if (is.data.frame(data)) data <- as.matrix(data)
  
  adjusted_data <- data
  for (i in 1:nrow(data)) {
    med <- median(data[i, ], na.rm = TRUE)
    mad_val <- mad(data[i, ], na.rm = TRUE)
    lower_threshold <- med - 4 * mad_val
    upper_threshold <- med + 4 * mad_val
    is_outlier <- data[i, ] <= lower_threshold | data[i, ] >= upper_threshold
    if (any(is_outlier)) {
      # For upper outliers
      if (any(is_outlier & data[i, ] >= upper_threshold)) {
        ref_val <- max(data[i, !is_outlier & data[i, ] < upper_threshold], na.rm = TRUE)
        interp_val <- med + 3 * mad_val
        adjusted_data[i, is_outlier & data[i, ] >= upper_threshold] <-
          approx(c(ref_val, max(data[i, is_outlier & data[i, ] >= upper_threshold])),
                 c(interp_val, upper_threshold),
                 xout = data[i, is_outlier & data[i, ] >= upper_threshold])$y
      } else {
        ref_val <- min(data[i, is_outlier & data[i, ] <= lower_threshold])
        interp_val <- med - 3 * mad_val
        adjusted_data[i, is_outlier & data[i, ] <= lower_threshold] <-
          approx(c(ref_val, min(data[i, !is_outlier & data[i, ] > lower_threshold], na.rm = TRUE)),
                 c(lower_threshold, interp_val),
                 xout = data[i, is_outlier & data[i, ] <= lower_threshold])$y
      }
    }
  }
  return(adjusted_data)
}

#' Correct for Drift in Data Using Autocorrelation Correction
#'
#' This function corrects for drift effects in metabolomics data by detrending based on run order within each plate segment.
#'
#' @param data A numeric matrix with rows representing metabolites and columns representing samples.
#' @param run_order An optional numeric vector representing the run order of the samples.
#' @param plates A numeric vector indicating the plate (or segment) assignment for each sample.
#' @param lag An integer specifying the lag to be used in the autocorrelation test.
#' @param test A character string specifying the autocorrelation test to use ("Ljung-Box" or "DW").
#' @param residualize A character string indicating the method for detrending ("mean" or "spline").
#' @param fdr_threshold A numeric value specifying the FDR threshold for significance.
#' @return A numeric matrix with drift corrected.
#' @examples
#' your_data_matrix <- matrix(rnorm(200), nrow = 20)
#' plates <- rep(1:4, length.out = ncol(your_data_matrix))
#' run_order <- seq_len(ncol(your_data_matrix))
#' drift_corrected <- autocorrelation_correct(your_data_matrix, run_order, plates)
#' @export
autocorrelation_correct <- function(data,
                                    run_order = NULL,
                                    plates,
                                    lag = 20,
                                    test = "Ljung-Box",
                                    residualize = "mean",
                                    fdr_threshold = 0.05) {
  if (!is.matrix(data)) stop("Data must be a numeric matrix.")
  if (!is.null(run_order) && length(run_order) != ncol(data)) {
    stop("Length of run_order must match number of columns in data.")
  }
  if (is.null(plates) || length(plates) != ncol(data)) {
    stop("Plates vector must be provided and its length must equal number of columns in data.")
  }
  
  detrended_data <- data
  unique_plates <- unique(plates)
  
  for (plate in unique_plates) {
    idx <- which(plates == plate)
    segment <- data[, idx, drop = FALSE]
    seg_run <- if (!is.null(run_order)) run_order[idx] else seq_along(idx)
    if (residualize == "spline") {
      spline_vals <- apply(segment, 1, function(x) {
        fit <- tryCatch(mgcv::gam(x ~ s(seg_run, bs = "cr")), error = function(e) NULL)
        if (is.null(fit)) rep(mean(x, na.rm = TRUE), length(x)) else fit$fitted.values
      })
      detrended_data[, idx] <- data[, idx] - t(spline_vals)
    } else if (residualize == "mean") {
      p_vals <- apply(segment, 1, function(x) {
        if (test == "Ljung-Box") {
          tryCatch(Box.test(x, lag = lag, type = "Ljung-Box")$p.value, error = function(e) NA)
        } else if (test == "DW") {
          tryCatch(lmtest::dwtest(x ~ seg_run)$p.value, error = function(e) NA)
        } else {
          stop("Invalid test method. Use 'Ljung-Box' or 'DW'.")
        }
      })
      p_vals <- p.adjust(p_vals, method = "fdr")
      correct_ids <- which(p_vals < fdr_threshold)
      if (length(correct_ids) > 0) {
        spline_vals <- apply(segment[correct_ids, , drop = FALSE], 1, function(x) {
          fit <- tryCatch(mgcv::gam(x ~ s(seg_run, bs = "cr")), error = function(e) NULL)
          if (is.null(fit)) rep(mean(x, na.rm = TRUE), length(x)) else fit$fitted.values
        })
        detrended_data[correct_ids, idx] <- data[correct_ids, idx] - t(spline_vals)
      }
      non_correct_ids <- setdiff(1:nrow(data), correct_ids)
      if (length(non_correct_ids) > 0) {
        detrended_data[non_correct_ids, idx] <- data[non_correct_ids, idx] - rowMeans(segment[non_correct_ids,], na.rm = TRUE)
      }
    }
  }
  return(detrended_data)
}

#' Perform ANOVA on Metabolites Across Plates
#'
#' This function performs an ANOVA test for each metabolite to assess plate effects.
#'
#' @param data A numeric matrix with rows representing metabolites and columns representing samples.
#' @param plates A numeric vector indicating the plate (or segment) assignment for each sample.
#' @return A matrix with F-statistics, p-values, and FDR-adjusted p-values for each metabolite.
#' @examples
#' your_data_matrix <- matrix(rnorm(200), nrow = 20)
#' plates <- rep(1:4, length.out = ncol(your_data_matrix))
#' aov_results <- anova_test(your_data_matrix, plates)
#' @export
anova_test <- function(data, plates) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("Data must be a numeric matrix.")
  }
  if (length(plates) != ncol(data)) {
    stop("Length of plates must match number of columns in data.")
  }
  if (length(unique(plates)) < 2) {
    # If only one plate is available, ANOVA is not applicable.
    n_metabolites <- nrow(data)
    results <- matrix(NA, nrow = n_metabolites, ncol = 3)
    colnames(results) <- c("F-statistic", "p-value", "p-adj")
    return(results)
  }
  n_metabolites <- nrow(data)
  results <- matrix(NA, nrow = n_metabolites, ncol = 2)
  colnames(results) <- c("F-statistic", "p-value")
  
  for (i in seq_len(n_metabolites)) {
    df <- data.frame(value = data[i, ], group = as.factor(plates))
    fit <- aov(value ~ group, data = df)
    sfit <- summary(fit)
    results[i, ] <- c(sfit[[1]]$`F value`[1], sfit[[1]]$`Pr(>F)`[1])
  }
  padj <- p.adjust(results[, "p-value"], method = "fdr")
  results <- cbind(results, padj)
  colnames(results) <- c("F-statistic", "p-value", "p-adj")
  return(results)
}

#' Residualize Data by Plate Based on ANOVA Results
#'
#' This function adjusts metabolite intensities for significant plate effects identified via ANOVA.
#'
#' @param data A numeric matrix with rows representing metabolites and columns representing samples.
#' @param plates A numeric vector indicating the plate (or segment) assignment for each sample.
#' @param fdr_threshold A numeric value specifying the FDR threshold for significance.
#' @return A numeric matrix with residualized intensities.
#' @examples
#' your_data_matrix <- matrix(rnorm(200), nrow = 20)
#' plates <- rep(1:4, length.out = ncol(your_data_matrix))
#' residualized_data <- residualize_by_plate_with_anova(your_data_matrix, plates)
#' @export
residualize_by_plate_with_anova <- function(data, plates, fdr_threshold = 0.05) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("Data must be a numeric matrix.")
  }
  if (length(plates) != ncol(data)) {
    stop("Length of plates must match number of columns in data.")
  }
  if (length(unique(plates)) < 2) {
    message("Only one plate detected. Skipping ANOVA-based residualization.")
    return(data)
  }
  
  aov_res <- anova_test(data, plates)
  sig_rows <- which(aov_res[, "p-adj"] < fdr_threshold)
  residualized <- data
  if (length(sig_rows) > 0) {
    for (plate in unique(plates)) {
      idx <- which(plates == plate)
      for (i in sig_rows) {
        plate_mean <- mean(data[i, idx], na.rm = TRUE)
        residualized[i, idx] <- data[i, idx] - plate_mean
      }
    }
  }
  return(residualized)
}

#' Scale Data by Plate
#'
#' This function scales the values for each metabolite within each plate by subtracting the plate mean and dividing by the plate standard deviation.
#'
#' @param data A numeric matrix with rows representing metabolites and columns representing samples.
#' @param plates A numeric vector indicating the plate (or segment) assignment for each sample.
#' @return A numeric matrix of scaled intensities.
#' @examples
#' your_data_matrix <- matrix(rnorm(200), nrow = 20)
#' plates <- rep(1:4, length.out = ncol(your_data_matrix))
#' scaled_data <- scale_by_plate(your_data_matrix, plates)
#' @export
scale_by_plate <- function(data, plates) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("Data must be a numeric matrix.")
  }
  if (length(plates) != ncol(data)) {
    stop("Length of plates must match number of columns in data.")
  }
  
  scaled_data <- data
  unique_plates <- unique(plates)
  for (plate in unique_plates) {
    idx <- which(plates == plate)
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
#' If plate information is not supplied, segments are automatically detected using an fkPELT-based approach and labeled as plates.
#'
#' The correction pipeline is as follows:
#' \enumerate{
#'   \item Median adjustment by dilution factor.
#'   \item Outlier adjustment using MAD.
#'   \item Drift correction using autocorrelation-based detrending.
#'   \item Batch effect correction using ANOVA-based residualization.
#'   \item Optional scaling by plate.
#' }
#'
#' When control samples are provided, the function can also auto-detect the optimal parameter settings (for plate detection,
#' autocorrelation test type and FDR threshold, ANOVA FDR threshold, and scaling) by selecting the combination that maximizes
#' the mean correlation between control samples.
#'
#' @param data A numeric matrix or data frame where rows represent metabolites and columns represent samples.
#' @param plates An optional numeric vector indicating plate assignments for each sample. If NULL, segments will be auto-detected.
#' @param run_order An optional numeric vector representing the run order of samples.
#' @param control_samples An optional numeric vector representing the columns corresponding to control samples. If provided,
#' these will be used for normalization and parameter tuning.
#' @param parameters Character string specifying whether to use fixed ("fixed") or auto-detected ("auto") parameters in presence of control samples.
#' @param fdr_threshold A numeric value specifying the FDR threshold for drift and batch corrections (default: 0.05).
#' @param median_adjustment A character string specifying the method for median adjustment ("shrink", "normalize", or "none").
#' @param residualize_non_autocorrelated A character string specifying the method for drift correction ("mean" or "spline").
#' @param test A character string specifying the autocorrelation test ("Ljung-Box" or "DW").
#' @param lag An integer specifying the lag for the autocorrelation test.
#' @param scale_by_plate Logical indicating whether to scale data by plate after corrections.
#' @return A numeric matrix of corrected intensities.
#' @examples
#' your_data_matrix <- matrix(rnorm(200), nrow = 20)
#' plates <- rep(1:4, length.out = ncol(your_data_matrix))
#' run_order <- seq_len(ncol(your_data_matrix))
#' corrected_data <- winn(your_data_matrix, plates = plates, run_order = run_order)
#' @export
winn <- function(data,
                 plates = NULL,
                 run_order = NULL,
                 control_samples = NULL,
                 parameters = c("fixed", "auto"),
                 fdr_threshold = 0.05,
                 median_adjustment = "shrink",
                 residualize_non_autocorrelated = "mean",
                 test = "Ljung-Box",
                 lag = 20,
                 scale_by_plate = FALSE) {
  parameters <- match.arg(parameters)
  
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.numeric(data)) stop("Data must be numeric.")
  n_samples <- ncol(data)
  if (!is.null(run_order) && length(run_order) != n_samples) {
    stop("Length of run_order must match number of columns in data.")
  }
  if (!is.null(control_samples)){
    if (any(is.na(match(control_samples, 1:n_samples)))) {
      stop("Control samples should refer to column numbers in the data matrix.")
    }
  }
  
  message("Starting Winn correction...")
  
  # Median adjustment using control samples if provided
  if (median_adjustment == "none") {
    message("Skipping median adjustment.")
    norm_data <- data
  } else {
    message("Performing median adjustment using method: ", median_adjustment)
    norm_data <- normalize_by_dilution_factor(data, processing = median_adjustment, control_samples = control_samples)
  }
  
  # Outlier adjustment
  message("Adjusting outliers using MAD...")
  norm_data <- adjust_outliers_mad(norm_data)
  
  # If control samples are provided and auto parameter detection is enabled, perform grid search
  if (!is.null(control_samples) && parameters == "auto") {
    message("Auto-detecting optimal parameters using control samples...")
    
    # Define grid for parameter search
    plate_options <- if (!is.null(plates)) c("provided") else c("auto")
    tests <- c("Ljung-Box", "DW")
    acorr_fdr_options <- c(0.1, 0.05, 0.01)
    anova_fdr_options <- c(0.1, 0.05, 0.01)
    scale_options <- c(TRUE, FALSE)
    
    best_mean_corr <- -Inf
    best_final_data <- NULL
    best_params <- list()
    
    for (plate_option in plate_options) {
      print(paste0("trying plate_option = ", plate_option))
      current_plates <- if (plate_option == "provided") {
        plates
      } else {
        .auto_detect_plates(norm_data)
      }
      
      for (current_test in tests) {
        print(paste0("trying test = ", current_test))
        for (acorr_fdr in acorr_fdr_options) {
          print(paste0("trying fdr = ", acorr_fdr))
          # Drift correction with current autocorrelation parameters
          drift_corrected <- autocorrelation_correct(norm_data,
                                                     run_order = run_order,
                                                     plates = current_plates,
                                                     lag = lag,
                                                     test = current_test,
                                                     residualize = residualize_non_autocorrelated,
                                                     fdr_threshold = acorr_fdr)
          for (anova_fdr in anova_fdr_options) {
            print(paste0("trying anova fdr = ", anova_fdr))
            # Batch effect correction with current ANOVA FDR
            batch_corrected <- residualize_by_plate_with_anova(drift_corrected, current_plates, fdr_threshold = anova_fdr)
            for (scale_opt in scale_options) {
              print(paste0("trying scale_option = ", scale_opt))
              final_data <- if (scale_opt) {
                scale_by_plate(batch_corrected, current_plates)
              } else {
                batch_corrected
              }
              
              mean_corr <- .mean_control_correlation(final_data, control_samples)
              print(paste0("mean correlation within controls = ", mean_corr))
              if (!is.na(mean_corr) && mean_corr > best_mean_corr) {
                best_mean_corr <- mean_corr
                best_final_data <- final_data
                best_params <- list(plate_option = plate_option,
                                    test = current_test,
                                    acorr_fdr = acorr_fdr,
                                    anova_fdr = anova_fdr,
                                    scale_by_plate = scale_opt)
              }
            }
          }
        }
      }
    }
    
    message("Optimal parameters selected based on control samples:")
    message("Plate option: ", best_params$plate_option)
    message("Autocorrelation test: ", best_params$test, " with FDR threshold: ", best_params$acorr_fdr)
    message("ANOVA FDR threshold: ", best_params$anova_fdr)
    message("Scale by plate: ", best_params$scale_by_plate)
    
    final_data <- best_final_data
  } else {
    # Use fixed parameters
    if (is.null(plates)) {
      message("No plate information provided. Auto-detecting segments using fkPELT...")
      plates <- .auto_detect_plates(norm_data)
      message("Auto-detected ", max(plates), " segments as plates.")
    } else {
      if (length(plates) != n_samples) {
        stop("Length of plates must match number of columns in data.")
      }
    }
    
    # Drift correction
    message("Correcting drift using autocorrelation test (", test, ") with FDR threshold: ", fdr_threshold, "...")
    drift_corrected <- autocorrelation_correct(norm_data,
                                               run_order = run_order,
                                               plates = plates,
                                               lag = lag,
                                               test = test,
                                               residualize = residualize_non_autocorrelated,
                                               fdr_threshold = fdr_threshold)
    
    # Batch effect correction via ANOVA
    message("Correcting batch effects using ANOVA-based residualization with FDR threshold: ", fdr_threshold, "...")
    batch_corrected <- residualize_by_plate_with_anova(drift_corrected, plates, fdr_threshold = fdr_threshold)
    
    # Optional scaling by plate
    if (scale_by_plate) {
      message("Scaling data by plate...")
      final_data <- scale_by_plate(batch_corrected, plates)
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

.auto_detect_plates <- function(data) {
  # Auto-detect segments as plates using fkPELT based on aggregated median signal
  agg_signal <- apply(data, 2, median, na.rm = TRUE)
  n <- length(agg_signal)
  num_knots <- floor(n / 60) + 2
  if (num_knots > 2) {
    knots <- seq(1, n, length.out = num_knots)
    knots <- knots[-c(1, length(knots))]
    knots <- ifelse(floor(knots) == knots, knots + 0.5, knots)
  } else {
    knots <- numeric(0)
  }
  change_points <- .fkPELT(agg_signal, knots)
  tau <- c(0, change_points, n)
  plates <- rep(NA, n)
  for (i in seq_along(tau[-length(tau)])) {
    plates[(tau[i] + 1):tau[i + 1]] <- i
  }
  return(plates)
}

.mean_control_correlation <- function(data, control_samples) {
  # Calculate mean pairwise correlation among control sample columns
  control_data <- data[, control_samples, drop = FALSE]
  if(ncol(control_data) < 2) return(NA)
  corr_matrix <- cor(control_data, use = "pairwise.complete.obs")
  lower_tri <- corr_matrix[lower.tri(corr_matrix)]
  return(mean(lower_tri, na.rm = TRUE))
}

.fkPELT <- function(data, knots) {
  if (is.null(data)) stop("Data cannot be NULL in .fkPELT")
  n <- length(data)
  penalty <- 3 * log(300) # 3*log(300) is the penalty constant that we use
  # in the PELT algorithm. Inspired by BIC criterion but modified.
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
      neglog[r] <- .fksplinecost(seg_data, knots, index1 = R[[t]][r] + 1, index2 = t)
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
      tot <- f[R[[t]][r] + 1] + log(300) + neglog[r]
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

.fksplinecost <- function(data, knots, index1 = 1, index2 = length(data)) {
  size <- length(data)
  if (size == 1) return(0)
  if (size < 5) {
    sd_val <- sd(data, na.rm = TRUE)
    mu_val <- mean(data, na.rm = TRUE)
    if (sd_val <= 1e-5) return(0)
    neglog <- 2 * (sum((data - mu_val)^2 / (2 * sd_val^2)) + size * log(sd_val * sqrt(2 * pi)))
    return(neglog)
  }
  cov <- index1:index2
  newknots <- knots[knots < index2 & knots > index1]
  fit <- tryCatch(lm(data ~ ns(cov, knots = newknots, intercept = TRUE)),
                  error = function(e) NULL)
  if (is.null(fit)) return(Inf)
  mu_val <- fit$fitted.values
  sd_val <- sd(data - mu_val, na.rm = TRUE)
  if (sd_val <= 1e-5) return(0)
  neglog <- 2 * (sum((data - mu_val)^2 / (2 * sd_val^2)) + size * log(sd_val * sqrt(2 * pi)))
  return(neglog)
}

###############################################################################
# Unit Tests (for interactive sessions)
###############################################################################
if (interactive()) {
  library(testthat)
  
  test_that("normalize_by_dilution_factor works (with and without control samples)", {
    set.seed(1)
    mat <- matrix(rnorm(100), nrow = 10)
    res1 <- normalize_by_dilution_factor(mat)
    res2 <- normalize_by_dilution_factor(mat, control_samples = 1:5)
    expect_equal(dim(res1), dim(mat))
    expect_equal(dim(res2), dim(mat))
  })
  
  test_that("winn works with auto-detected plates", {
    set.seed(1)
    mat <- matrix(rnorm(200), nrow = 20)
    res <- winn(mat)
    expect_equal(dim(res), dim(mat))
  })
  
  test_that("winn works with provided plates", {
    set.seed(1)
    mat <- matrix(rnorm(400), nrow = 20)
    plates <- rep(1:4, each = 5)
    run_order <- seq_len(ncol(mat))
    res <- winn(mat, plates = plates, run_order = run_order)
    expect_equal(dim(res), dim(mat))
  })
  
  test_that("winn auto detection with control samples works", {
    set.seed(1)
    mat <- matrix(rnorm(400), nrow = 20)
    plates <- rep(1:4, each = 5)
    run_order <- seq_len(ncol(mat))
    control_samples <- 1:4
    res <- winn(mat, plates = plates, run_order = run_order, control_samples = control_samples, parameters = "auto")
    expect_equal(dim(res), dim(mat))
    # Optionally, one might test that the mean correlation among controls is above a minimal threshold.
  })
}
