#' Run WiNN as Explicit Ablation Stages
#'
#' A thin, fixed-parameter orchestration layer around the package's existing
#' outlier, drift, batch, and dilution-factor functions. It is intended for
#' controlled mechanistic ablations and returns every intermediate matrix and
#' gate diagnostic without changing the standard [winn()] API.
#'
#' @param data Numeric feature-by-sample matrix or data frame.
#' @param batch Supplied batch or segment label for every sample.
#' @param run_order Numeric acquisition order for every sample.
#' @param control_samples Optional control-sample column indices. Use `NULL` for
#' QC-agnostic ablations.
#' @param parameters Must be `"fixed"`; automatic tuning is deliberately not
#' available through this interface.
#' @param use_outlier_shrinkage Whether to apply [adjust_outliers_mad()].
#' @param drift_gate One of `"none"`, `"selective"`, or `"all"`.
#' @param batch_gate One of `"none"`, `"selective"`, or `"all"`.
#' @param pqn_mode One of `"none"` or `"shrink"`.
#' @param fdr_threshold Fixed FDR threshold for drift and batch gates.
#' @param test Autocorrelation test, normally `"Ljung-Box"`.
#' @param lag Optional fixed lag; `NULL` uses the package's adaptive lag.
#' @param spline_method Drift smoother, normally `"conservative"`.
#' @param return_intermediates Whether to return all stage matrices.
#' @param return_diagnostics Whether to return gate and magnitude diagnostics.
#'
#' @return A list containing the final matrix, intermediate matrices,
#' diagnostics, stage runtimes, and the exact fixed configuration.
#' @export
winn_ablation <- function(data,
                          batch,
                          run_order = NULL,
                          control_samples = NULL,
                          parameters = "fixed",
                          use_outlier_shrinkage = TRUE,
                          drift_gate = c("selective", "none", "all"),
                          batch_gate = c("selective", "none", "all"),
                          pqn_mode = c("shrink", "none"),
                          fdr_threshold = 0.05,
                          test = "Ljung-Box",
                          lag = NULL,
                          spline_method = "conservative",
                          return_intermediates = TRUE,
                          return_diagnostics = TRUE) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame.")
  }
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    stop("Data must be numeric.")
  }
  if (!identical(parameters, "fixed")) {
    stop("winn_ablation() supports fixed parameters only.")
  }
  if (missing(batch) || is.null(batch) || length(batch) != ncol(data)) {
    stop("A supplied batch vector matching the sample count is required.")
  }
  if (!is.null(run_order) && length(run_order) != ncol(data)) {
    stop("Length of run_order must match number of columns in data.")
  }
  if (!is.null(control_samples) &&
      any(is.na(match(control_samples, seq_len(ncol(data)))))) {
    stop("control_samples must contain valid column indices.")
  }
  if (any(is.na(data) | is.infinite(data))) {
    warning("Data contains NA or infinite values. Results may be unreliable.")
  }
  drift_gate <- match.arg(drift_gate)
  batch_gate <- match.arg(batch_gate)
  pqn_mode <- match.arg(pqn_mode)

  raw <- data
  stage_runtimes <- stats::setNames(rep(0, 4L), c("outlier", "drift", "batch", "pqn"))

  outlier_start <- proc.time()[["elapsed"]]
  post_outlier <- if (isTRUE(use_outlier_shrinkage)) {
    adjust_outliers_mad(raw)
  } else {
    raw
  }
  stage_runtimes[["outlier"]] <- proc.time()[["elapsed"]] - outlier_start

  outlier_delta <- abs(log1p(post_outlier) - log1p(raw))
  outlier_changed <- outlier_delta > 1e-12
  outlier_diagnostics <- list(
    entries_changed = sum(outlier_changed, na.rm = TRUE),
    proportion_entries_changed = mean(outlier_changed, na.rm = TRUE),
    features_changed = sum(rowSums(outlier_changed, na.rm = TRUE) > 0L),
    samples_changed = sum(colSums(outlier_changed, na.rm = TRUE) > 0L),
    median_absolute_log1p_change = stats::median(outlier_delta, na.rm = TRUE),
    p90_absolute_log1p_change = as.numeric(stats::quantile(
      outlier_delta,
      probs = 0.9,
      na.rm = TRUE,
      names = FALSE
    ))
  )

  drift_start <- proc.time()[["elapsed"]]
  drift_result <- autocorrelation_correct(
    post_outlier,
    run_order = run_order,
    batch = batch,
    lag = lag,
    test = test,
    detrend = "mean",
    fdr_threshold = fdr_threshold,
    spline_method = spline_method,
    gate = drift_gate,
    return_diagnostics = TRUE
  )
  post_drift <- drift_result$data
  stage_runtimes[["drift"]] <- proc.time()[["elapsed"]] - drift_start

  batch_start <- proc.time()[["elapsed"]]
  batch_result <- anova_batch_correction(
    post_drift,
    batch = batch,
    fdr_threshold = fdr_threshold,
    gate = batch_gate,
    return_diagnostics = TRUE
  )
  post_batch <- batch_result$data
  stage_runtimes[["batch"]] <- proc.time()[["elapsed"]] - batch_start

  pqn_start <- proc.time()[["elapsed"]]
  if (identical(pqn_mode, "none")) {
    final <- post_batch
    sample_ids <- colnames(post_batch)
    if (is.null(sample_ids)) {
      sample_ids <- paste0("sample_", seq_len(ncol(post_batch)))
    }
    pqn_diagnostics <- data.frame(
      sample_id = sample_ids,
      dilution_factor = NA_real_,
      lower_threshold = NA_real_,
      upper_threshold = NA_real_,
      altered = FALSE,
      scaling_factor = 1,
      reference_spectrum_type = "not_applied",
      processing = "none",
      stringsAsFactors = FALSE
    )
  } else {
    pqn_result <- normalize_by_dilution_factor(
      post_batch,
      processing = "shrink",
      control_samples = control_samples,
      return_diagnostics = TRUE
    )
    final <- pqn_result$data
    pqn_diagnostics <- pqn_result$diagnostics
  }
  stage_runtimes[["pqn"]] <- proc.time()[["elapsed"]] - pqn_start

  dimnames(final) <- dimnames(raw)
  intermediates <- list(
    raw = raw,
    post_outlier = post_outlier,
    post_drift = post_drift,
    post_batch = post_batch,
    final = final
  )
  diagnostics <- list(
    outlier = outlier_diagnostics,
    drift = drift_result$diagnostics,
    batch = batch_result$diagnostics,
    pqn = pqn_diagnostics,
    matrix_state = do.call(rbind, lapply(names(intermediates), function(stage) {
      value <- intermediates[[stage]]
      data.frame(
        stage = stage,
        n_features = nrow(value),
        n_samples = ncol(value),
        na_count = sum(is.na(value)),
        nan_count = sum(is.nan(value)),
        infinite_count = sum(is.infinite(value)),
        minimum = suppressWarnings(min(value, na.rm = TRUE)),
        maximum = suppressWarnings(max(value, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))
  )

  list(
    data = final,
    intermediates = if (isTRUE(return_intermediates)) intermediates else NULL,
    diagnostics = if (isTRUE(return_diagnostics)) diagnostics else NULL,
    stage_runtime_sec = stage_runtimes,
    total_runtime_sec = sum(stage_runtimes),
    configuration = list(
      parameters = "fixed",
      supplied_batch = TRUE,
      control_samples_supplied = !is.null(control_samples),
      use_outlier_shrinkage = isTRUE(use_outlier_shrinkage),
      drift_gate = drift_gate,
      batch_gate = batch_gate,
      pqn_mode = pqn_mode,
      fdr_threshold = fdr_threshold,
      test = test,
      lag = if (is.null(lag)) "adaptive" else as.integer(lag),
      spline_method = spline_method,
      remove_batch_effects = "anova",
      scale_by_batch = FALSE
    )
  )
}
