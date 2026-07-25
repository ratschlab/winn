test_that("standard spline correction records complete diagnostics", {
  set.seed(20260723)
  n_features <- 6L
  n_samples <- 36L
  batch <- rep(1:3, each = 12L)
  run_order <- seq_len(n_samples)
  log_data <- matrix(
    rnorm(n_features * n_samples, mean = 8, sd = 0.05),
    nrow = n_features
  )
  within_batch_drift <- rep(seq(-0.5, 0.5, length.out = 12L), 3L)
  log_data <- log_data + matrix(
    within_batch_drift,
    nrow = n_features,
    ncol = n_samples,
    byrow = TRUE
  )
  data <- pmax(expm1(log_data), 0)
  rownames(data) <- paste0("feature_", seq_len(n_features))
  colnames(data) <- paste0("sample_", seq_len(n_samples))

  result <- expect_no_error(autocorrelation_correct(
    data,
    run_order = run_order,
    batch = batch,
    detrend = "spline",
    spline_method = "standard",
    gate = "all",
    return_diagnostics = TRUE
  ))

  expect_identical(dim(result$data), dim(data))
  expect_identical(dimnames(result$data), dimnames(data))
  expect_true(all(is.finite(result$data)))
  expect_true(all(result$data >= 0))

  corrected <- result$diagnostics$actual_correction
  expect_true(any(corrected))
  expect_true(all(!is.na(result$diagnostics$smoother_type[corrected])))
  expect_true(all(nzchar(result$diagnostics$fallback_path[corrected])))
  expect_true(all(is.finite(result$diagnostics$basis_dimension[corrected])))
  expect_true(all(result$diagnostics$basis_dimension[corrected] >= 2))
  expect_true(all(!is.na(result$diagnostics$warning_status)))
  expect_true(all(!is.na(result$diagnostics$warning_message)))
})
