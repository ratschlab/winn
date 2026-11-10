## Regression tests for smooth terms in drift correction (winn 0.1.6).
## winn <= 0.1.5 wrote smooth terms as `mgcv::s()` inside model formulas. mgcv's
## formula parser does not recognise the namespace-qualified call as a smooth,
## so every spline fit errored and drift correction silently fell back to loess
## (conservative mode) or to a linear fit (standard mode).

simulate_drift <- function(n_features = 6L, n_per_batch = 40L, n_batches = 2L, seed = 20261110L) {
  set.seed(seed)
  n_samples <- n_per_batch * n_batches
  run_order <- seq_len(n_samples)
  batch <- rep(seq_len(n_batches), each = n_per_batch)
  drift <- rep(sin(seq(0, pi, length.out = n_per_batch)), n_batches)
  log_data <- matrix(rnorm(n_features * n_samples, mean = 8, sd = 0.05), nrow = n_features)
  log_data <- log_data + matrix(0.6 * drift, nrow = n_features, ncol = n_samples, byrow = TRUE)
  data <- pmax(expm1(log_data), 0)
  rownames(data) <- paste0("feature_", seq_len(n_features))
  colnames(data) <- paste0("sample_", seq_len(n_samples))
  list(data = data, run_order = run_order, batch = batch)
}

test_that("conservative spline fits an mgcv P-spline instead of falling back", {
  set.seed(20261110)
  x <- seq_len(120)
  y <- sin(x / 20) + rnorm(length(x), sd = 0.2)
  fit <- .fit_conservative_spline(y, x)
  expect_identical(attr(fit, "winn_smoother_type"), "conservative_ps_reml")
  expect_identical(attr(fit, "winn_fallback_path"), "ps_reml_gamma_1.4")
  expect_s3_class(fit, "gam")
})

test_that("conservative drift correction never uses the loess or linear fallback on smooth drift", {
  sim <- simulate_drift()
  result <- autocorrelation_correct(
    sim$data,
    run_order = sim$run_order,
    batch = sim$batch,
    detrend = "spline",
    spline_method = "conservative",
    gate = "all",
    return_diagnostics = TRUE
  )
  corrected <- result$diagnostics$actual_correction
  expect_true(any(corrected))
  expect_true(all(result$diagnostics$smoother_type[corrected] == "conservative_ps_reml"))
  expect_false(any(grepl("loess|linear", result$diagnostics$fallback_path[corrected])))
})

test_that("standard drift correction fits an mgcv cubic regression spline", {
  sim <- simulate_drift()
  result <- autocorrelation_correct(
    sim$data,
    run_order = sim$run_order,
    batch = sim$batch,
    detrend = "spline",
    spline_method = "standard",
    gate = "all",
    return_diagnostics = TRUE
  )
  corrected <- result$diagnostics$actual_correction
  expect_true(any(corrected))
  expect_true(all(result$diagnostics$smoother_type[corrected] == "standard_gam_cr"))
})

test_that("winn() reports conservative splines for every drift-corrected feature", {
  sim <- simulate_drift(n_features = 8L)
  result <- suppressMessages(winn(
    log1p(sim$data),
    batch = sim$batch,
    run_order = sim$run_order,
    return_details = TRUE
  ))
  drift <- result$stage_decisions$drift
  corrected <- drift$actual_correction %in% TRUE
  expect_true(any(corrected))
  expect_true(all(drift$smoother_type[corrected] == "conservative_ps_reml"))
})
