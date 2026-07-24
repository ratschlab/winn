make_ablation_fixture <- function() {
  set.seed(42)
  n_features <- 18L
  n_samples <- 48L
  batch <- rep(1:3, each = 16L)
  run_order <- seq_len(n_samples)
  base <- matrix(rnorm(n_features * n_samples, 8, 0.12), nrow = n_features)
  drift <- rep(seq(-0.35, 0.35, length.out = 16L), 3L)
  base[1:6, ] <- base[1:6, ] + matrix(
    drift,
    nrow = 6L,
    ncol = n_samples,
    byrow = TRUE
  )
  base[7:12, batch == 2L] <- base[7:12, batch == 2L] + 0.45
  base[7:12, batch == 3L] <- base[7:12, batch == 3L] - 0.30
  data <- pmax(expm1(base), 0)
  rownames(data) <- paste0("feature_", seq_len(n_features))
  colnames(data) <- paste0("sample_", seq_len(n_samples))
  list(data = data, batch = batch, run_order = run_order)
}

test_that("complete ablation reproduces fixed winn without controls", {
  fixture <- make_ablation_fixture()
  standard <- suppressMessages(winn(
    fixture$data,
    batch = fixture$batch,
    run_order = fixture$run_order,
    control_samples = NULL,
    parameters = "fixed",
    fdr_threshold = 0.05,
    median_adjustment = "shrink",
    detrend_non_autocorrelated = "mean",
    spline_method = "conservative",
    remove_batch_effects = "anova",
    test = "Ljung-Box",
    lag = NULL,
    scale_by_batch = FALSE
  ))
  staged <- winn_ablation(
    fixture$data,
    batch = fixture$batch,
    run_order = fixture$run_order,
    control_samples = NULL,
    parameters = "fixed",
    use_outlier_shrinkage = TRUE,
    drift_gate = "selective",
    batch_gate = "selective",
    pqn_mode = "shrink"
  )
  expect_equal(staged$data, standard, tolerance = 1e-12)
})

test_that("ablation preserves dimensions, identifiers, order, and non-negativity", {
  fixture <- make_ablation_fixture()
  result <- winn_ablation(
    fixture$data,
    batch = fixture$batch,
    run_order = fixture$run_order
  )
  expect_identical(dim(result$data), dim(fixture$data))
  expect_identical(rownames(result$data), rownames(fixture$data))
  expect_identical(colnames(result$data), colnames(fixture$data))
  expect_true(all(is.finite(result$data)))
  expect_true(all(result$data >= 0))
})

test_that("forced-all drift corrects every testable profile", {
  fixture <- make_ablation_fixture()
  outlier <- adjust_outliers_mad(fixture$data)
  result <- autocorrelation_correct(
    outlier,
    run_order = fixture$run_order,
    batch = fixture$batch,
    gate = "all",
    return_diagnostics = TRUE
  )
  diagnostics <- result$diagnostics
  expect_true(all(diagnostics$actual_correction[diagnostics$testable]))
  expect_true(all(diagnostics$decision_reason[diagnostics$testable] == "forced_all"))
  expect_equal(
    sum(diagnostics$actual_correction),
    sum(diagnostics$testable)
  )
})

test_that("forced-all batch uses mean-preserving residualization for every eligible feature", {
  fixture <- make_ablation_fixture()
  result <- anova_batch_correction(
    fixture$data,
    batch = fixture$batch,
    gate = "all",
    return_diagnostics = TRUE
  )
  expect_true(all(result$diagnostics$actual_correction[result$diagnostics$eligible]))
  corrected_log <- log1p(result$data)
  raw_log <- log1p(fixture$data)
  expect_equal(rowMeans(corrected_log), rowMeans(raw_log), tolerance = 1e-10)
})

test_that("disabled stages are exact pass-throughs", {
  fixture <- make_ablation_fixture()
  outlier <- adjust_outliers_mad(fixture$data)
  no_drift <- autocorrelation_correct(
    outlier,
    run_order = fixture$run_order,
    batch = fixture$batch,
    gate = "none",
    return_diagnostics = TRUE
  )
  expect_identical(no_drift$data, outlier)

  no_batch <- anova_batch_correction(
    no_drift$data,
    batch = fixture$batch,
    gate = "none",
    return_diagnostics = TRUE
  )
  expect_equal(no_batch$data, no_drift$data, tolerance = 1e-12)

  no_pqn <- winn_ablation(
    fixture$data,
    batch = fixture$batch,
    run_order = fixture$run_order,
    drift_gate = "selective",
    batch_gate = "selective",
    pqn_mode = "none"
  )
  expect_identical(no_pqn$data, no_pqn$intermediates$post_batch)
})

test_that("diagnostic counts reproduce the saved selection masks", {
  fixture <- make_ablation_fixture()
  result <- winn_ablation(
    fixture$data,
    batch = fixture$batch,
    run_order = fixture$run_order,
    drift_gate = "selective",
    batch_gate = "selective",
    pqn_mode = "shrink"
  )
  drift <- result$diagnostics$drift
  batch <- result$diagnostics$batch
  pqn <- result$diagnostics$pqn
  expect_equal(
    sum(drift$actual_correction),
    sum(drift$decision_reason == "selected")
  )
  expect_equal(
    sum(batch$actual_correction),
    sum(batch$decision_reason == "selected")
  )
  expect_equal(sum(pqn$altered), sum(pqn$scaling_factor != 1))
})
