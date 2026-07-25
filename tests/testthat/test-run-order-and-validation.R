make_permutation_fixture <- function() {
  set.seed(20260726)
  n_features <- 12L
  n_samples <- 36L
  batch <- rep(letters[1:3], each = 12L)
  run_order <- seq_len(n_samples)
  log_data <- matrix(
    rnorm(n_features * n_samples, mean = 7, sd = 0.08),
    nrow = n_features
  )
  drift <- rep(seq(-0.3, 0.3, length.out = 12L), 3L)
  log_data[1:5, ] <- log_data[1:5, ] + matrix(
    drift,
    nrow = 5L,
    ncol = n_samples,
    byrow = TRUE
  )
  data <- pmax(expm1(log_data), 0)
  rownames(data) <- paste0("feature_", seq_len(n_features))
  colnames(data) <- paste0("sample_", seq_len(n_samples))
  list(
    data = data,
    batch = batch,
    run_order = run_order,
    controls = c(4L, 8L, 16L, 20L, 28L, 32L)
  )
}

test_that("drift correction is invariant to a joint sample permutation", {
  fixture <- make_permutation_fixture()
  set.seed(7)
  permutation <- sample(seq_len(ncol(fixture$data)))

  ordered <- autocorrelation_correct(
    fixture$data,
    run_order = fixture$run_order,
    batch = fixture$batch,
    lag = 4L,
    gate = "all",
    spline_method = "standard"
  )
  permuted <- autocorrelation_correct(
    fixture$data[, permutation, drop = FALSE],
    run_order = fixture$run_order[permutation],
    batch = fixture$batch[permutation],
    lag = 4L,
    gate = "all",
    spline_method = "standard"
  )

  expect_equal(permuted[, order(permutation), drop = FALSE], ordered, tolerance = 1e-10)
})

test_that("the full fixed pipeline is invariant to joint sample permutation", {
  fixture <- make_permutation_fixture()
  set.seed(9)
  permutation <- sample(seq_len(ncol(fixture$data)))
  permuted_controls <- match(fixture$controls, permutation)

  ordered <- suppressMessages(winn(
    fixture$data,
    batch = fixture$batch,
    run_order = fixture$run_order,
    control_samples = fixture$controls,
    parameters = "fixed",
    lag = 4L,
    spline_method = "standard"
  ))
  permuted <- suppressMessages(winn(
    fixture$data[, permutation, drop = FALSE],
    batch = fixture$batch[permutation],
    run_order = fixture$run_order[permutation],
    control_samples = permuted_controls,
    parameters = "fixed",
    lag = 4L,
    spline_method = "standard"
  ))

  expect_equal(permuted[, order(permutation), drop = FALSE], ordered, tolerance = 1e-10)
})

test_that("run order and batch metadata are validated before correction", {
  fixture <- make_permutation_fixture()
  duplicated_order <- fixture$run_order
  duplicated_order[2] <- duplicated_order[1]

  expect_error(
    autocorrelation_correct(
      fixture$data,
      run_order = duplicated_order,
      batch = fixture$batch
    ),
    "unique within each batch"
  )
  expect_error(
    autocorrelation_correct(
      fixture$data,
      run_order = replace(fixture$run_order, 1, NA_real_),
      batch = fixture$batch
    ),
    "finite numeric vector"
  )
  expect_error(
    autocorrelation_correct(
      fixture$data,
      run_order = fixture$run_order,
      batch = replace(fixture$batch, 1, NA_character_)
    ),
    "must not contain missing"
  )
})

test_that("log-based entry points reject negative and non-finite intensities", {
  fixture <- make_permutation_fixture()
  negative <- fixture$data
  negative[1, 1] <- -0.25
  nonfinite <- fixture$data
  nonfinite[1, 1] <- Inf

  expect_error(
    winn(negative, batch = fixture$batch, run_order = fixture$run_order),
    "non-negative"
  )
  expect_error(
    winn(nonfinite, batch = fixture$batch, run_order = fixture$run_order),
    "only finite"
  )
})

test_that("zero handling and dilution-factor failures are explicit", {
  data <- rbind(
    zero_reference = c(0, 0, 0, 0),
    quantitative = c(2, 4, 6, 8)
  )
  diagnostics <- normalize_by_dilution_factor(
    data,
    processing = "normalize",
    return_diagnostics = TRUE
  )
  expect_equal(unique(diagnostics$diagnostics$zero_reference_features_omitted), 1)
  expect_equal(diagnostics$data[1, ], rep(0, 4))

  zero_sample <- matrix(c(1, 2, 0, 0), nrow = 2)
  expect_error(
    normalize_by_dilution_factor(zero_sample, processing = "normalize"),
    "finite and positive"
  )
  expect_error(
    normalize_by_dilution_factor(data, processing = "invalid"),
    "arg"
  )
})

test_that("auto mode requires controls and details are machine-readable", {
  fixture <- make_permutation_fixture()
  expect_error(
    winn(
      fixture$data,
      batch = fixture$batch,
      run_order = fixture$run_order,
      parameters = "auto"
    ),
    "requires at least two control"
  )

  result <- suppressMessages(winn(
    fixture$data,
    batch = fixture$batch,
    run_order = fixture$run_order,
    control_samples = fixture$controls,
    parameters = "fixed",
    lag = 4L,
    spline_method = "standard",
    return_details = TRUE
  ))
  expect_s3_class(result, "winn_result")
  expect_identical(dim(result$data), dim(fixture$data))
  expect_equal(result$selected_parameters$ljung_box_fitdf, 0L)
  expect_named(result$stage_decisions, c("outlier", "drift", "batch", "dilution"))
  expect_true(all(result$stage_decisions$drift$ljung_box_fitdf == 0L))
})

test_that("scaling rejects undefined within-batch standard deviations", {
  data <- matrix(seq_len(12), nrow = 3)
  expect_error(
    scale_by_batch(data, batch = c("singleton", "b", "b", "b")),
    "at least two samples"
  )

  constant <- matrix(5, nrow = 2, ncol = 4)
  scaled <- scale_by_batch(constant, batch = rep(c("a", "b"), each = 2))
  expect_equal(scaled, matrix(0, nrow = 2, ncol = 4))
})
