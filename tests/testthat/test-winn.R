test_that("adaptive lag enforces the heuristic and minimum degrees of freedom", {
  expect_equal(
    winn:::.resolve_autocorrelation_lag(lag = NULL, n_obs = 50, model_df = 1L),
    10L
  )
  expect_equal(
    winn:::.resolve_autocorrelation_lag(lag = NULL, n_obs = 15, model_df = 1L),
    4L
  )
  expect_equal(
    winn:::.resolve_autocorrelation_lag(lag = 2, n_obs = 20, model_df = 1L),
    4L
  )
  expect_true(is.na(
    winn:::.resolve_autocorrelation_lag(lag = NULL, n_obs = 4, model_df = 1L)
  ))
})

test_that("fkPELT penalties resolve for named strategies", {
  set.seed(1)
  mat <- matrix(abs(rnorm(48, mean = 100, sd = 5)), nrow = 6)
  agg <- apply(mat, 2, median, na.rm = TRUE)

  mbic <- winn:::.resolve_pelt_penalty(
    agg_signal = agg,
    pelt_penalty = "mbic"
  )
  expect_equal(mbic, 3 * log(length(agg)))

  bic <- winn:::.resolve_pelt_penalty(
    agg_signal = agg,
    pelt_penalty = "bic"
  )
  expect_true(is.numeric(bic) && length(bic) == 1L && bic > 0)
  expect_equal(winn:::.resolve_pelt_penalty(agg_signal = agg, pelt_penalty = NULL), mbic)
})

test_that("winn supports adaptive lag without replicate tuning inputs", {
  set.seed(2)
  mat <- matrix(abs(rnorm(120, mean = 50, sd = 10)), nrow = 12)
  batch <- rep(1:2, each = 5)

  corrected <- expect_no_error(
    suppressMessages(winn(
      mat[, 1:10, drop = FALSE],
      batch = batch,
      run_order = seq_len(10),
      lag = NULL
    ))
  )

  expect_equal(dim(corrected), c(12, 10))
  expect_true(all(is.finite(corrected)))
})

test_that("winn rejects the removed automatic penalty alias", {
  set.seed(3)
  mat <- matrix(abs(rnorm(60, mean = 25, sd = 4)), nrow = 6)

  expect_error(
    winn(
      mat,
      batch = rep(1:2, each = 5),
      run_order = seq_len(10),
      pelt_penalty = "auto"
    ),
    "pelt_penalty"
  )
})
