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

test_that("vectorized Ljung-Box p-values match Box.test", {
  set.seed(11)
  seg <- matrix(abs(rnorm(60, mean = 20, sd = 3)), nrow = 6)
  lag <- 4L
  fitdf <- 1L

  expected <- apply(seg, 1, function(x) {
    Box.test(log1p(x), lag = lag, type = "Ljung-Box", fitdf = fitdf)$p.value
  })

  observed <- winn:::.compute_ljung_box_pvalues(
    segment = seg,
    segment_lag = lag,
    model_df = fitdf
  )

  expect_equal(observed, expected, tolerance = 1e-10)
})

test_that("threshold candidate materializers collapse duplicate states", {
  prepared_drift <- list(
    mode = "mean",
    data = matrix(c(1, 2, 3, 4), nrow = 2),
    batches = list(list(
      idx = 1:2,
      adjusted_p = c(0.009, 0.5),
      candidate_ids = 1L,
      detrended = matrix(c(10, 20), nrow = 1)
    ))
  )

  drift_candidates <- winn:::.materialize_autocorrelation_candidates(
    prepared_drift,
    thresholds = c(0.1, 0.05, 0.01)
  )
  expect_length(drift_candidates, 1)

  prepared_batch <- list(
    original = matrix(c(log1p(1), log1p(2), log1p(3), log1p(4)), nrow = 2),
    corrected_all = matrix(c(log1p(5), log1p(6), log1p(3), log1p(4)), nrow = 2),
    padj = c(0.009, 0.5)
  )

  batch_candidates <- winn:::.materialize_anova_candidates(
    prepared_batch,
    thresholds = c(0.1, 0.05, 0.01)
  )
  expect_length(batch_candidates, 1)
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

test_that("fixed mode requires a single test but auto mode can opt into DW", {
  set.seed(4)
  mat <- matrix(abs(rnorm(60, mean = 40, sd = 6)), nrow = 6)
  batch <- rep(1:2, each = 5)
  run_order <- seq_len(10)

  expect_error(
    winn(
      mat,
      batch = batch,
      run_order = run_order,
      test = c("Ljung-Box", "DW")
    ),
    "length 1"
  )

  corrected <- expect_no_error(
    suppressMessages(winn(
      mat,
      batch = batch,
      run_order = run_order,
      control_samples = c(1, 6),
      parameters = "auto",
      test = c("Ljung-Box", "DW")
    ))
  )

  expect_equal(dim(corrected), c(6, 10))
})

test_that("vectorized ANOVA correction matches the formula-based result", {
  set.seed(5)
  mat <- matrix(abs(rnorm(80, mean = 30, sd = 5)), nrow = 8)
  batch <- factor(rep(1:2, each = 5))

  reference <- {
    z <- log1p(mat)
    pvals <- vapply(seq_len(nrow(z)), function(i) {
      df <- data.frame(value = z[i, ], batch = batch)
      summary(aov(value ~ batch, data = df))[[1]]$`Pr(>F)`[1]
    }, numeric(1))
    padj <- p.adjust(pvals, method = "fdr")
    corrected <- z
    sig <- which(padj < 0.05)
    if (length(sig) > 0) {
      overall_means <- rowMeans(z, na.rm = TRUE)
      for (i in sig) {
        batch_means <- tapply(z[i, ], batch, mean, na.rm = TRUE)
        shift <- batch_means - overall_means[i]
        corrected[i, ] <- z[i, ] - shift[as.character(batch)]
      }
    }
    pmax(expm1(corrected), 0)
  }

  expect_equal(
    anova_batch_correction(mat, batch, fdr_threshold = 0.05),
    reference,
    tolerance = 1e-8
  )
})
