test_that("ComBat preserves matrix structure and reduces simulated batch offsets", {
  skip_if_not_installed("sva")
  set.seed(5)
  batch <- rep(1:2, each = 6L)
  log_data <- matrix(rnorm(20L * 12L, mean = 5, sd = 0.1), nrow = 20L)
  log_data[, batch == 2L] <- log_data[, batch == 2L] + 1
  data <- expm1(log_data)
  dimnames(data) <- list(paste0("feature_", seq_len(nrow(data))),
                         paste0("sample_", seq_len(ncol(data))))

  corrected <- suppressMessages(combat_batch_correction(data, batch))

  expect_identical(dim(corrected), dim(data))
  expect_identical(dimnames(corrected), dimnames(data))
  expect_true(all(is.finite(corrected)))
  expect_true(all(corrected >= 0))
  batch_offset <- function(x) {
    mean(abs(rowMeans(log1p(x[, batch == 1L, drop = FALSE])) -
               rowMeans(log1p(x[, batch == 2L, drop = FALSE]))))
  }
  expect_lt(batch_offset(corrected), batch_offset(data) / 2)
})
