test_that("upper-only MAD outliers are shrunk to the upper threshold", {
  input <- matrix(
    c(10, 11, 9, 10, 10, 100),
    nrow = 1L,
    dimnames = list("upper_feature", paste0("sample_", seq_len(6L)))
  )
  centre <- median(input[1L, ])
  spread <- mad(input[1L, ])

  adjusted <- adjust_outliers_mad(input)

  expect_equal(adjusted[1L, 6L], centre + 4 * spread)
  expect_identical(adjusted[1L, seq_len(5L)], input[1L, seq_len(5L)])
})

test_that("lower-only MAD outliers are shrunk to the lower threshold", {
  input <- matrix(
    c(-100, 10, 11, 9, 10, 10),
    nrow = 1L,
    dimnames = list("lower_feature", paste0("sample_", seq_len(6L)))
  )
  centre <- median(input[1L, ])
  spread <- mad(input[1L, ])

  adjusted <- adjust_outliers_mad(input)

  expect_equal(adjusted[1L, 1L], centre - 4 * spread)
  expect_identical(adjusted[1L, 2:6], input[1L, 2:6])
})

test_that("simultaneous MAD tails are shrunk independently", {
  input <- matrix(
    c(-100, 9, 10, 10, 11, 100),
    nrow = 1L,
    dimnames = list("two_tailed_feature", paste0("sample_", seq_len(6L)))
  )
  centre <- median(input[1L, ])
  spread <- mad(input[1L, ])

  adjusted <- adjust_outliers_mad(input)

  expect_equal(adjusted[1L, 1L], centre - 4 * spread)
  expect_equal(adjusted[1L, 6L], centre + 4 * spread)
  expect_identical(adjusted[1L, 2:5], input[1L, 2:5])
})

test_that("both tails use thresholds computed once from original values", {
  input <- matrix(
    c(-100, -50, 9, 9, 10, 10, 10, 11, 11, 50, 100, NA_real_),
    nrow = 1L,
    dimnames = list("one_pass_feature", paste0("sample_", seq_len(12L)))
  )
  original_finite <- input[1L, is.finite(input[1L, ])]
  centre <- median(original_finite)
  spread <- mad(original_finite)
  lower_threshold <- centre - 4 * spread
  upper_threshold <- centre + 4 * spread

  adjusted <- adjust_outliers_mad(input)

  expect_equal(adjusted[1L, 1L], lower_threshold)
  expect_equal(adjusted[1L, 11L], upper_threshold)
  expect_true(all(adjusted[1L, 1:2] <= centre - 3 * spread))
  expect_true(all(adjusted[1L, 10:11] >= centre + 3 * spread))
  expect_identical(adjusted[1L, 3:9], input[1L, 3:9])
  expect_true(is.na(adjusted[1L, 12L]))
})

test_that("zero-MAD and non-finite-only features pass through unchanged", {
  input <- rbind(
    constant = rep(5, 6L),
    missing = rep(NA_real_, 6L),
    mixed_nonfinite = c(NA_real_, Inf, -Inf, NA_real_, Inf, -Inf)
  )

  expect_silent(adjusted <- adjust_outliers_mad(input))
  expect_identical(adjusted, input)
})

test_that("independent MAD adjustment preserves unflagged data and structure", {
  input <- rbind(
    upper = c(10, 11, 9, 10, 10, 100),
    lower = c(-100, 10, 11, 9, 10, 10),
    both = c(-100, 9, 10, 10, 11, 100),
    unchanged = c(2, 3, 4, 5, 6, 7)
  )
  colnames(input) <- paste0("sample_", c(6, 2, 5, 1, 4, 3))

  flagged <- t(apply(input, 1L, function(values) {
    centre <- median(values)
    spread <- mad(values)
    values <= centre - 4 * spread | values >= centre + 4 * spread
  }))
  adjusted <- adjust_outliers_mad(input)

  expect_identical(dim(adjusted), dim(input))
  expect_identical(dimnames(adjusted), dimnames(input))
  expect_identical(adjusted[!flagged], input[!flagged])
  expect_identical(adjusted["unchanged", ], input["unchanged", ])
})
