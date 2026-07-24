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

test_that("simultaneous MAD tails use the intentional upper-priority policy", {
  input <- matrix(
    c(-100, 9, 10, 10, 11, 100),
    nrow = 1L,
    dimnames = list("two_tailed_feature", paste0("sample_", seq_len(6L)))
  )
  centre <- median(input[1L, ])
  spread <- mad(input[1L, ])

  adjusted <- adjust_outliers_mad(input)

  expect_equal(adjusted[1L, 6L], centre + 4 * spread)
  expect_identical(adjusted[1L, 1L], input[1L, 1L])
  expect_identical(adjusted[1L, 2:5], input[1L, 2:5])
})

test_that("zero and undefined MAD retain their established edge behavior", {
  constant <- matrix(
    rep(5, 6L),
    nrow = 1L,
    dimnames = list("constant_feature", paste0("sample_", seq_len(6L)))
  )
  expect_warning(
    constant_adjusted <- adjust_outliers_mad(constant),
    "no non-missing arguments to max"
  )
  expect_identical(constant_adjusted, constant)

  undefined <- matrix(
    rep(NA_real_, 6L),
    nrow = 1L,
    dimnames = list("missing_feature", paste0("sample_", seq_len(6L)))
  )
  expect_error(
    adjust_outliers_mad(undefined),
    "missing value where TRUE/FALSE needed"
  )
})

test_that("directional MAD adjustment preserves unflagged data and matrix structure", {
  input <- rbind(
    upper = c(10, 11, 9, 10, 10, 100),
    lower = c(-100, 10, 11, 9, 10, 10),
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
