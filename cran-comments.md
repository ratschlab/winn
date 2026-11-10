# CRAN submission: winn 0.1.6

This is a patch release of winn. Version 0.1.5 has been on CRAN since
2026-09-27, and its CRAN checks are currently OK on all flavors. Thank you for
reviewing the update.

## Reason for the update

Version 0.1.6 fixes a bug that silently disabled spline smoothing in drift
correction, which is the package's main correction step.

- Drift-smoothing formulas wrote the smooth term as `mgcv::s()`. mgcv's formula
  parser does not recognise the namespace-qualified call as a smooth term, so
  every spline fit errored. The error was caught internally, and drift
  correction fell back to loess (the default conservative mode) or to a linear
  fit (standard mode) without warning.
- Smooth terms now use `s()`, so the documented splines are fitted.
- New regression tests require the documented smoother in both modes and for
  drift corrections made through `winn()`.

Results change only for drift-corrected features. No exported function or
argument changed.

## Test environments

- Local: macOS arm64 (aarch64-apple-darwin23), R 4.6.1, HTML Tidy 5.8.0.
  The source archive was built with `R CMD build .` and checked with
  `R CMD check --as-cran --timings winn_0.1.6.tar.gz`.
- GitHub Actions: Windows and macOS (R release), and Ubuntu (R release,
  R-devel, and R oldrel-1), using `R CMD check --no-manual --as-cran` with
  warnings treated as failures.

## R CMD check results

0 errors | 0 warnings | 0 notes

All 105 test expectations passed, with no failures, warnings, or skips.
Examples, including the additional `--run-donttest` pass, vignette rebuilding,
and the PDF and HTML manual checks also passed.

"Spectrometry" is retained in the title as the correct scientific term. If the
incoming spelling check flags it, it is a false positive.

## Reverse dependencies

There are no reverse dependencies on CRAN.
