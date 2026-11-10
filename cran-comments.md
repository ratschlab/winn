# CRAN submission: winn 0.1.6

This is a patch release of winn, which is currently on CRAN as version 0.1.5.

## Changes

- Fixes a bug in drift correction. Drift-smoothing formulas wrote the smooth
  term as `mgcv::s()`. mgcv's formula parser does not recognise the
  namespace-qualified call as a smooth term, so every spline fit errored. The
  error was caught internally, and drift correction silently fell back to
  loess (default conservative mode) or to a linear fit (standard mode).
  Smooth terms now use `s()`, and the documented splines are fitted.
- Adds regression tests for both smoother modes and for `winn()`.

Results change only for features that are drift-corrected. No exported
interface changed.

## Test environments

- local macOS (aarch64), R 4.6: `R CMD check --as-cran`
- GitHub Actions: macOS release, Windows release, Ubuntu devel, release and
  oldrel-1

## R CMD check results

0 errors | 0 warnings | 0 notes (local `--as-cran`; see the pull request for
the CI results)

## Reverse dependencies

There are no reverse dependencies on CRAN.
