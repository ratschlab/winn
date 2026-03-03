#' winn: White Noise Normalization for Metabolomics
#'
#' The `winn` package provides a transparent, stepwise workflow for correcting
#' technical artifacts in LC-MS and GC-MS metabolomics data. WINN is designed
#' for real study data that exhibit drift, batch effects, dilution variability,
#' and occasional outliers. Each correction step can be used independently or
#' combined in the full workflow via [winn()].
#'
#' Core steps in the WINN pipeline:
#'
#' 1. Robust outlier adjustment (MAD-based)
#' 2. Drift correction by run order with a white-noise gatekeeper
#' 3. Batch correction (ANOVA mean-shift or ComBat)
#' 4. Median or PQN normalization for dilution effects
#' 5. Optional per-batch scaling
#'
#' If batch labels are not provided, WINN can automatically detect batch
#' segments using an fkPELT change-point method on the median sample signal. If
#' pooled QC samples are provided, WINN can automatically tune several
#' parameters to maximize QC consistency.
#'
#' @section Data format:
#' WINN expects a metabolites-by-samples matrix. Run order, batch annotations,
#' and QC index vectors should all align with `ncol(data)`.
#'
#' @section Dependencies:
#' - `mgcv` for spline-based drift detrending
#' - `lmtest` for the Durbin-Watson autocorrelation test
#' - `sva` for ComBat batch correction
#'
#' @section Recommended usage:
#' - Provide `batch` and `run_order` whenever possible.
#' - Use QC samples and `parameters = "auto"` for robust tuning.
#' - Inspect QC CV and correlation before and after correction.
#'
#' @importFrom stats Box.test aov approx cor lm loess mad median model.matrix
#' @importFrom stats p.adjust sd var
#' @keywords internal
"_PACKAGE"
