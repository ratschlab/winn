# winn

`winn` implements White Noise Normalization for metabolomics LC-MS data. The
package applies a decision-guided correction pipeline: each metabolite is first
tested for serial dependence, and only metabolites that fail the white-noise
gatekeeper are drift-corrected. This keeps stable features untouched while
still correcting drift, batch structure, dilution effects, and outliers when
they are present.

The package supports both QC-informed tuning and QC-agnostic workflows. If
batch labels are unavailable, `winn()` can infer batch-like segments with an
fkPELT-style change-point routine before downstream correction.

## Installation

The current release is installed from GitHub. The package imports the
Bioconductor package `sva` for the ComBat correction path.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("sva")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("ratschlab/winn")
```

After a CRAN release is published, `install.packages("winn")` will also work.

## Reproducible Example

The example below creates a small simulated dataset and runs the standard WiNN
wrapper with fixed parameters. It is intentionally lightweight and matches the
current exported API.

```r
library(winn)

set.seed(1)

n_met <- 20L
n_samples <- 48L
batch <- rep(seq_len(4), each = 12)
run_order <- seq_len(n_samples)
qc_idx <- seq(4, n_samples, by = 8)
study_idx <- setdiff(seq_len(n_samples), qc_idx)

base_log <- matrix(rnorm(n_met * n_samples, mean = 8, sd = 0.15), nrow = n_met)
true_log <- base_log

signal_metabolites <- 1:4
study_signal <- as.numeric(scale(sin(study_idx / 6) + rnorm(length(study_idx), sd = 0.2)))
true_log[signal_metabolites, study_idx] <- true_log[signal_metabolites, study_idx] +
  0.35 * matrix(
    study_signal,
    nrow = length(signal_metabolites),
    ncol = length(study_idx),
    byrow = TRUE
  )

pooled_qc <- rowMeans(true_log[, study_idx, drop = FALSE])
true_log[, qc_idx] <- pooled_qc

dilution <- exp(rnorm(n_samples, sd = 0.05))
drift <- rep(seq(-0.18, 0.18, length.out = 12), times = 4)
batch_shift <- rep(c(-0.12, 0.04, 0.10, -0.06), each = 12)
noise <- matrix(rnorm(n_met * n_samples, sd = 0.08), nrow = n_met)

observed_log <- true_log +
  matrix(log(dilution), nrow = n_met, ncol = n_samples, byrow = TRUE) +
  matrix(drift + batch_shift, nrow = n_met, ncol = n_samples, byrow = TRUE) +
  noise

observed_intensity <- pmax(expm1(observed_log), 0)

corrected_intensity <- winn(
  observed_intensity,
  batch = batch,
  run_order = run_order,
  control_samples = qc_idx,
  parameters = "fixed",
  remove_batch_effects = "anova",
  median_adjustment = "shrink",
  lag = NULL,
  scale_by_batch = FALSE
)

dim(corrected_intensity)
```

## Stepwise Pipeline

Each correction step is also available as a standalone function.

```r
step1 <- adjust_outliers_mad(observed_intensity)

step2 <- autocorrelation_correct(
  step1,
  run_order = run_order,
  batch = batch,
  lag = NULL,
  test = "Ljung-Box",
  detrend = "mean"
)

step3 <- anova_batch_correction(step2, batch = batch)

step4 <- normalize_by_dilution_factor(
  step3,
  processing = "shrink",
  control_samples = qc_idx
)

step5 <- scale_by_batch(step4, batch = batch)  # optional
```

When pooled QCs are available and you want the package to search across
parameter settings, use `parameters = "auto"`. That mode is slower because it
evaluates multiple combinations using the package's QC quality score.

```r
corrected_auto <- winn(
  observed_intensity,
  batch = batch,
  run_order = run_order,
  control_samples = qc_idx,
  parameters = "auto",
  lag = NULL,
  scale_by_batch = FALSE
)
```

If batch labels are missing, set `batch = NULL` and provide `run_order` so WiNN
can segment the run automatically.

```r
corrected_qc_agnostic <- winn(
  observed_intensity,
  batch = NULL,
  run_order = run_order,
  lag = NULL,
  pelt_penalty = "mbic"
)
```

## Practical Notes

- Use `lag = NULL` unless you have a study-specific reason to force a fixed
  Ljung-Box lag.
- Use `remove_batch_effects = "anova"` for a lightweight mean-only correction;
  switch to `remove_batch_effects = "combat"` when you want empirical Bayes
  adjustment and `sva` is available.
- Keep `scale_by_batch = FALSE` when absolute abundance differences matter for
  downstream interpretation.
- The vignette provides a larger end-to-end example:
  `vignette("winn_tutorial", package = "winn")`.

## Citation

If you use `winn`, cite the WiNN method paper:

Demler O, Giulianini F, MacFarlane C, Tanna T, and collaborators (2024).
"WiNNbeta: Batch and drift correction method by white noise normalization for
metabolomic studies." arXiv preprint, arXiv:2404.07906.

