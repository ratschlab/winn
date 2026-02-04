# winn — White Noise Normalization for Metabolomics

**winn** is an R package for correcting technical artifacts in LC‑MS/GC‑MS metabolomics. It is designed to be practical, transparent, and robust: each correction step is explicit and can be switched on/off. The core function, `winn()`, applies a multi‑step pipeline to reduce **drift**, **batch effects**, **outliers**, and **sample‑to‑sample dilution**, while preserving biological signal.

This package was built for real-world workflows where run order, batch structure, and QC samples are available but noisy. It works well on both targeted and untargeted datasets.

---

## What WINN Does (High‑Level)

For a data matrix of **metabolites × samples**:

1. **Outlier adjustment** using a robust MAD rule (per metabolite).
2. **Drift correction** using autocorrelation‑guided detrending by run order (per batch/segment).
3. **Batch correction** with either ANOVA mean‑shift or ComBat (log1p space).
4. **Median/PQN adjustment** to correct dilution or global scaling.
5. **Optional per‑batch scaling** (z‑score within batch).

If batch labels are not provided, WINN can **auto‑detect batch segments** using an fkPELT change‑point approach on the median sample signal.

---

## Installation

```r
# CRAN dependencies
install.packages(c("mgcv", "lmtest"))

# Bioconductor dependency for ComBat
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("sva")

# Install winn from GitHub
# devtools::install_github("ratschlab/winn")
```

If you do not plan to use ComBat, you can skip the `sva` installation and use `remove_batch_effects = "anova"`.

---

## Quick Start

```r
library(winn)

# data: metabolites x samples
# batch: length = ncol(data)
# run_order: length = ncol(data)

corrected <- winn(
  data,
  batch = batch,
  run_order = run_order,
  median_adjustment = "shrink",
  remove_batch_effects = "anova"
)
```

### With QC-based auto parameter selection
If you have pooled QC samples, WINN can tune parameters automatically to maximize QC consistency and reduce variability.

```r
qc_idx <- c(5, 15, 25, 35)  # column indices of QC samples

corrected_auto <- winn(
  data,
  batch = batch,
  run_order = run_order,
  control_samples = qc_idx,
  parameters = "auto"
)
```

---

## Recommended Workflow

1. **Provide batch + run order** if known.
2. **Use QC samples** and `parameters = "auto"` whenever possible.
3. Inspect QC consistency before/after (CV, correlation).
4. If results appear over‑corrected, relax `fdr_threshold` or disable scaling.

---

## Key Parameters (Practical Guide)

- `remove_batch_effects = "anova"` (default): fast and conservative mean‑shift correction.
- `remove_batch_effects = "combat"`: stronger correction, needs `sva`; best when batch effects dominate.
- `median_adjustment = "shrink"`: mild dilution correction; good default.
- `median_adjustment = "normalize"`: full PQN; can be more aggressive.
- `scale_by_batch = TRUE`: z‑scores within batch; use cautiously if absolute levels matter.
- `pelt_penalty`: optional override for auto batch segmentation (higher = fewer segments).

---

## Interpretation and Limitations

WINN is meant to correct **technical artifacts** while preserving biological variation. However:

- If run order correlates with phenotype, drift correction can remove true biology.
- If batch is confounded with biological groups, ComBat or ANOVA can dampen real effects.
- Over‑aggressive scaling may erase absolute concentration information.

For these reasons, **diagnostic checks (QC CV, correlation, PCA)** should be part of every analysis.

---

## Tutorial

A full tutorial with simulated drift + batch effects is available in the vignette:

- `vignettes/winn_tutorial.Rmd`

It demonstrates QC‑based auto parameter selection and provides quantitative metrics for signal preservation.

---

## Citation

If you use WINN in a manuscript, please cite the GitHub repository and include the package version used. A formal paper citation can be added once available.

---

## License

GPL‑3

---

## Contact

Issues and suggestions are welcome via GitHub.
