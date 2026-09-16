# CRAN submission comments

## Test environments

* Local macOS Tahoe 26.6.2, R 4.6.1 (aarch64):
  `R CMD check --as-cran --timings winn_0.1.4.tar.gz`
* HTML manual validation used HTML Tidy 5.8.0.

## R CMD check results

There were 0 ERRORs, 0 WARNINGs, and 1 NOTE: `New submission`.

All 90 test expectations passed without failures, warnings, or skips.
Examples, vignette rebuilding, and PDF and HTML manual checks passed.

This is the first CRAN submission of winn, version 0.1.4. The optional `sva`
dependency is a Bioconductor package used only for the ComBat correction path;
the default ANOVA path does not require it.

Tanmay Tanna is the sole package author and maintainer, with contact address
tanmay@tanmaytanna.com. The bibliography credits the full author list of the
underlying method paper separately from package authorship.
