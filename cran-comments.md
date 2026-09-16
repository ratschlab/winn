# CRAN resubmission: winn 0.1.5

This is a resubmission following the incoming checks for winn 0.1.4 dated
2026-09-16. Thank you for reviewing the package.

## Changes addressing the incoming checks

1. The Debian check reported 5.143 seconds for the ComBat example. Loading
   the optional Bioconductor package sva accounts for most of the runtime.
   The example is now wrapped in `\donttest{}` with an explanatory comment,
   following the CRAN Cookbook guidance for examples exceeding five seconds.
   It remains executable with `example()` and is exercised by the additional
   `--run-donttest` pass of `R CMD check --as-cran`. A small ComBat integration
   test also exercises batch correction in the routine test suite.

2. The Title and Description fields have been revised. The title is now
   "White Noise Normalization for Mass Spectrometry Profiling Data".
   "Spectrometry" is the correct scientific term for the analytical technique
   and has been retained in the title. If flagged, it is a spelling false
   positive. All other previously flagged words have been removed from
   DESCRIPTION. The complete, accurate method citation is retained in
   inst/CITATION and the README.

The "New submission" NOTE is expected because winn is not yet on CRAN.
The version has been incremented from 0.1.4 to 0.1.5 for this resubmission.

Reference for the example change:
<https://contributor.r-project.org/cran-cookbook/general_issues.html#structuring-of-examples>

## Validation

The source tarball was produced with `R CMD build .` and checked with
`R CMD check --as-cran --timings winn_0.1.5.tar.gz` on macOS arm64, R 4.6.1,
using HTML Tidy 5.8.0.

There were 0 ERRORs, 0 WARNINGs, and 1 NOTE: "New submission". The normal
examples and the additional examples with `--run-donttest` both passed.
All 95 test expectations passed, with no failures, warnings, or skips.
Vignette rebuilding and PDF and HTML manual checks passed.

The CRAN additional incoming spelling check is not enabled in the local
environment. The retained scientific term "Spectrometry" is explained above;
all other previously flagged words have been removed from DESCRIPTION.
