# CRAN submission comments

## Test environments

* Local macOS Tahoe 26.6.2, R 4.6.1 (aarch64): `R CMD check --as-cran`

## R CMD check results

There were no ERRORs or WARNINGs.

The local check reports the expected `New submission` note. The optional `sva`
dependency is a Bioconductor package used only for the ComBat correction path;
the default ANOVA path does not require it. A full local check may additionally
report that HTML manual validation was skipped when the installed HTML Tidy
version is too old.
