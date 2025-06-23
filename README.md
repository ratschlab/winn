# winn (White Noise Normalization) package for pre-processing metabolomics data


Winn provides functions for metabolomics data correction using a multi-step pipeline:
1. **Outlier Adjustment** using the MAD method.
2. **Drift Correction** using autocorrelation-based detrending.
3. **Batch Effect Correction** using ANOVA-based correction or ComBat.
4. **Median Adjustment** by dilution factor-based normalization (PQN) or shrinking.
5. **Optional Scaling** by batch.

If no batch information is supplied, the package auto-detects segments (as batchs) using an fkPELT-based approach.

## Installation

To install the package and all dependencies (including Bioconductor packages):

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor dependencies
BiocManager::install("sva")

# Install winn from GitHub
devtools::install_github("ratschlab/winn")


## Authors

* [**Tanmay Tanna**](https://github.com/TanmayTanna)

