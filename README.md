# Winn package


Winn provides functions for metabolomics data correction using a multi-step pipeline:
1. **Median Adjustment** by dilution factor.
2. **Outlier Adjustment** using the MAD method.
3. **Drift Correction** using autocorrelation-based detrending.
4. **Batch Effect Correction** using ANOVA-based residualization.
5. **Optional Scaling** by plate.

If no plate information is supplied, the package auto-detects segments (as plates) using an fkPELT-based approach.

## Authors

* [**Tanmay Tanna**](https://github.com/TanmayTanna)

