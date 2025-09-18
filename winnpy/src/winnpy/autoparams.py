from __future__ import annotations
import numpy as np
import pandas as pd
from .outlier import adjust_outliers_mad
from .drift import autocorrelation_correct
from .batch import anova_batch_correction, combat_batch_correction
from .median_adjust import normalize_by_dilution_factor
from .scale import scale_by_batch_func
from .utils import auto_detect_batch, mean_control_correlation

def auto_select(
    data: pd.DataFrame,
    batch,
    run_order,
    control_samples,
    detrend_non_autocorrelated: str = "mean",
    remove_batch_effects: str = "anova",
    lag: int = 20,
    # New: can pass a fixed spline_method; if you want to grid search, pass None
    spline_method: str | None = None,
):
    """
    Returns (best_data, best_params); grid covers the same as the R version:
    - median/shrink
    - acorr FDR ∈ {0.1,0.05,0.01}
    - anova FDR ∈ {0.1,0.05,0.01}
    - scale ∈ {True, False}
    - spline_method ∈ {"conservative","standard"} (if not specified)
    """
    X1 = adjust_outliers_mad(data)
    batch_options = ["provided"] if batch is not None else ["auto"]
    tests = ["Ljung-Box"]
    norms = ["shrink", "normalize"]
    acorr_fdrs = [0.1, 0.05, 0.01]
    anova_fdrs = [0.1, 0.05, 0.01]
    scales = [True, False]
    spline_methods = [spline_method] if spline_method in ("conservative", "standard") else ["conservative", "standard"]

    best_score = -np.inf
    best = None
    for bopt in batch_options:
        current_batch = np.asarray(batch) if bopt == "provided" else auto_detect_batch(X1)
        for tname in tests:
            for spm in spline_methods:
                for acfdr in acorr_fdrs:
                    X2 = autocorrelation_correct(
                        X1, run_order, current_batch, lag, tname,
                        detrend_non_autocorrelated, acfdr, spline_method=spm
                    )
                    for afdr in anova_fdrs:
                        if remove_batch_effects == "anova":
                            X3 = anova_batch_correction(X2, current_batch, fdr_threshold=afdr)
                        else:
                            X3 = combat_batch_correction(X2, current_batch)
                        for nm in norms:
                            X4 = normalize_by_dilution_factor(X3, processing=nm, control_samples=control_samples)
                            for sc in scales:
                                X5 = scale_by_batch_func(X4, current_batch) if sc else X4
                                ctrl = X5.loc[:, control_samples]
                                with np.errstate(divide="ignore", invalid="ignore"):
                                    sdr = (ctrl.std(axis=1)/ctrl.mean(axis=1)).replace([np.inf, -np.inf], np.nan).mean()
                                mean_corr = mean_control_correlation(X5, control_samples)
                                score = (mean_corr if np.isfinite(mean_corr) else -np.inf) - (sdr if np.isfinite(sdr) else 0.0)
                                if score > best_score:
                                    best_score = score
                                    best = (X5, dict(batch_option=bopt, test=tname, acorr_fdr=acfdr,
                                                     anova_fdr=afdr, median=nm, scale=sc, spline_method=spm))
    if best is None:
        raise RuntimeError("Auto-parameter search failed.")
    return best
