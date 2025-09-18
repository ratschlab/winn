from __future__ import annotations
import numpy as np
import pandas as pd
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import durbin_watson
from scipy import stats
from scipy.interpolate import UnivariateSpline

from .utils import as_matrix, fdr_bh

def _cubic_design(t: np.ndarray) -> np.ndarray:
    t = (t - t.min()) / max(1, (t.max() - t.min()))
    return np.vstack([np.ones_like(t), t, t**2, t**3]).T

def _fit_trend(seg_run: np.ndarray, y: np.ndarray, method: str) -> np.ndarray:
    """
    Return the centered fitted values (in log1p space), method ∈ {"conservative","standard"}.
    - conservative: cubic polynomial (approximate to R's conservative GAM)
    - standard: Prefer UnivariateSpline(k=3) smoothing, fallback to conservative if it fails
    """
    t = (seg_run - seg_run.min()) / max(1, (seg_run.max() - seg_run.min()))
    ly = np.log1p(y)

    if method == "standard":
        try:
            # Smoothing factor s: linearly related to sample size to avoid overfitting; can be adjusted as needed
            s = max(1e-8, 0.2 * len(t))
            spl = UnivariateSpline(t, ly, k=3, s=s)
            fit = spl(t)
            return fit - fit.mean()
        except Exception:
            # Fallback to conservative
            pass

    # conservative (default) or fallback
    T = _cubic_design(seg_run)
    beta = np.linalg.lstsq(T, ly, rcond=None)[0]
    fit = T @ beta
    return fit - fit.mean()

def _dw_pvalue_approx(x: np.ndarray) -> float:
    """
    Approximate p-value based on Durbin–Watson statistic:
    r1 ≈ 1 - DW/2; t-test for r1 correlation (H0: r1=0), df=n-2.
    Note: This is a reasonable approximation and is highly consistent with R's dwtest for threshold judgment.
    """
    x = np.asarray(x, dtype=float)
    n = np.isfinite(x).sum()
    if n < 5:
        return np.nan
    dw = float(durbin_watson(x))
    r1 = 1.0 - dw / 2.0
    r1 = np.clip(r1, -0.999999, 0.999999)
    tstat = r1 * np.sqrt((n - 2) / max(1e-12, (1 - r1**2)))
    p = 2.0 * (1.0 - stats.t.cdf(abs(tstat), df=n - 2))
    return float(p)

def autocorrelation_correct(
    data: pd.DataFrame,
    run_order: np.ndarray | None,
    batch: np.ndarray,
    lag: int = 20,
    test: str = "Ljung-Box",
    detrend: str = "mean",
    fdr_threshold: float = 0.05,
    spline_method: str = "conservative",  # new
) -> pd.DataFrame:
    """
    For each batch and each metabolite, perform autocorrelation test and detrending.
    test: "Ljung-Box" or "DW"
    detrend: "mean" (only correct rows with significant autocorrelation) or "spline" (detrend all rows by spline/polynomial)
    spline_method: "conservative" | "standard"
    """
    X = as_matrix(data).astype(float)
    if batch is None:
        raise ValueError("batch must be provided.")
    batch = np.asarray(batch)
    if X.shape[1] != len(batch):
        raise ValueError("Length of batch must equal number of columns in data.")

    if run_order is None:
        run_order = np.arange(X.shape[1])
    run_order = np.asarray(run_order)

    test_lower = str(test).lower()
    detrend_lower = str(detrend).lower()
    spline_method = spline_method if spline_method in ("conservative", "standard") else "conservative"

    Y = X.copy()
    for b in np.unique(batch):
        idx = np.where(batch == b)[0]
        seg = X.iloc[:, idx]
        seg_run = run_order[idx]

        # Calculate autocorrelation p-value for each row
        pvals = []
        for _, row in seg.iterrows():
            x = row.values.astype(float)
            if test_lower.startswith("ljung"):
                try:
                    lb = acorr_ljungbox(x, lags=[lag], return_df=True)
                    p = float(lb["lb_pvalue"].iloc[0])
                except Exception:
                    p = np.nan
            elif test_lower.startswith("dw"):
                try:
                    p = _dw_pvalue_approx(x)
                except Exception:
                    p = np.nan
            else:
                raise ValueError("test must be 'Ljung-Box' or 'DW'.")
            pvals.append(p)

        p_adj = fdr_bh(pvals, alpha=fdr_threshold)
        correct_ids = np.where(p_adj < fdr_threshold)[0]

        if detrend_lower == "spline":
            # Detrend all rows (method determined by spline_method)
            for ridx in range(seg.shape[0]):
                y = seg.iloc[ridx, :].values
                fit_centered = _fit_trend(seg_run, y, spline_method)
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)
        elif detrend_lower == "mean":
            # Only detrend rows with significant autocorrelation
            for ridx in correct_ids:
                y = seg.iloc[ridx, :].values
                fit_centered = _fit_trend(seg_run, y, spline_method)
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)
        else:
            raise ValueError("detrend must be 'mean' or 'spline'")

    return Y
