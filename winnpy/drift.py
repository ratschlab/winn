from __future__ import annotations
import numpy as np
import pandas as pd
from statsmodels.stats.diagnostic import acorr_ljungbox
from .utils import as_matrix, fdr_bh

def _cubic_design(t: np.ndarray) -> np.ndarray:
    t = (t - t.min()) / max(1, (t.max() - t.min()))
    return np.vstack([np.ones_like(t), t, t**2, t**3]).T

def autocorrelation_correct(
    data: pd.DataFrame,
    run_order: np.ndarray | None,
    batch: np.ndarray,
    lag: int = 20,
    test: str = "Ljung-Box",
    detrend: str = "mean",
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    X = as_matrix(data).astype(float)
    if batch is None:
        raise ValueError("batch must be provided.")
    batch = np.asarray(batch)
    if X.shape[1] != len(batch):
        raise ValueError("Length of batch must equal number of columns in data.")

    if run_order is None:
        run_order = np.arange(X.shape[1])
    run_order = np.asarray(run_order)

    Y = X.copy()
    for b in np.unique(batch):
        idx = np.where(batch == b)[0]
        seg = X.iloc[:, idx]
        seg_run = run_order[idx]

        # Compute p-values (only supports Ljung-Box)
        pvals = []
        for _, row in seg.iterrows():
            x = row.values.astype(float)
            if test.lower().startswith("ljung"):
                try:
                    lb = acorr_ljungbox(x, lags=[lag], return_df=True)
                    p = float(lb["lb_pvalue"].iloc[0])
                except Exception:
                    p = np.nan
            else:
                raise NotImplementedError("DW p-value not implemented.")
            pvals.append(p)
        p_adj = fdr_bh(pvals, alpha=fdr_threshold)
        correct_ids = np.where(p_adj < fdr_threshold)[0]

        T = _cubic_design(seg_run)

        if detrend == "spline":
            # Detrend all rows using spline approximation (same as R)
            for ridx in range(seg.shape[0]):
                y = seg.iloc[ridx, :].values
                beta = np.linalg.lstsq(T, np.log1p(y), rcond=None)[0]
                fit = T @ beta
                fit_centered = fit - fit.mean()
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)

        elif detrend == "mean":
            # Only detrend rows with significant autocorrelation, keep others unchanged
            for ridx in correct_ids:
                y = seg.iloc[ridx, :].values
                beta = np.linalg.lstsq(T, np.log1p(y), rcond=None)[0]
                fit = T @ beta
                fit_centered = fit - fit.mean()
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)
        else:
            raise ValueError("detrend must be 'mean' or 'spline'")

    return Y
