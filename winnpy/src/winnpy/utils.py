from __future__ import annotations
from typing import Sequence, Union, Optional
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# optional: ruptures for PELT
try:
    import ruptures as rpt
    _HAS_RUPTURES = True
except Exception:
    _HAS_RUPTURES = False

ArrayLike = Union[pd.Series, np.ndarray, list]

def as_matrix(data: Union[pd.DataFrame, np.ndarray]) -> pd.DataFrame:
    if isinstance(data, np.ndarray):
        return pd.DataFrame(data)
    if isinstance(data, pd.DataFrame):
        return data.copy()
    raise TypeError("data must be a pandas.DataFrame or numpy.ndarray")

def fdr_bh(pvals: Sequence[float], alpha: float) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    ok = np.isfinite(pvals)
    padj = np.ones_like(pvals)
    if ok.any():
        padj[ok] = multipletests(pvals[ok], alpha=alpha, method="fdr_bh")[1]
    return padj

def row_cv(df: pd.DataFrame) -> pd.Series:
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = sd / mu
        cv.replace([np.inf, -np.inf], np.nan, inplace=True)
    return cv

def mean_control_correlation(data: pd.DataFrame, control_cols: Sequence[Union[int, str]]) -> float:
    ctrl = data.loc[:, control_cols]
    if ctrl.shape[1] < 2:
        return np.nan
    corr = ctrl.corr()
    tril = corr.values[np.tril_indices_from(corr.values, k=-1)]
    return float(np.nanmean(tril))

def auto_detect_batch(
    data: Union[pd.DataFrame, np.ndarray],
    model: str = "rbf",
    pen: Optional[float] = None
) -> np.ndarray:
    """
    PELT approximates fkPELT: finds change points on the sequence of column medians; falls back to single batch if ruptures is not available.
    """
    X = as_matrix(data)
    agg = X.median(axis=0, skipna=True).values
    n = len(agg)
    if not _HAS_RUPTURES or n < 5:
        return np.ones(n, dtype=int)
    if pen is None:
        pen = 3.0 * np.log(300.0)
    algo = rpt.Pelt(model=model).fit(agg)
    chg = algo.predict(pen=pen)  # right endpoints (includes n)
    labels = np.empty(n, dtype=int)
    start = 0
    for seg_id, end in enumerate(chg, start=1):
        labels[start:end] = seg_id
        start = end
    return labels
