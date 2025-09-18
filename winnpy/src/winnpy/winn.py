from __future__ import annotations
from typing import Optional, Sequence, Union
import numpy as np
import pandas as pd

from scipy.stats import median_abs_deviation
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import durbin_watson
from scipy import stats
from scipy.interpolate import UnivariateSpline

# Optional dependency: ruptures (batch auto-segmentation)
try:
    import ruptures as rpt
    _HAS_RUPTURES = True
except Exception:
    _HAS_RUPTURES = False

# Optional dependency: ComBat
_pycombat = None
try:
    # pip install combat
    from combat.pycombat import pycombat as _pycombat
except Exception:
    pass

ArrayLike = Union[pd.Series, np.ndarray, list]
Number = Union[int, float]

def _as_matrix(data: Union[pd.DataFrame, np.ndarray]) -> pd.DataFrame:
    if isinstance(data, np.ndarray):
        return pd.DataFrame(data)
    if isinstance(data, pd.DataFrame):
        return data.copy()
    raise TypeError("data must be a pandas.DataFrame or numpy.ndarray")

# ---------------------------
# normalize_by_dilution_factor
# ---------------------------
def normalize_by_dilution_factor(
    data: Union[pd.DataFrame, np.ndarray],
    processing: str = "shrink",
    control_samples: Optional[Sequence[Union[int, str]]] = None,
) -> pd.DataFrame:
    X = _as_matrix(data)
    if control_samples is not None:
        ref = X.loc[:, control_samples].median(axis=1, skipna=True)
    else:
        ref = X.median(axis=1, skipna=True)
    quot = X.divide(ref.replace(0, np.nan), axis=0)
    dilution = quot.median(axis=0, skipna=True)
    if processing == "normalize":
        scale = dilution.replace(0, np.nan)
        return X.divide(scale, axis=1)
    if processing == "shrink":
        mu = dilution.mean()
        sd = dilution.std(ddof=1)
        low, high = mu - sd, mu + sd
        scale = dilution.clip(lower=low, upper=high).replace(0, np.nan)
        return X.divide(scale, axis=1)
    raise ValueError("processing must be 'shrink' or 'normalize'")

# ---------------------------
# adjust_outliers_mad
# ---------------------------
def adjust_outliers_mad(data: Union[pd.DataFrame, np.ndarray]) -> pd.DataFrame:
    X = _as_matrix(data).astype(float)
    Y = X.copy()
    for i, row in X.iterrows():
        x = row.values.astype(float)
        med = np.nanmedian(x)
        mad = median_abs_deviation(x, nan_policy="omit", scale="normal")
        if not np.isfinite(mad) or mad == 0:
            continue
        lower, upper = med - 4 * mad, med + 4 * mad
        lo_tgt, hi_tgt = med - 3 * mad, med + 3 * mad
        mask_hi = x >= upper
        mask_lo = x <= lower
        try:
            ref_hi = np.nanmax(x[~mask_hi & (x < upper)])
        except ValueError:
            ref_hi = upper
        try:
            ref_lo = np.nanmin(x[~mask_lo & (x > lower)])
        except ValueError:
            ref_lo = lower
        if np.any(mask_hi):
            x_hi = x[mask_hi]
            max_out = np.nanmax(x_hi)
            if np.isfinite(ref_hi) and (max_out - ref_hi) > 1e-12:
                y_hi = hi_tgt + (x_hi - ref_hi) * (upper - hi_tgt) / (max_out - ref_hi)
                x[mask_hi] = y_hi
            else:
                x[mask_hi] = upper
        if np.any(mask_lo):
            x_lo = x[mask_lo]
            min_out = np.nanmin(x_lo)
            if np.isfinite(ref_lo) and (ref_lo - min_out) > 1e-12:
                y_lo = lower + (x_lo - min_out) * (lo_tgt - lower) / (ref_lo - min_out)
                x[mask_lo] = y_lo
            else:
                x[mask_lo] = lower
        Y.loc[i, :] = x
    return Y

# ---------------------------
# helpers for drift
# ---------------------------
def _cubic_design(t: np.ndarray) -> np.ndarray:
    t = (t - t.min()) / max(1, (t.max() - t.min()))
    return np.vstack([np.ones_like(t), t, t**2, t**3]).T

def _fit_trend(seg_run: np.ndarray, y: np.ndarray, method: str) -> np.ndarray:
    t = (seg_run - seg_run.min()) / max(1, (seg_run.max() - seg_run.min()))
    ly = np.log1p(y)
    if method == "standard":
        try:
            s = max(1e-8, 0.2 * len(t))
            spl = UnivariateSpline(t, ly, k=3, s=s)
            fit = spl(t)
            return fit - fit.mean()
        except Exception:
            pass
    T = _cubic_design(seg_run)
    beta = np.linalg.lstsq(T, ly, rcond=None)[0]
    fit = T @ beta
    return fit - fit.mean()

def _dw_pvalue_approx(x: np.ndarray) -> float:
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

# ---------------------------
# autocorrelation_correct
# ---------------------------
def autocorrelation_correct(
    data: Union[pd.DataFrame, np.ndarray],
    run_order: Optional[ArrayLike] = None,
    batch: ArrayLike = None,
    lag: int = 20,
    test: str = "Ljung-Box",
    detrend: str = "mean",
    fdr_threshold: float = 0.05,
    spline_method: str = "conservative",  # New: pass through
) -> pd.DataFrame:
    """
    test: "Ljung-Box" or "DW"
    detrend: "mean" | "spline"
    spline_method: "conservative" (cubic polynomial) | "standard" (prefer spline, fallback if fails)
    """
    X = _as_matrix(data).astype(float)
    if batch is None:
        raise ValueError("batch must be provided (vector length = number of columns).")
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

        pvals = np.array(pvals, dtype=float)
        ok = np.isfinite(pvals)
        p_adj = np.ones_like(pvals)
        if ok.any():
            p_adj[ok] = multipletests(pvals[ok], alpha=fdr_threshold, method="fdr_bh")[1]
        correct_ids = np.where(p_adj < fdr_threshold)[0]

        if detrend_lower == "spline":
            for ridx in range(seg.shape[0]):
                y = seg.iloc[ridx, :].values
                fit_centered = _fit_trend(seg_run, y, spline_method)
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)
        elif detrend_lower == "mean":
            for ridx in correct_ids:
                y = seg.iloc[ridx, :].values
                fit_centered = _fit_trend(seg_run, y, spline_method)
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)
        else:
            raise ValueError("detrend must be 'mean' or 'spline'")

    return Y

# ---------------------------
# anova_batch_correction (mean-only)
# ---------------------------
from statsmodels.api import OLS, add_constant
def anova_batch_correction(
    data: Union[pd.DataFrame, np.ndarray],
    batch: ArrayLike,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    X = _as_matrix(data).astype(float)
    batch = pd.Categorical(batch)
    if X.shape[1] != len(batch):
        raise ValueError("Length of batch must equal number of columns in data.")
    Z = np.log1p(X)
    B = pd.get_dummies(batch, drop_first=True).astype(float)
    Xdesign = add_constant(B.values, has_constant="add")
    pvals = []
    for _, row in Z.iterrows():
        y = row.values
        model = OLS(y, Xdesign).fit()
        pvals.append(float(model.f_pvalue) if np.isfinite(model.f_pvalue) else 1.0)
    pvals = np.array(pvals, dtype=float)
    ok = np.isfinite(pvals)
    padj = np.ones_like(pvals)
    if ok.any():
        padj[ok] = multipletests(pvals[ok], alpha=fdr_threshold, method="fdr_bh")[1]
    Y = Z.copy()
    overall_means = Z.mean(axis=1)
    df_batch = pd.Series(batch).astype(str).values
    sig_idx = np.where(padj < fdr_threshold)[0]
    if len(sig_idx) > 0:
        for ridx in sig_idx:
            vals = Z.iloc[ridx, :].values
            batch_means = pd.Series(vals).groupby(df_batch).mean()
            shift = batch_means - overall_means.iloc[ridx]
            Y.iloc[ridx, :] = vals - pd.Index(df_batch).map(shift).values
    return pd.DataFrame(np.expm1(Y.values), index=X.index, columns=X.columns)

# ---------------------------
# ComBat
# ---------------------------
def combat_batch_correction(
    data: Union[pd.DataFrame, np.ndarray],
    batch: ArrayLike,
    par_prior: bool = True,
    mean_only: bool = False,
    ref_batch: Optional[Union[int, str]] = None,
) -> pd.DataFrame:
    if _pycombat is None:
        raise ImportError("pyComBat not found. Install `combat` (pip install combat) or use ANOVA.")
    X = _as_matrix(data).astype(float)
    batch = pd.Categorical(batch)
    Z = np.log1p(X)
    dat = Z.values
    corrected = _pycombat(dat=dat, batch=batch.codes, par_prior=par_prior, mean_only=mean_only, ref_batch=ref_batch)
    return pd.DataFrame(np.expm1(corrected), index=X.index, columns=X.columns)

# ---------------------------
# scale_by_batch
# ---------------------------
def scale_by_batch_func(
    data: Union[pd.DataFrame, np.ndarray],
    batch: ArrayLike,
) -> pd.DataFrame:
    X = _as_matrix(data).astype(float)
    batch = np.asarray(batch)
    if X.shape[1] != len(batch):
        raise ValueError("Length of batch must equal number of columns in data.")
    Y = X.copy()
    for b in np.unique(batch):
        idx = np.where(batch == b)[0]
        seg = X.iloc[:, idx]
        mu = seg.mean(axis=1)
        sd = seg.std(axis=1, ddof=1).replace(0, 1.0)
        Y.iloc[:, idx] = seg.subtract(mu, axis=0).divide(sd, axis=0)
    return Y

# ---------------------------
# helpers: batch auto-detect, control corr
# ---------------------------
def _auto_detect_batch(data: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
    X = _as_matrix(data)
    agg = X.median(axis=0, skipna=True).values
    n = len(agg)
    if not _HAS_RUPTURES or n < 5:
        return np.ones(n, dtype=int)
    algo = rpt.Pelt(model="rbf").fit(agg)
    pen = 3.0 * np.log(300.0)
    chg = algo.predict(pen=pen)
    labels = np.empty(n, dtype=int)
    start = 0
    for seg_id, end in enumerate(chg, start=1):
        labels[start:end] = seg_id
        start = end
    return labels

def _mean_control_correlation(data: pd.DataFrame, control_samples: Sequence[Union[int, str]]) -> float:
    ctrl = data.loc[:, control_samples]
    if ctrl.shape[1] < 2:
        return np.nan
    corr = ctrl.corr()
    tril = corr.values[np.tril_indices_from(corr.values, k=-1)]
    return float(np.nanmean(tril))

# ---------------------------
# winn (pipeline)
# ---------------------------
def winn(
    data: Union[pd.DataFrame, np.ndarray],
    batch: Optional[ArrayLike] = None,
    run_order: Optional[ArrayLike] = None,
    control_samples: Optional[Sequence[Union[int, str]]] = None,
    parameters: str = "fixed",
    fdr_threshold: float = 0.05,
    median_adjustment: str = "shrink",
    detrend_non_autocorrelated: str = "mean",
    remove_batch_effects: str = "anova",
    test: str = "Ljung-Box",
    lag: int = 20,
    scale_by_batch: bool = False,
    spline_method: str = "conservative",   # New: pass through
) -> pd.DataFrame:
    X0 = _as_matrix(data).astype(float)
    X1 = adjust_outliers_mad(X0)

    if control_samples is not None and parameters == "auto":
        tests = ["Ljung-Box"]  # DW approximation is available, but keep consistent with R's default
        norms = ["shrink", "normalize"]
        acorr_fdrs = [0.1, 0.05, 0.01]
        anova_fdrs = [0.1, 0.05, 0.01]
        scales = [True, False]
        spline_methods = ["conservative", "standard"]

        best_score = -np.inf
        best = None
        batch_options = ["provided"] if batch is not None else ["auto"]

        for bopt in batch_options:
            current_batch = np.asarray(batch) if bopt == "provided" else _auto_detect_batch(X1)
            for tname in tests:
                for spm in spline_methods:
                    for acfdr in acorr_fdrs:
                        X2 = autocorrelation_correct(
                            X1, run_order=run_order, batch=current_batch, lag=lag,
                            test=tname, detrend=detrend_non_autocorrelated,
                            fdr_threshold=acfdr, spline_method=spm,
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
                                        sdr = (ctrl.std(axis=1) / ctrl.mean(axis=1)).replace([np.inf, -np.inf], np.nan).mean()
                                    mean_corr = _mean_control_correlation(X5, control_samples)
                                    score = (mean_corr if np.isfinite(mean_corr) else -np.inf) - (sdr if np.isfinite(sdr) else 0.0)
                                    if score > best_score:
                                        best_score = score
                                        best = (X5, dict(batch_option=bopt, test=tname, acorr_fdr=acfdr,
                                                         anova_fdr=afdr, median=nm, scale=sc, spline_method=spm))
        if best is None:
            raise RuntimeError("Auto-parameter search failed to produce a result.")
        return best[0]

    current_batch = np.asarray(batch) if batch is not None else _auto_detect_batch(X1)
    X2 = autocorrelation_correct(
        X1, run_order=run_order, batch=current_batch, lag=lag,
        test=test, detrend=detrend_non_autocorrelated, fdr_threshold=fdr_threshold,
        spline_method=spline_method,
    )
    if remove_batch_effects == "anova":
        X3 = anova_batch_correction(X2, current_batch, fdr_threshold=fdr_threshold)
    else:
        X3 = combat_batch_correction(X2, current_batch)
    if median_adjustment.lower() == "none":
        X4 = X3
    else:
        X4 = normalize_by_dilution_factor(X3, processing=median_adjustment, control_samples=control_samples)
    X5 = scale_by_batch_func(X4, current_batch) if scale_by_batch else X4
    return X5
