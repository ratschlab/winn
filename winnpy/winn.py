# src/winnpy/winn.py
from __future__ import annotations
from typing import Optional, Sequence, Union
import numpy as np
import pandas as pd

from scipy.stats import median_abs_deviation
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.api import OLS, add_constant

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
    """
    processing: "shrink" | "normalize"
      - normalize: PQN style, scale by the sample's "median of quotients"
      - shrink: If the sample's dilution factor deviates from the mean by >1 SD, shrink to (mean±sd)
    """
    X = _as_matrix(data)
    if control_samples is not None:
        cols = control_samples
        ref = X.loc[:, cols].median(axis=1, skipna=True)  # Row median: reference spectrum
    else:
        ref = X.median(axis=1, skipna=True)

    # quotients: for each sample, this row / reference spectrum
    quot = X.divide(ref.replace(0, np.nan), axis=0)
    dilution = quot.median(axis=0, skipna=True)  # Median of quotients for each column (sample)

    if processing == "normalize":
        scale = dilution.replace(0, np.nan)
        Y = X.divide(scale, axis=1)
        return Y

    if processing == "shrink":
        mu = dilution.mean()
        sd = dilution.std(ddof=1)
        low, high = mu - sd, mu + sd
        scale = dilution.copy()
        # Only shrink when outside 1SD range
        scale = scale.where(~(scale < low), other=low)
        scale = scale.where(~(scale > high), other=high)
        scale = scale.replace(0, np.nan)
        Y = X.divide(scale, axis=1)
        return Y

    raise ValueError("processing must be 'shrink' or 'normalize'")


# ---------------------------
# adjust_outliers_mad
# ---------------------------
def adjust_outliers_mad(data: Union[pd.DataFrame, np.ndarray]) -> pd.DataFrame:
    """
    For each metabolite (row), adjust outliers using MAD.
    Rule matches R: threshold is med ± 4*MAD; out-of-bounds values are linearly "pulled back" to med ± 3*MAD.
    Note: MAD is multiplied by 1.4826 as in R by default.
    """
    X = _as_matrix(data).astype(float)
    Y = X.copy()
    for i, row in X.iterrows():
        x = row.values.astype(float)
        med = np.nanmedian(x)
        mad = median_abs_deviation(x, nan_policy="omit", scale="normal")  # ×1.4826
        if not np.isfinite(mad) or mad == 0:
            continue
        lower, upper = med - 4 * mad, med + 4 * mad
        lo_tgt, hi_tgt = med - 3 * mad, med + 3 * mad

        mask_hi = x >= upper
        mask_lo = x <= lower
        # Reference point: non-outlier and close to the threshold on each side
        try:
            ref_hi = np.nanmax(x[~mask_hi & (x < upper)])
        except ValueError:
            ref_hi = upper
        try:
            ref_lo = np.nanmin(x[~mask_lo & (x > lower)])
        except ValueError:
            ref_lo = lower

        # Linearly pull back to [lo/hi_tgt, lower/upper]
        # Upper side: map [ref_hi, max_out] to [hi_tgt, upper]
        if np.any(mask_hi):
            x_hi = x[mask_hi]
            max_out = np.nanmax(x_hi)
            if np.isfinite(ref_hi) and (max_out - ref_hi) > 1e-12:
                y_hi = hi_tgt + (x_hi - ref_hi) * (upper - hi_tgt) / (max_out - ref_hi)
                x[mask_hi] = y_hi
            else:
                x[mask_hi] = upper

        # Lower side: map [min_out, ref_lo] to [lower, lo_tgt]
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
) -> pd.DataFrame:
    """
    For each batch and each metabolite: first test for autocorrelation (LB or DW), if significant, use spline/LOESS to detrend; otherwise, use detrend rule (mean/keep).
    Note: DW p-value is not officially implemented in Python, only Ljung-Box is supported here (DW branch raises NotImplementedError).
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

    Y = X.copy()
    for b in np.unique(batch):
        idx = np.where(batch == b)[0]
        seg = X.iloc[:, idx]
        seg_run = run_order[idx]

        # Compute p-values for each row
        pvals = []
        for i, row in seg.iterrows():
            x = row.values.astype(float)
            if test.lower().startswith("ljung"):
                # LB test on original sequence
                # Only take the last lag's p-value, similar to R Box.test(..., lag=lag)
                try:
                    lb = acorr_ljungbox(x, lags=[lag], return_df=True)
                    p = float(lb["lb_pvalue"].iloc[0])
                except Exception:
                    p = np.nan
            else:
                raise NotImplementedError("DW test p-value is not implemented in Python skeleton.")
            pvals.append(p)
        pvals = np.array(pvals, dtype=float)
        # FDR (BH)
        ok = np.isfinite(pvals)
        p_adj = np.ones_like(pvals)
        if ok.any():
            p_adj[ok] = multipletests(pvals[ok], alpha=fdr_threshold, method="fdr_bh")[1]
        correct_ids = np.where(p_adj < fdr_threshold)[0]

        # Detrend
        if detrend == "spline":
            # Detrend all rows using cubic polynomial approximation (same as R)
            t = (seg_run - seg_run.min()) / max(1, (seg_run.max() - seg_run.min()))
            T = np.vstack([np.ones_like(t), t, t**2, t**3]).T  # cubic polynomial approximation
            for ridx in range(seg.shape[0]):  # ← Key: process all rows
                y = seg.iloc[ridx, :].values
                beta = np.linalg.lstsq(T, np.log1p(y), rcond=None)[0]
                fit = T @ beta
                fit_centered = fit - fit.mean()
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)

        elif detrend == "mean":
            # Only process rows with significant autocorrelation; keep others unchanged (same as R)
            t = (seg_run - seg_run.min()) / max(1, (seg_run.max() - seg_run.min()))
            T = np.vstack([np.ones_like(t), t, t**2, t**3]).T
            for ridx in correct_ids:
                y = seg.iloc[ridx, :].values
                beta = np.linalg.lstsq(T, np.log1p(y), rcond=None)[0]
                fit = T @ beta
                fit_centered = fit - fit.mean()
                Y.iloc[ridx, idx] = np.expm1(np.log1p(y) - fit_centered)
        else:
            raise ValueError("detrend must be 'mean' or 'spline'")

    return Y


# ---------------------------
# anova_batch_correction (mean-only)
# ---------------------------
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
    for i, row in Z.iterrows():
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
    """
    Call pyComBat (raise error if not installed). Keep processing in log1p space, revert with expm1.
    """
    if _pycombat is None:
        raise ImportError("pyComBat not found. Install `combat` (pip install combat) or use ANOVA.")
    X = _as_matrix(data).astype(float)
    batch = pd.Categorical(batch)
    Z = np.log1p(X)
    # pycombat(dat: genes×samples, batch: samples)
    dat = Z.values
    corrected = _pycombat(dat=dat, batch=batch.codes, par_prior=par_prior, mean_only=mean_only, ref_batch=ref_batch)
    return pd.DataFrame(np.expm1(corrected), index=X.index, columns=X.columns)


# ---------------------------
# scale_by_batch (z-score per batch per metabolite)
# ---------------------------
def scale_by_batch_func(
    data: Union[pd.DataFrame, np.ndarray],
    batch: ArrayLike,
) -> pd.DataFrame:
    """
    As in R: for each batch and each metabolite, (x - batch mean)/batch std; if sd=0, use 1 instead.
    """
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
# Helpers: batch auto-detect, control corr
# ---------------------------
def _auto_detect_batch(data: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
    """
    Use ruptures' PELT to detect change points on the sequence of "median total intensity" per sample, approximating R's fkPELT.
    Penalty is set to 3*log(300) (corresponds to R constant); if ruptures is not installed, fallback to single batch.
    """
    X = _as_matrix(data)
    agg = X.median(axis=0, skipna=True).values  # Median of each column
    n = len(agg)
    if not _HAS_RUPTURES or n < 5:
        return np.ones(n, dtype=int)
    algo = rpt.Pelt(model="rbf").fit(agg)
    # Penalty can be further tuned to match R's number of change points
    pen = 3.0 * np.log(300.0)
    chg = algo.predict(pen=pen)  # The last n is included in the boundaries
    # chg is the right endpoint index for each segment, e.g. [k1, k2, ..., n]
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
    corr = ctrl.corr()  # Pearson
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
) -> pd.DataFrame:
    """
    Python reimplementation (skeleton) of the Winn pipeline.
    """
    X0 = _as_matrix(data).astype(float)

    # 1) MAD outlier
    X1 = adjust_outliers_mad(X0)

    # Automatic parameter search (only if control_samples provided and parameters="auto")
    if control_samples is not None and parameters == "auto":
        # Simplified grid, can be expanded as needed
        tests = ["Ljung-Box"]  # DW p-value not supported in Python
        norms = ["shrink", "normalize"]
        acorr_fdrs = [0.1, 0.05, 0.01]
        anova_fdrs = [0.1, 0.05, 0.01]
        scales = [True, False]
        best_score = -np.inf
        best = None

        batch_options = ["provided"] if batch is not None else ["auto"]

        for bopt in batch_options:
            current_batch = np.asarray(batch) if bopt == "provided" else _auto_detect_batch(X1)
            for tname in tests:
                for acfdr in acorr_fdrs:
                    X2 = autocorrelation_correct(
                        X1, run_order=run_order, batch=current_batch, lag=lag,
                        test=tname, detrend=detrend_non_autocorrelated, fdr_threshold=acfdr,
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
                                # Scoring: mean correlation between control samples - SDR
                                ctrl = X5.loc[:, control_samples]
                                with np.errstate(divide="ignore", invalid="ignore"):
                                    sdr = (ctrl.std(axis=1) / ctrl.mean(axis=1)).replace([np.inf, -np.inf], np.nan).mean()
                                mean_corr = _mean_control_correlation(X5, control_samples)
                                score = (mean_corr if np.isfinite(mean_corr) else -np.inf) - (sdr if np.isfinite(sdr) else 0.0)
                                if score > best_score:
                                    best_score = score
                                    best = (X5, dict(batch_option=bopt, test=tname, acorr_fdr=acfdr,
                                                     anova_fdr=afdr, median=nm, scale=sc))
        if best is None:
            raise RuntimeError("Auto-parameter search failed to produce a result.")
        return best[0]

    # Fixed parameters
    current_batch = np.asarray(batch) if batch is not None else _auto_detect_batch(X1)

    X2 = autocorrelation_correct(
        X1, run_order=run_order, batch=current_batch, lag=lag,
        test=test, detrend=detrend_non_autocorrelated, fdr_threshold=fdr_threshold,
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
