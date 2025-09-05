from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
from .utils import as_matrix

def adjust_outliers_mad(data: pd.DataFrame) -> pd.DataFrame:
    """
    Threshold is med ± 4*MAD (MAD×1.4826); out-of-bounds values are linearly "pulled back" to med ± 3*MAD.
    Rows = metabolites, columns = samples.
    """
    X = as_matrix(data).astype(float)
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

        # Reference point
        try:
            ref_hi = np.nanmax(x[~mask_hi & (x < upper)])
        except ValueError:
            ref_hi = upper
        try:
            ref_lo = np.nanmin(x[~mask_lo & (x > lower)])
        except ValueError:
            ref_lo = lower

        # Pull back upper side
        if np.any(mask_hi):
            x_hi = x[mask_hi]
            max_out = np.nanmax(x_hi)
            if np.isfinite(ref_hi) and (max_out - ref_hi) > 1e-12:
                y_hi = hi_tgt + (x_hi - ref_hi) * (upper - hi_tgt) / (max_out - ref_hi)
                x[mask_hi] = y_hi
            else:
                x[mask_hi] = upper

        # Pull back lower side
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
