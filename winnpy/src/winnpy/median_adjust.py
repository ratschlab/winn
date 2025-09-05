from __future__ import annotations
import numpy as np
import pandas as pd
from .utils import as_matrix

def normalize_by_dilution_factor(
    data: pd.DataFrame,
    processing: str = "shrink",
    control_samples = None,
) -> pd.DataFrame:
    X = as_matrix(data)
    if control_samples is not None:
        ref = X.loc[:, control_samples].median(axis=1, skipna=True)
    else:
        ref = X.median(axis=1, skipna=True)

    quot = X.divide(ref.replace(0, np.nan), axis=0)
    dilution = quot.median(axis=0, skipna=True)

    if processing == "normalize":
        return X.divide(dilution.replace(0, np.nan), axis=1)

    if processing == "shrink":
        mu = dilution.mean()
        sd = dilution.std(ddof=1)
        low, high = mu - sd, mu + sd
        scale = dilution.clip(lower=low, upper=high).replace(0, np.nan)
        return X.divide(scale, axis=1)

    raise ValueError("processing must be 'shrink' or 'normalize'")
