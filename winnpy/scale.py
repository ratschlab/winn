from __future__ import annotations
import numpy as np
import pandas as pd
from .utils import as_matrix

def scale_by_batch_func(data: pd.DataFrame, batch: np.ndarray) -> pd.DataFrame:
    X = as_matrix(data).astype(float)
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
