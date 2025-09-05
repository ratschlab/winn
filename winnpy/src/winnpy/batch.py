from __future__ import annotations
import numpy as np
import pandas as pd
from statsmodels.api import OLS, add_constant

# optional: pyComBat
_pycombat = None
try:
    from combat.pycombat import pycombat as _pycombat
except Exception:
    pass

from .utils import as_matrix, fdr_bh

def anova_batch_correction(
    data: pd.DataFrame,
    batch: np.ndarray,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Perform one-way ANOVA in log1p space (equivalent to y ~ C(batch)), only impute batch mean differences for significant rows."""
    X = as_matrix(data).astype(float)
    batch = pd.Categorical(batch)
    if X.shape[1] != len(batch):
        raise ValueError("Length of batch must equal number of columns in data.")

    Z = np.log1p(X)
    # drop_first=True to avoid perfect collinearity; overall F-test as batch effect significance
    B = pd.get_dummies(batch, drop_first=True)
    Xdesign = add_constant(B.values.astype(float), has_constant="add")  # ← 关键：astype(float)
    pvals = []
    for _, row in Z.iterrows():
        y = row.values
        model = OLS(y, Xdesign).fit()
        pvals.append(float(model.f_pvalue) if np.isfinite(model.f_pvalue) else 1.0)


    padj = fdr_bh(pvals, alpha=fdr_threshold)

    Y = Z.copy()
    overall_means = Z.mean(axis=1)
    df_batch = pd.Series(batch).astype(str).values
    sig_idx = np.where(padj < fdr_threshold)[0]
    for ridx in sig_idx:
        vals = Z.iloc[ridx, :].values
        batch_means = pd.Series(vals).groupby(df_batch).mean()
        shift = batch_means - overall_means.iloc[ridx]
        Y.iloc[ridx, :] = vals - pd.Index(df_batch).map(shift).values
    return pd.DataFrame(np.expm1(Y.values), index=X.index, columns=X.columns)

def combat_batch_correction(
    data: pd.DataFrame,
    batch: np.ndarray,
    par_prior: bool = True,
    mean_only: bool = False,
    ref_batch = None,
) -> pd.DataFrame:
    if _pycombat is None:
        raise ImportError("pyComBat not found. Install `combat` (pip install combat) or use ANOVA.")
    X = as_matrix(data).astype(float)
    batch = pd.Categorical(batch)
    Z = np.log1p(X)
    Zdf = pd.DataFrame(Z.values, index=Z.index, columns=Z.columns)  # ← 关键：DataFrame
    corrected_df = _pycombat(
        data=Zdf,
        batch=batch.codes,
        par_prior=par_prior,
        mean_only=mean_only,
        ref_batch=ref_batch,
    )
    corrected = corrected_df.values if hasattr(corrected_df, "values") else np.asarray(corrected_df)
    return pd.DataFrame(np.expm1(corrected), index=X.index, columns=X.columns)

