"""Exceedance correlations (Ang & Chen 2002) with the Hong-Tu-Zhou (2003)
asymmetry test. Replaces exceedence_correl.m (fixed), corrcoef12.m,
newey_west.m."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2


def newey_west_cov(Xi: np.ndarray, lags: int | None = None) -> np.ndarray:
    Tn = Xi.shape[0]
    if lags is None:
        lags = int(np.floor(4 * (Tn / 100.0) ** (2.0 / 9.0)))
    Xi = Xi - Xi.mean(axis=0)
    om = Xi.T @ Xi / Tn
    for l in range(1, lags + 1):
        w = 1.0 - l / (lags + 1.0)
        G = Xi[l:].T @ Xi[:-l] / Tn
        om += w * (G + G.T)
    return om


def exceedance_correl(X, Y, qc=None, min_n: int = 3) -> dict:
    """Lower tail: corr(X, Y | X <= Qx(1-q), Y <= Qy(1-q));
    Upper tail: corr(X, Y | X >= Qx(q),  Y >= Qy(q)),  q in [0.5, 0.95].

    Returns {"table": DataFrame(threshold, lower, upper, n_lower, n_upper),
             "teststat", "pval"}.
    """
    if qc is None:
        qc = np.arange(0.5, 0.951, 0.05)
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    ok = np.isfinite(X) & np.isfinite(Y)
    X, Y = X[ok], Y[ok]
    Tn = len(X)

    rows, Xi_cols, rhodiff = [], [], []
    for q in qc:
        lo = (X <= np.quantile(X, 1 - q)) & (Y <= np.quantile(Y, 1 - q))
        up = (X >= np.quantile(X, q)) & (Y >= np.quantile(Y, q))
        n_lo, n_up = int(lo.sum()), int(up.sum())
        r_lo = r_up = np.nan
        if n_lo >= min_n and n_up >= min_n:
            r_lo = float(np.corrcoef(X[lo], Y[lo])[0, 1])
            r_up = float(np.corrcoef(X[up], Y[up])[0, 1])
            zx_up = (X - X[up].mean()) / X[up].std(ddof=1)
            zy_up = (Y - Y[up].mean()) / Y[up].std(ddof=1)
            zx_lo = (X - X[lo].mean()) / X[lo].std(ddof=1)
            zy_lo = (Y - Y[lo].mean()) / Y[lo].std(ddof=1)
            xi = (Tn / n_up * (zx_up * zy_up - r_up) * up
                  - Tn / n_lo * (zx_lo * zy_lo - r_lo) * lo)
            Xi_cols.append(xi)
            rhodiff.append(r_up - r_lo)
        rows.append({"threshold": round(float(q), 4), "lower": r_lo,
                     "upper": r_up, "n_lower": n_lo, "n_upper": n_up})

    teststat = pval = np.nan
    if rhodiff:
        Om = newey_west_cov(np.column_stack(Xi_cols))
        rd = np.asarray(rhodiff)
        try:
            teststat = float(Tn * rd @ np.linalg.solve(Om, rd))
            pval = float(chi2.sf(teststat, df=len(rd)))
        except np.linalg.LinAlgError:
            pass

    return {"table": pd.DataFrame(rows), "teststat": teststat, "pval": pval}
