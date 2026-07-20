"""Direct-projection quantile regressions with moving-block bootstrap.

Replaces: QRboot.m, run_static_qr.m (rq.m -> statsmodels QuantReg)."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.regression.quantile_regression import QuantReg
from statsmodels.tools.sm_exceptions import IterationLimitWarning

from .utils import block_boot_idx, trailing_mean


def _fit_qr(y: np.ndarray, Z: np.ndarray, tau: float) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IterationLimitWarning)
        res = QuantReg(y, Z).fit(q=tau, max_iter=5000, p_tol=1e-7)
    return np.asarray(res.params)


def qr_direct(y, X=None, horizons=range(1, 13),
              taus=np.arange(0.05, 0.951, 0.05),
              include_lag_y: bool = True, target: str = "point",
              n_boot: int = 0, block_length: int | None = None,
              var_names=None, seed: int | None = None) -> dict:
    """For each horizon h, regress the target at t+h on conditioning
    variables at t.

    target = "point"   : y(t+h)                     (coefficient figures)
    target = "average" : mean(y(t+1..t+h))          (Adrian et al.; density/GaR)

    Returns {h: {"coef", "se", "ols", "yq", "boot", "n", "taus",
                 "var_names", "target"}} where coef/se are (k x ntau) frames
    and yq is a (T x ntau) array of fitted quantiles aligned so row t is the
    fitted quantile of the target dated t (first h rows NaN).
    """
    assert target in ("point", "average")
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=float)
    n = len(y)
    taus = np.asarray(taus, dtype=float)

    cols = [np.ones(n)]
    names = ["Intercept"]
    if X is not None:
        X = np.atleast_2d(np.asarray(X, dtype=float))
        if X.shape[0] != n:
            X = X.T
        for j in range(X.shape[1]):
            cols.append(X[:, j])
            names.append(f"X{j + 1}")
    if include_lag_y:
        cols.append(y)
        names.append("y_lag")
    if var_names is not None:
        names = list(var_names)
    Z = np.column_stack(cols)
    k = Z.shape[1]

    out = {}
    for h in horizons:
        yt = y if target == "point" else trailing_mean(y, h)
        yh = yt[h:]
        Zh = Z[:-h]
        ok = np.isfinite(yh) & np.all(np.isfinite(Zh), axis=1)
        yh_f, Zh_f = yh[ok], Zh[ok]
        nh = len(yh_f)

        coef = np.column_stack([_fit_qr(yh_f, Zh_f, t) for t in taus])
        ols = np.linalg.lstsq(Zh_f, yh_f, rcond=None)[0]

        yq = np.full((n, len(taus)), np.nan)
        target_rows = np.arange(h, n)[ok]
        yq[target_rows] = Zh_f @ coef

        se = np.full_like(coef, np.nan)
        boot = None
        if n_boot > 0:
            bl = max(12, h) if block_length is None else block_length
            boot = np.full((n_boot, k, len(taus)), np.nan)
            for b in range(n_boot):
                idx = block_boot_idx(nh, bl, rng)
                for j, t in enumerate(taus):
                    try:
                        boot[b, :, j] = _fit_qr(yh_f[idx], Zh_f[idx], t)
                    except Exception:
                        pass
            se = np.nanstd(boot, axis=0)

        out[h] = {
            "h": h,
            "coef": pd.DataFrame(coef, index=names,
                                 columns=[f"tau{t:.2f}" for t in taus]),
            "se": pd.DataFrame(se, index=names,
                               columns=[f"tau{t:.2f}" for t in taus]),
            "ols": pd.Series(ols, index=names),
            "yq": yq, "boot": boot, "n": nh, "taus": taus,
            "var_names": names, "target": target,
        }
    return out


def coef_table(res: dict, var: str) -> pd.DataFrame:
    """Tidy coefficient/SE table across horizons and quantiles for one
    regressor, with normal-approximation p-values."""
    rows = []
    for h, r in res.items():
        for j, tau in enumerate(r["taus"]):
            c = r["coef"].loc[var].iloc[j]
            s = r["se"].loc[var].iloc[j]
            rows.append({"h": h, "tau": round(float(tau), 4),
                         "coef": c, "se": s})
    tab = pd.DataFrame(rows)
    from scipy.stats import norm
    tab["pval"] = 2 * norm.sf(np.abs(tab["coef"] / tab["se"]))
    tab["sig"] = pd.cut(tab["pval"], [-np.inf, .01, .05, .10, np.inf],
                        labels=["***", "**", "*", ""])
    return tab
