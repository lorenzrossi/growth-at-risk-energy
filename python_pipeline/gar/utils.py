"""Shared utilities for the GaR project (Python port).

Replaces: netIncrease.m, blockBootIdx.m."""

from __future__ import annotations

import numpy as np
import pandas as pd


def net_increase(P: np.ndarray, n_months: int = 12) -> np.ndarray:
    """Hamilton (2003) net price increase on the 100*log scale.

    NI[t] = max(0, 100*log(P[t]) - max(100*log(P[t-12:t])))

    Single canonical construction for both NGPI and NOPI. ``P`` must be
    strictly positive PRICE LEVELS (not logs, not differences).
    """
    P = np.asarray(P, dtype=float)
    if np.any(P[~np.isnan(P)] <= 0):
        raise ValueError("P must contain strictly positive price levels.")
    logP = 100.0 * np.log(P)
    n = len(logP)
    NI = np.zeros(n)
    for t in range(n_months, n):
        NI[t] = max(0.0, logP[t] - np.nanmax(logP[t - n_months:t]))
    return NI


def block_boot_idx(n: int, block_length: int = 12,
                   rng: np.random.Generator | None = None) -> np.ndarray:
    """Moving-block bootstrap indices (Kunsch, 1989).

    Use instead of iid resampling whenever residuals are serially dependent
    (direct multi-horizon projections in particular). Recommended
    block_length: max(12, h) for monthly data at horizon h.
    """
    if rng is None:
        rng = np.random.default_rng()
    block_length = int(min(max(1, round(block_length)), n))
    n_blocks = int(np.ceil(n / block_length))
    starts = rng.integers(0, n - block_length + 1, n_blocks)
    idx = np.concatenate([np.arange(s, s + block_length) for s in starts])
    return idx[:n]


def log_diff(x: np.ndarray) -> np.ndarray:
    """100*log first difference (monthly log growth in percent)."""
    x = np.asarray(x, dtype=float)
    out = np.full_like(x, np.nan)
    out[1:] = 100.0 * np.diff(np.log(x))
    return out


def trailing_mean(y: np.ndarray, h: int) -> np.ndarray:
    """Trailing h-month mean, aligned like MATLAB filter(ones(1,h)/h, 1, y):
    out[t] = mean(y[t-h+1 : t+1]), NaN for the first h-1 observations."""
    y = np.asarray(y, dtype=float)
    if h == 1:
        return y.copy()
    return pd.Series(y).rolling(h).mean().to_numpy()


def subsample_window(df: pd.DataFrame, start: str = "2004-01-01",
                     end: str = "2023-12-31") -> pd.DataFrame:
    m = (df["Time"] >= pd.Timestamp(start)) & (df["Time"] <= pd.Timestamp(end))
    return df.loc[m].reset_index(drop=True)
