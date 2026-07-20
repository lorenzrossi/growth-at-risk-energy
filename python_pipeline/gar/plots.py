"""Matplotlib figure builders. Replaces PlotQRbands.m, PlotPredictiveTS.m and
the per-script plotting code."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .utils import trailing_mean


def _save(fig, filename):
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filename, bbox_inches="tight")
    plt.close(fig)
    return filename


def plot_coef_by_horizon(res, var, tau, var_label=None, filename=None):
    """Coefficient across horizons at fixed tau, 68/95% block-bootstrap bands."""
    var_label = var_label or var
    hs, c, s = [], [], []
    for h, r in sorted(res.items()):
        j = int(np.argmin(np.abs(r["taus"] - tau)))
        hs.append(h)
        c.append(r["coef"].loc[var].iloc[j])
        s.append(r["se"].loc[var].iloc[j])
    hs, c, s = map(np.asarray, (hs, c, s))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.fill_between(hs, c - 1.96 * s, c + 1.96 * s, color="0.8", label="95% CI")
    ax.fill_between(hs, c - s, c + s, color="0.6", label="68% CI")
    ax.plot(hs, c, "r-", lw=2, label="Estimate")
    ax.axhline(0, ls="--", c="k", lw=1)
    ax.set(xlabel="Horizon (months)", ylabel=f"{var_label} coefficient",
           title=f"{var_label} coefficients (Q {tau:.2f})")
    ax.legend()
    return _save(fig, filename) if filename else fig


def plot_coef_by_quantile(res, var, h, var_label=None, filename=None):
    var_label = var_label or var
    r = res[h]
    taus, c = r["taus"], r["coef"].loc[var].to_numpy()
    s = r["se"].loc[var].to_numpy()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if np.all(np.isfinite(s)):
        ax.fill_between(taus, c - 1.96 * s, c + 1.96 * s, color="0.8")
        ax.fill_between(taus, c - s, c + s, color="0.6")
    ax.plot(taus, c, "r-", lw=2)
    ax.axhline(0, ls="--", c="k", lw=1)
    ax.set(xlabel=r"$\tau$", ylabel=f"{var_label} coefficient",
           title=f"{var_label} coefficients across quantiles (h = {h})")
    return _save(fig, filename) if filename else fig


def plot_qr_scatter(res, y, x, h, x_label="X", y_label="IP growth",
                    taus_show=(0.05, 0.5, 0.95), filename=None):
    r = res[h]
    n = len(y)
    yv = y if r["target"] == "point" else trailing_mean(np.asarray(y), h)
    xs = np.asarray(x)[: n - h]
    ys = np.asarray(yv)[h:]
    trows = np.arange(h, n)
    ok = np.isfinite(xs) & np.isfinite(ys)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(xs[ok], ys[ok], s=22, alpha=0.6, color="steelblue")
    i_lo, i_hi = np.nanargmin(xs), np.nanargmax(xs)
    for tq in taus_show:
        j = int(np.argmin(np.abs(r["taus"] - tq)))
        yf = r["yq"][trows, j]
        ax.plot([xs[i_lo], xs[i_hi]], [yf[i_lo], yf[i_hi]],
                lw=2, label=f"Q{round(100 * tq):02d}",
                ls="-" if tq == 0.5 else "--")
    ax.set(xlabel=x_label, ylabel=f"{y_label} {h} month(s) ahead",
           title=f"Quantile regressions (h = {h})")
    ax.legend()
    return _save(fig, filename) if filename else fig


def plot_fan(r, time, y_real, title="Predicted distribution", filename=None):
    taus = r["taus"]
    j = {q: int(np.argmin(np.abs(taus - q))) for q in (.05, .25, .5, .75, .95)}
    yq = r["yq"]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.fill_between(time, yq[:, j[.05]], yq[:, j[.95]], color="0.85")
    ax.fill_between(time, yq[:, j[.25]], yq[:, j[.75]], color="0.65")
    ax.plot(time, yq[:, j[.5]], "k-", lw=1.2, label="Median")
    ax.plot(time, y_real, "r-", lw=0.8, label="Realised")
    ax.set(ylabel="IPI growth", title=title)
    ax.legend()
    return _save(fig, filename) if filename else fig


def plot_es_longrise(match, time, title="Expected Shortfall and Longrise",
                     filename=None):
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(time, match["shortfall"], "-", lw=1.2, label="Shortfall")
    ax.plot(time, match["longrise"], "--", lw=1.2, label="Longrise")
    ax.set(title=title)
    ax.legend()
    return _save(fig, filename) if filename else fig


def plot_inverse_cdf(match_full, match_iponly, r_full, date_index, h,
                     labels=("Energy-augmented", "IP only", "Raw QR"),
                     filename=None):
    taus = match_full["taus"]
    jt = date_index + h
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(taus, match_full["QST"][jt], lw=2, label=labels[0])
    ax.plot(taus, match_iponly["QST"][jt], "--", lw=2, label=labels[1])
    ax.plot(taus, r_full["yq"][jt], "-*", lw=1, label=labels[2])
    ax.set(xlabel=r"$\tau$", ylabel="IPI growth",
           title=f"Fitted inverse CDF (h = {h})")
    ax.grid(alpha=.3)
    ax.legend()
    return _save(fig, filename) if filename else fig


def plot_exceedance(tab_full, tab_ex=None, title="Exceedance correlations",
                    filename=None):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(tab_full["threshold"], tab_full["lower"], "-o", label="Lower (full)")
    ax.plot(tab_full["threshold"], tab_full["upper"], "-s", label="Upper (full)")
    if tab_ex is not None:
        ax.plot(tab_ex["threshold"], tab_ex["lower"], "--o", label="Lower (ex-2020)")
        ax.plot(tab_ex["threshold"], tab_ex["upper"], "--s", label="Upper (ex-2020)")
    ax.set(xlabel="Quantile threshold", ylabel="Exceedance correlation",
           title=title)
    ax.grid(alpha=.3)
    ax.legend()
    return _save(fig, filename) if filename else fig
