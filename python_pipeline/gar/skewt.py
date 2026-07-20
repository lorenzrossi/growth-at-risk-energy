"""Azzalini-Capitanio skew-t distribution and quantile-matching (step 2 of
Adrian-Boyarchenko-Giannone).

Replaces: Step2match_improved.m, QuantilesInterpolation_improved.m and the
azzalini/ MATLAB suite (dskt/pskt/qskt).

The quantile function is obtained by inverting a grid-based CDF (trapezoidal
integration of the exact density); accurate to ~1e-3 on the standardized
scale, ample for quantile matching. Fitting exploits that for fixed
(alpha, nu) the optimal (xi, omega) solves a LINEAR least-squares problem,
so only alpha requires numerical optimisation -- faster and more robust than
the 3-parameter search in the MATLAB code.
"""

from __future__ import annotations

import numpy as np
from scipy import stats, optimize
from scipy.integrate import cumulative_trapezoid


def skewt_pdf_std(z: np.ndarray, alpha: float, nu: float) -> np.ndarray:
    """Density of the standardized (xi=0, omega=1) Azzalini skew-t."""
    z = np.asarray(z, dtype=float)
    w = alpha * z * np.sqrt((nu + 1.0) / (nu + z ** 2))
    return 2.0 * stats.t.pdf(z, nu) * stats.t.cdf(w, nu + 1.0)


def skewt_pdf(x, xi, omega, alpha, nu):
    return skewt_pdf_std((np.asarray(x, float) - xi) / omega, alpha, nu) / omega


class _StdQuantiles:
    """Cache of standardized skew-t quantile functions on a z-grid."""

    def __init__(self, z_lo: float = -80.0, z_hi: float = 80.0, n: int = 8001):
        self.z = np.linspace(z_lo, z_hi, n)
        self._cache: dict[tuple[float, float], np.ndarray] = {}

    def cdf_grid(self, alpha: float, nu: float) -> np.ndarray:
        key = (round(float(alpha), 6), float(nu))
        if key not in self._cache:
            pdf = skewt_pdf_std(self.z, alpha, nu)
            cdf = cumulative_trapezoid(pdf, self.z, initial=0.0)
            cdf /= cdf[-1]
            self._cache[key] = cdf
        return self._cache[key]

    def ppf(self, p, alpha: float, nu: float) -> np.ndarray:
        cdf = self.cdf_grid(alpha, nu)
        return np.interp(np.asarray(p, float), cdf, self.z)


_STD = _StdQuantiles()


def skewt_ppf(p, xi, omega, alpha, nu):
    """Quantile function of the Azzalini skew-t."""
    return xi + omega * _STD.ppf(p, alpha, nu)


def fit_skewt_quantiles(q_targ: np.ndarray, taus: np.ndarray,
                        q_match=(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
                        nu_grid=range(2, 21)) -> dict:
    """Fit (xi, omega, alpha, nu) so the skew-t quantiles match q_targ.

    Profiles over integer nu (as in the AER replication files); for each
    (alpha, nu) the optimal (xi, omega >= 0) is a linear LS fit of
    q_sel on [1, z_tau(alpha, nu)]. Crossed target quantiles are rearranged
    (Chernozhukov-Fernandez-Val-Galichon 2010) and flagged.
    """
    q_targ = np.asarray(q_targ, dtype=float)
    taus = np.asarray(taus, dtype=float)

    crossed = bool(np.any(np.diff(q_targ) < 0))
    if crossed:
        q_targ = np.sort(q_targ)

    sel = [int(np.argmin(np.abs(taus - q))) for q in q_match]
    p_sel = taus[sel]
    q_sel = q_targ[sel]

    def inner_ls(alpha: float, nu: float):
        zq = _STD.ppf(p_sel, alpha, nu)
        A = np.column_stack([np.ones_like(zq), zq])
        coef, *_ = np.linalg.lstsq(A, q_sel, rcond=None)
        xi, omega = coef
        omega = max(omega, 1e-4)          # scale must be positive
        resid = q_sel - (xi + omega * zq)
        return float(resid @ resid), xi, omega

    best = {"ssq": np.inf}
    for nu in nu_grid:
        r = optimize.minimize_scalar(lambda a: inner_ls(a, nu)[0],
                                     bounds=(-60.0, 60.0), method="bounded",
                                     options={"xatol": 1e-4})
        ssq, xi, omega = inner_ls(r.x, nu)
        if ssq < best["ssq"]:
            best = {"xi": xi, "omega": omega, "alpha": float(r.x),
                    "nu": float(nu), "ssq": ssq}
    best["crossed"] = crossed
    return best


def step2_match(YQ: np.ndarray, yq_unc: np.ndarray, taus: np.ndarray,
                yy: np.ndarray | None = None, alpha_es: float = 0.05,
                verbose: bool = False) -> dict:
    """Fit a skew-t to each row of conditional quantiles; compute densities,
    quantiles, moments, entropy, GaR and Expected Shortfall / Longrise.

    Mirrors Step2match's Res struct. Rows containing NaN are skipped.
    """
    if yy is None:
        yy = np.arange(-30.0, 30.0 + 1e-9, 0.02)
    taus = np.asarray(taus, dtype=float)
    Tn, ntau = YQ.shape
    delta = yy[1] - yy[0]
    jq50 = int(np.argmin(np.abs(taus - 0.5)))
    p_grid = np.arange(0.01, alpha_es + 1e-12, 0.01)   # ES integral, delta1=0.01

    un = fit_skewt_quantiles(yq_unc, taus)
    p_unc = skewt_pdf(yy, un["xi"], un["omega"], un["alpha"], un["nu"])

    out = {
        "par": np.full((Tn, 4), np.nan), "PST": np.full((Tn, len(yy)), np.nan),
        "QST": np.full((Tn, ntau), np.nan), "yy": yy,
        "moments": np.full((Tn, 4), np.nan), "entropy": np.full((Tn, 2), np.nan),
        "shortfall": np.full(Tn, np.nan), "longrise": np.full(Tn, np.nan),
        "gar": np.full(Tn, np.nan), "n_crossed": 0, "unc": un, "taus": taus,
    }

    for t in range(Tn):
        q_targ = YQ[t]
        if np.any(np.isnan(q_targ)):
            continue
        ft = fit_skewt_quantiles(q_targ, taus)
        if not np.isfinite(ft["ssq"]):
            continue
        if ft["crossed"]:
            out["n_crossed"] += 1
        xi, om, al, nu = ft["xi"], ft["omega"], ft["alpha"], ft["nu"]
        out["par"][t] = [xi, om, al, nu]

        d = skewt_pdf(yy, xi, om, al, nu)
        out["PST"][t] = d
        out["QST"][t] = skewt_ppf(taus, xi, om, al, nu)

        m1 = np.sum(yy * d) * delta
        v = np.sum((yy - m1) ** 2 * d) * delta
        out["moments"][t] = [m1, v,
                             np.sum((yy - m1) ** 3 * d) * delta / v ** 1.5,
                             np.sum((yy - m1) ** 4 * d) * delta / v ** 2]

        med = out["QST"][t, jq50]
        with np.errstate(divide="ignore", invalid="ignore"):
            lg = np.log(d) - np.log(p_unc)
        lg[~np.isfinite(lg)] = 0.0
        out["entropy"][t] = [np.sum(lg * d * (yy < med)) * delta,
                             np.sum(lg * d * (yy > med)) * delta]

        out["shortfall"][t] = np.mean(skewt_ppf(p_grid, xi, om, al, nu))
        out["longrise"][t] = np.mean(skewt_ppf(1.0 - p_grid, xi, om, al, nu))
        out["gar"][t] = skewt_ppf(alpha_es, xi, om, al, nu)

        if verbose and (t + 1) % 24 == 0:
            print(f"  skew-t fit: {t + 1} / {Tn}")

    return out
