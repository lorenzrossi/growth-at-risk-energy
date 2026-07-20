#!/usr/bin/env python3
"""Numerical validation of the Python port against reference values verified
during the code review (which matched the magnitudes in the MATLAB draft).

Run from py_port:  python validate_core.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from gar import build_datasets, qr_direct, fit_skewt_quantiles, skewt_ppf

ROOT = Path(__file__).resolve().parent
FILES = {
    "aggregate": ROOT / "input" / "DataIP_GaR.csv",
    "brent": ROOT / "input" / "DataOIL.csv",
    "panel": ROOT / "input" / "country_sector.csv",
    "eu_sectors": {},
}
TAUS = np.round(np.arange(0.05, 0.951, 0.05), 4)

results = []


def check(name, value, ref, tol):
    ok = np.isfinite(value) and abs(value - ref) <= tol
    print(f"{name:<55s} py = {value:+.4f} | ref = {ref:+.4f} | "
          f"{'PASS' if ok else 'FAIL'}")
    results.append(ok)


def main():
    np.random.seed(123)
    dl = build_datasets(FILES)
    d = dl["aggregate"]

    # 1. unconditional correlations
    ok = d["Diff_IP"].notna() & d["Diff_GP"].notna() & d["Diff_BP"].notna()
    check("corr(Diff_IP, Diff_GP), 2004-2023",
          d.loc[ok, "Diff_IP"].corr(d.loc[ok, "Diff_GP"]), 0.077, 0.02)
    check("corr(Diff_IP, Diff_BP), 2004-2023",
          d.loc[ok, "Diff_IP"].corr(d.loc[ok, "Diff_BP"]), 0.499, 0.02)

    # 2. NGPI construction
    check("max(NGPI) aggregate (100*log scale)",
          float(d["NGPI"].max()), 37.49, 0.5)

    # 3. aggregate Q05 NGPI coefficients
    okn = d["Diff_IP"].notna() & d["NGPI"].notna()
    y = d.loc[okn, "Diff_IP"].to_numpy()
    x = d.loc[okn, "NGPI"].to_numpy()
    res = qr_direct(y, x, horizons=[1, 12], taus=TAUS, include_lag_y=True,
                    target="point", var_names=["Intercept", "NGPI", "lag"])
    j05 = int(np.argmin(np.abs(TAUS - 0.05)))
    check("Aggregate Q05 NGPI coefficient, h = 1",
          res[1]["coef"].loc["NGPI"].iloc[j05], -0.073, 0.03)
    check("Aggregate Q05 NGPI coefficient, h = 12",
          res[12]["coef"].loc["NGPI"].iloc[j05], -0.196, 0.04)

    # 4. Germany energy sector Q95 NGPI at h = 12
    gkey = next((k for k in dl if k.startswith("Germany.")
                 and "energy" in k.lower()), None)
    if gkey:
        g = dl[gkey]
        okg = g["Diff_IP"].notna() & g["NGPI"].notna()
        resg = qr_direct(g.loc[okg, "Diff_IP"].to_numpy(),
                         g.loc[okg, "NGPI"].to_numpy(),
                         horizons=[12], taus=TAUS, include_lag_y=True,
                         target="point",
                         var_names=["Intercept", "NGPI", "lag"])
        j95 = int(np.argmin(np.abs(TAUS - 0.95)))
        check("Germany energy Q95 NGPI coefficient, h = 12",
              resg[12]["coef"].loc["NGPI"].iloc[j95], -0.145, 0.04)
    else:
        print("Germany energy dataset not found -- skipped check 4")

    # 5. averaged target equals point target at h = 1
    res_avg = qr_direct(y, x, horizons=[1], taus=TAUS, include_lag_y=True,
                        target="average", var_names=["Intercept", "NGPI", "lag"])
    check("Averaged == point target at h = 1 (Q05 NGPI)",
          res_avg[1]["coef"].loc["NGPI"].iloc[j05],
          res[1]["coef"].loc["NGPI"].iloc[j05], 1e-8)

    # 6. skew-t round trip: fit on known quantiles, recover them
    true = dict(xi=-0.5, omega=2.0, alpha=-3.0, nu=6.0)
    q_true = skewt_ppf(TAUS, **true)
    ft = fit_skewt_quantiles(q_true, TAUS)
    q_fit = skewt_ppf(TAUS, ft["xi"], ft["omega"], ft["alpha"], ft["nu"])
    check("Skew-t round trip: max |q_fit - q_true|",
          float(np.max(np.abs(q_fit - q_true))), 0.0, 0.05)

    print(f"\n{sum(results)} of {len(results)} checks passed.")
    if not all(results):
        raise SystemExit("Validation failed -- inspect the FAIL lines above.")


if __name__ == "__main__":
    main()
