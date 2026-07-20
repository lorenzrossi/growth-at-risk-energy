#!/usr/bin/env python3
"""Master script for the GaR project (Python port).

Runs the full pipeline from raw CSVs to figures with a fixed seed.
Usage:  python run_all.py  (from the py_port folder)
"""

from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from gar import (build_datasets, qr_direct, coef_table, step2_match,
                 exceedance_correl, trailing_mean)
from gar.plots import (plot_coef_by_horizon, plot_coef_by_quantile,
                       plot_qr_scatter, plot_fan, plot_es_longrise,
                       plot_inverse_cdf, plot_exceedance)

# ---------------- configuration (EDIT PATHS HERE) ----------------
ROOT = Path(__file__).resolve().parent
CFG = {
    "files": {
        "aggregate": ROOT / "input" / "DataIP_GaR.csv",
        "brent": ROOT / "input" / "DataOIL.csv",
        "panel": ROOT / "input" / "country_sector.csv",
        "eu_sectors": {"energy": ROOT / "input" / "EU27_MIG_energy.csv"},
    },
    "horizons": list(range(1, 13)),
    "horizons_density": [1, 3, 6, 12],   # extend to 1..12 for the final run
    "taus": np.round(np.arange(0.05, 0.951, 0.05), 4),
    "n_boot": 1000,
    "use_cache": True,
    "seed": 123,
}
OUT = ROOT / "output"
CACHE = ROOT / "data" / "cache"
SPECS = [("Diff_GP", "Nat. Gas Price Change"), ("NGPI", "NGPI"),
         ("Diff_BP", "Brent Price Change"), ("NOPI", "NOPI")]


def clean_xy(d: pd.DataFrame, var: str):
    ok = d["Diff_IP"].notna() & d[var].notna()
    return (d.loc[ok, "Diff_IP"].to_numpy(), d.loc[ok, var].to_numpy(),
            d.loc[ok, "Time"].to_numpy())


def stage_01_exceedance(dl):
    fig = OUT / "fig_exceedance"; tab = OUT / "tab_exceedance"
    tab.mkdir(parents=True, exist_ok=True)
    pairs = [("aggregate", "Diff_GP", "EU27_IP_vs_Gas"),
             ("aggregate", "Diff_BP", "EU27_IP_vs_Oil")]
    for geo in ("Germany", "France", "Italy"):
        key = f"{geo}.Manufacturing"
        if key in dl:
            pairs += [(key, "Diff_GP", f"{geo}_IP_vs_Gas"),
                      (key, "Diff_BP", f"{geo}_IP_vs_Oil")]
    for ds, var, name in pairs:
        d = dl[ds]
        y, x, tm = clean_xy(d, var)
        ex = pd.DatetimeIndex(tm).year != 2020
        full = exceedance_correl(y, x)
        ex20 = exceedance_correl(y[ex], x[ex])
        full["table"].to_csv(tab / f"{name}_full.csv", index=False)
        ex20["table"].to_csv(tab / f"{name}_ex2020.csv", index=False)
        plot_exceedance(full["table"], ex20["table"],
                        title=f"{name}  (HTZ p = {full['pval']:.3f})",
                        filename=fig / f"{name}.pdf")
        j90 = (full["table"]["threshold"] - 0.90).abs().idxmin()
        print(f"01: {name:24s} HTZ p={full['pval']:.3f} | "
              f"n@0.9 lower={full['table'].loc[j90, 'n_lower']} "
              f"upper={full['table'].loc[j90, 'n_upper']}")


def stage_02_exploratory(dl):
    fig = OUT / "fig_exploratory"; tab = OUT / "tab_exploratory"
    tab.mkdir(parents=True, exist_ok=True)
    datasets = {"aggregate": "aggregate"}
    if "EU27.energy" in dl:
        datasets["energy"] = "EU27.energy"
    for ds, key in datasets.items():
        d = dl[key]
        for var, label in SPECS:
            if var not in d.columns:
                continue
            y, x, _ = clean_xy(d, var)
            res = qr_direct(y, x, horizons=CFG["horizons"], taus=CFG["taus"],
                            include_lag_y=True, target="point",
                            n_boot=CFG["n_boot"], seed=CFG["seed"],
                            var_names=["Intercept", var, "Diff_IP_lag"])
            coef_table(res, var).to_csv(tab / f"{ds}_{var}_coefs.csv",
                                        index=False)
            for tau in (0.05, 0.95):
                plot_coef_by_horizon(res, var, tau, label,
                    filename=fig / f"{ds}_{var}_Q{round(100*tau):02d}_vs_horizon.pdf")
            for h in [hh for hh in (1, 3, 6, 12) if hh in res]:
                plot_coef_by_quantile(res, var, h, label,
                    filename=fig / f"{ds}_{var}_coef_by_quantile_H{h}.pdf")
                plot_qr_scatter(res, y, x, h, x_label=label,
                    filename=fig / f"{ds}_{var}_scatter_H{h}.pdf")
            print(f"02: {ds} / {var:8s} done (n={res[1]['n']})")


def stage_03_density(dl):
    fig = OUT / "fig_density"
    CACHE.mkdir(parents=True, exist_ok=True)
    focus_dates = [np.datetime64("2020-04-01"), np.datetime64("2022-03-01")]
    datasets = {"AGG": "aggregate"}
    if "EU27.energy" in dl:
        datasets["ENERGY"] = "EU27.energy"
    for ds, key in datasets.items():
        d = dl[key]
        for var, label in SPECS:
            if var not in d.columns:
                continue
            y, x, tm = clean_xy(d, var)
            res_m = qr_direct(y, x, horizons=CFG["horizons_density"],
                              taus=CFG["taus"], include_lag_y=True,
                              target="average",
                              var_names=["Intercept", var, "Diff_IP_lag"])
            res_i = qr_direct(y, None, horizons=CFG["horizons_density"],
                              taus=CFG["taus"], include_lag_y=True,
                              target="average")
            for h in CFG["horizons_density"]:
                cache = CACHE / f"match_{var}_{ds}_H{h}.pkl"
                if cache.exists() and CFG["use_cache"]:
                    m = pickle.loads(cache.read_bytes())
                else:
                    yq_unc = np.nanquantile(trailing_mean(y, h), CFG["taus"])
                    m = {"full": step2_match(res_m[h]["yq"], yq_unc, CFG["taus"]),
                         "iponly": step2_match(res_i[h]["yq"], yq_unc, CFG["taus"])}
                    cache.write_bytes(pickle.dumps(m))
                r = res_m[h]
                plot_es_longrise(m["full"], tm,
                    title=f"ES and Longrise, {var} ({ds}, h = {h})",
                    filename=fig / f"ES_{var}_{ds}_H{h}.pdf")
                plot_fan(r, tm, trailing_mean(y, h),
                    title=f"Predicted distribution, {var} ({ds}, h = {h})",
                    filename=fig / f"Fan_{var}_{ds}_H{h}.pdf")
                for fd in focus_dates:
                    jt = np.flatnonzero(tm == fd)
                    if len(jt) == 1 and jt[0] + h < len(tm):
                        plot_inverse_cdf(m["full"], m["iponly"], r, int(jt[0]), h,
                            labels=(f"IP and {label}", "IP only", "Raw QR"),
                            filename=fig / f"InvCDF_{var}_{ds}_H{h}_"
                                     f"{pd.Timestamp(fd):%Y_%m}.pdf")
                pd.DataFrame({"Time": tm, "GaR05": m["full"]["gar"],
                              "ES": m["full"]["shortfall"],
                              "Longrise": m["full"]["longrise"],
                              "skew": m["full"]["moments"][:, 2]}
                             ).to_csv(fig / f"GaR_{var}_{ds}_H{h}.csv",
                                      index=False)
                print(f"03: {ds} / {var} h={h} done "
                      f"(crossed: {m['full']['n_crossed']})")


def stage_04_subsample(dl):
    fig = OUT / "fig_subsample"; tab = OUT / "tab_subsample"
    tab.mkdir(parents=True, exist_ok=True)
    import matplotlib.pyplot as plt
    datasets = {"aggregate": ("aggregate", 0.05)}
    if "EU27.energy" in dl:
        datasets["energy"] = ("EU27.energy", 0.95)
    for ds, (key, tau_show) in datasets.items():
        d = dl[key]
        yr = np.asarray(pd.DatetimeIndex(d["Time"]).year)
        masks = {"full": np.ones(len(d), bool),
                 "pre2020": yr <= 2019,
                 "ex2020": yr != 2020}
        for var, label in SPECS[:3]:
            if var not in d.columns:
                continue
            tabs = []
            for mn, mask in masks.items():
                dd = d.loc[mask].reset_index(drop=True)
                y, x, _ = clean_xy(dd, var)
                res = qr_direct(y, x, horizons=CFG["horizons"],
                                taus=CFG["taus"], include_lag_y=True,
                                target="point", n_boot=CFG["n_boot"],
                                seed=CFG["seed"],
                                var_names=["Intercept", var, "Diff_IP_lag"])
                t = coef_table(res, var); t["sample"] = mn
                tabs.append(t)
            tab_all = pd.concat(tabs, ignore_index=True)
            tab_all.to_csv(tab / f"{ds}_{var}_subsamples.csv", index=False)

            df = tab_all[np.isclose(tab_all["tau"], tau_show)]
            fig_, ax = plt.subplots(figsize=(7.5, 4.5))
            for mn, g in df.groupby("sample"):
                ax.plot(g["h"], g["coef"], lw=2, label=mn)
                ax.fill_between(g["h"], g["coef"] - 1.96 * g["se"],
                                g["coef"] + 1.96 * g["se"], alpha=.12)
            ax.axhline(0, ls="--", c="k", lw=1)
            ax.set(xlabel="Horizon (months)", ylabel=f"{label} coefficient",
                   title=f"{label}, Q {tau_show:.2f}: subsample comparison ({ds})")
            ax.legend()
            (fig / f"{ds}_{var}_Q{round(100*tau_show):02d}_subsamples.pdf"
             ).parent.mkdir(parents=True, exist_ok=True)
            fig_.savefig(fig / f"{ds}_{var}_Q{round(100*tau_show):02d}_subsamples.pdf",
                         bbox_inches="tight")
            plt.close(fig_)
            print(f"04: {ds} / {var} done")


def stage_05_countries(dl):
    fig = OUT / "fig_countries"; tab = OUT / "tab_countries"
    tab.mkdir(parents=True, exist_ok=True)
    for geo in ("Germany", "Italy", "France"):
        keys = [k for k in dl if k.startswith(f"{geo}.")]
        for sec_label, pat in (("Manufacturing", "Manufacturing"),
                               ("Energy", "energy")):
            key = next((k for k in keys if pat.lower() in k.lower()), None)
            if key is None:
                continue
            d = dl[key]
            for var, label in SPECS:
                if var not in d.columns:
                    continue
                y, x, _ = clean_xy(d, var)
                res = qr_direct(y, x, horizons=CFG["horizons"],
                                taus=CFG["taus"], include_lag_y=True,
                                target="point", n_boot=CFG["n_boot"],
                                seed=CFG["seed"],
                                var_names=["Intercept", var, "Diff_IP_lag"])
                coef_table(res, var).to_csv(
                    tab / f"{geo}_{sec_label}_{var}_coefs.csv", index=False)
                for tau in (0.05, 0.95):
                    plot_coef_by_horizon(res, var, tau,
                        f"{geo} {sec_label} - {label}",
                        filename=fig / f"{geo}_{sec_label}_{var}_"
                                       f"Q{round(100*tau):02d}.pdf")
                print(f"05: {geo} / {sec_label} / {var} done")


def main():
    np.random.seed(CFG["seed"])
    print("=== 00: data construction ===")
    dl = build_datasets(CFG["files"])
    (ROOT / "data").mkdir(exist_ok=True)
    print(f"    {len(dl)} datasets: {', '.join(sorted(dl))}")

    print("=== 01: exceedance correlations ===")
    stage_01_exceedance(dl)
    print("=== 02: exploratory quantile regressions ===")
    stage_02_exploratory(dl)
    print("=== 03: skew-t densities / GaR ===")
    stage_03_density(dl)
    print("=== 04: subsample comparison ===")
    stage_04_subsample(dl)
    print("=== 05: country analysis ===")
    stage_05_countries(dl)
    print("=== run_all complete: see output/ ===")


if __name__ == "__main__":
    main()
