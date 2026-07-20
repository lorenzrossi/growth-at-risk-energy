"""Data construction. Replaces ReadData.m, ReadData_OIL.m,
Country_Read_Data.m.

All growth rates and NGPI/NOPI are recomputed from raw IP and price LEVELS;
precomputed columns in the CSVs (including the corrupted NGPI column produced
by the old ReadData.m) are ignored by design. Oil variables are merged into
every IP dataset by DATE (never by row order)."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .utils import log_diff, net_increase, subsample_window


def _read_time_csv(path) -> pd.DataFrame:
    df = pd.read_csv(path)
    tcol = next(c for c in ("Time", "Date", "TIME_PERIOD") if c in df.columns)
    df["Time"] = pd.to_datetime(df[tcol])
    return df


def _prep_series(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values("Time").reset_index(drop=True)
    out = pd.DataFrame({"Time": df["Time"], "IP": df["IP"]})
    out["Diff_IP"] = log_diff(df["IP"].to_numpy())
    if "GAS_PRICE" in df.columns:
        out["GAS_PRICE"] = df["GAS_PRICE"]
        out["Diff_GP"] = log_diff(df["GAS_PRICE"].to_numpy())
        out["NGPI"] = net_increase(df["GAS_PRICE"].to_numpy())
    return subsample_window(out)


def _key(geo: str, sector: str) -> str:
    geo_c = re.sub(r"\W+", "_", geo)
    sec_c = re.sub(r"\W+", "_", sector)
    return f"{geo_c}.{sec_c}"


def build_datasets(files: dict) -> dict[str, pd.DataFrame]:
    """files: {"aggregate": path, "brent": path, "panel": path,
               "eu_sectors": {name: path, ...}}"""
    data = {}

    agg = _read_time_csv(files["aggregate"])
    data["aggregate"] = _prep_series(agg)

    oil_raw = _read_time_csv(files["brent"]).sort_values("Time")
    oil = pd.DataFrame({
        "Time": oil_raw["Time"],
        "BRENT_PRICE": oil_raw["BRENT_PRICE"],
        "Diff_BP": log_diff(oil_raw["BRENT_PRICE"].to_numpy()),
        "NOPI": net_increase(oil_raw["BRENT_PRICE"].to_numpy()),
    })
    data["oil"] = subsample_window(oil)

    panel = _read_time_csv(files["panel"])
    for geo in panel["geo"].unique():
        for sec in panel["nace_r2"].unique():
            sub = panel[(panel["geo"] == geo) & (panel["nace_r2"] == sec)]
            if len(sub) < 24:
                continue
            data[_key(geo, sec)] = _prep_series(sub)

    for name, path in files.get("eu_sectors", {}).items():
        if not Path(path).exists():
            print(f"  [warn] EU sector file not found, skipping: {path}")
            continue
        data[f"EU27.{name}"] = _prep_series(_read_time_csv(path))

    oil_cols = data["oil"][["Time", "BRENT_PRICE", "Diff_BP", "NOPI"]]
    for name in list(data):
        if name == "oil":
            continue
        data[name] = (data[name].merge(oil_cols, on="Time", how="left")
                      .sort_values("Time").reset_index(drop=True))

    return data
