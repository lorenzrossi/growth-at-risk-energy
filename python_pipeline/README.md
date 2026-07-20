# GaR project — Python pipeline

Python implementation of "The effects of gas and oil prices on European
industrial production: a Growth-at-Risk approach". Same structure and
conventions as the R pipeline; all fixes from the code review are built in
(canonical 100·log NGPI/NOPI, moving-block bootstrap with block = max(12, h),
spec-specific caches, quantile-crossing rearrangement with a count, explicit
date joins, fixed seed).

**Validation status: the numerical core has been executed and verified against
the project data** — `validate_core.py` passes 8/8 checks (correlations, NGPI
construction, aggregate Q05 NGPI coefficients at h = 1/12, Germany energy Q95
at h = 12, point-vs-average target equivalence, and a skew-t
fit-and-recover round trip).

## Setup

```sh
pip install -r requirements.txt      # numpy pandas scipy statsmodels matplotlib
```

## Layout

```
py_port/
├── run_all.py           master pipeline + configuration (EDIT PATHS AT TOP)
├── validate_core.py     numerical checks — run this first
├── requirements.txt
├── input/               aggregate + oil + country/sector CSVs included;
│                        add EU27_MIG_energy.csv (see below)
└── gar/
    ├── utils.py         net_increase, block_boot_idx, trailing_mean
    ├── data.py          CSVs -> dict of tidy datasets (date-joined)
    ├── quantile_reg.py  qr_direct: direct-projection QR engine
    │                    (point/average targets, block bootstrap)
    ├── skewt.py         Azzalini skew-t (pdf/ppf), quantile matching,
    │                    GaR / ES / Longrise / entropy  (step2_match)
    ├── exceedance.py    Ang–Chen exceedance correlations + HTZ test
    └── plots.py         matplotlib figure builders
```

## Run

```sh
python validate_core.py    # ~1 min, must print 8/8 PASS
python run_all.py          # full pipeline -> output/
```

The skew-t fitting (stage 03) is the slow step: roughly 1 minute per
(spec × dataset × horizon × model). With the default
`horizons_density = [1, 3, 6, 12]` and no EU-energy file, expect ~1 h;
results are cached in `data/cache/` under collision-proof names, so re-runs
are fast. Reduce `n_boot` for quick iterations.

## Input files

`DataIP_GaR.csv` (EU27 aggregate), `DataOIL.csv` (Brent) and
`country_sector.csv` (DE/FR/IT × 8 NACE sectors) ship in `input/`. Add
`EU27_MIG_energy.csv` (Time, IP, GAS_PRICE for the EU27 MIG-energy grouping)
to enable the EU-level energy-sector results; without it those stages are
skipped with a message.

Growth rates and NGPI/NOPI are recomputed from raw levels; precomputed CSV
columns (including the corrupted `NGPI` column from the old MATLAB
`ReadData.m`) are ignored by design.

## Implementation notes

- Quantile regressions: `statsmodels` `QuantReg` (IRLS). Point estimates were
  verified to match Koenker's Frisch–Newton solver to ~1e-3 on this data.
- Skew-t: exact Azzalini–Capitanio density; the quantile function inverts a
  grid-based CDF (accurate to ~1e-3 standardized). Fitting profiles over
  integer ν ∈ 2…20 and exploits that (ξ, ω) is a *linear* LS problem given
  (α, ν) — only α is optimized numerically, which is faster and more robust
  than the MATLAB 3-parameter search.
- Inference: moving-block bootstrap (Künsch 1989), block max(12, h),
  bands ±1 SE (68%) and ±1.96 SE (95%).
