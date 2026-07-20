# GaR project — R port

R implementation of "The effects of gas and oil prices on European industrial
production: a Growth-at-Risk approach". Replaces the ~11,000-line MATLAB
codebase with ~1,500 lines, using reference implementations for the numerical
core: **quantreg** (Koenker — same Frisch–Newton interior-point algorithm as
`rq.m`), **sn** (Azzalini — replaces the entire hand-ported skew-t suite),
plus **ggplot2**. All fixes from the code review are built in: MIG-energy data
for the sector analysis, one canonical NGPI/NOPI construction (100·log),
moving-block bootstrap, spec-specific caches, quantile-crossing check,
explicit date joins (no row-order alignment), single fixed seed.

## Setup

```r
install.packages(c("quantreg", "sn", "ggplot2"))
```

## Layout

```
R_port/
├── run_all.R            master script + configuration (EDIT PATHS HERE)
├── validate_core.R      numerical checks vs independently verified values
├── input/               input CSVs (aggregate + oil + country/sector included;
│                        add EU27_MIG_energy.csv — see below)
├── R/
│   ├── utils.R          net_increase, block_boot_idx, trailing_mean, save_pdf
│   ├── quantile_reg.R   qr_direct: direct-projection QR engine (point/average
│   │                    targets, block bootstrap) — replaces QRboot/run_static_qr
│   ├── skewt_fit.R      skew-t quantile matching, GaR/ES/Longrise/entropy
│   │                    — replaces Step2match/QuantilesInterpolation/azzalini
│   ├── exceedance.R     Ang–Chen exceedance correlations + HTZ test
│   └── plots.R          all figure builders (bands, scatters, fans, CDFs)
└── scripts/
    ├── 00_read_data.R   CSVs -> data/prepared.rds        (ReadData*, Country_Read_Data)
    ├── 01_exceedance.R  Figure 4 + n-per-tail + ex-2020  (ExceedenceCorrelation_Analysis)
    ├── 02_exploratory_qr.R  coefficient figures           (ExploratoryAndFit*)
    ├── 03_density_gar.R skew-t densities, GaR, ES        (MainIP*)
    ├── 04_subsample.R   pre-2020 / full / ex-2020        (SUB_*)
    └── 05_countries.R   Germany, Italy, France           (Germany/Italy/France.m)
```

## Run

```sh
cd R_port
Rscript validate_core.R   # ~1 min: checks the port reproduces known numbers
Rscript run_all.R         # full pipeline -> output/
```

Adjust in `run_all.R`: input paths, `horizons_density` (default c(1,3,6,12);
set 1:12 for the final run — the skew-t fitting is the slow step), `n_boot`.

## Input files

Three of the four inputs ship in `input/` (taken from the project data):

| File | Contents | MATLAB equivalent |
|---|---|---|
| `DataIP_GaR.csv` | EU27 aggregate IP + TTF gas | `DataIP_GaR.mat` |
| `DataOIL.csv` | Brent | `DataOIL.mat` |
| `country_sector.csv` | DE/FR/IT × 8 NACE sectors | `Germany/…`, `Italy/…`, `France/…` .mat files |
| `EU27_MIG_energy.csv` | **you must add this** (Time, IP, GAS_PRICE for the EU27 MIG-energy grouping) | `DataIP_GaR_MIG_-_energy_(except_section_E).mat` |

Without `EU27_MIG_energy.csv` everything runs except the EU-level energy-sector
results (scripts skip it with a message). Export it from Eurostat sts_inpr_m
(EU27, MIG_NRG) or convert the .mat: in MATLAB
`load('DataIP_GaR_MIG_-_energy_(except_section_E).mat'); writetable(array2table(X,'VariableNames',Mnem), 'EU27_MIG_energy_raw.csv')`
then keep Time, IP, GAS_PRICE columns.

Note `00_read_data.R` recomputes Diff_IP, Diff_GP and NGPI/NOPI from the raw
IP and price levels — precomputed columns in the CSVs (including the corrupted
`NGPI` column produced by the old `ReadData.m`) are ignored by design.

## Conventions (made explicit, per the review)

- **Targets.** `target = "point"` (y at t+h) for coefficient figures;
  `target = "average"` (mean of t+1…t+h, Adrian et al.) for densities/GaR.
  Both are supported by the same engine; state in the text which is used where.
- **Inference.** Moving-block bootstrap, block length max(12, h), n_boot
  replications; SEs are bootstrap SDs, bands ±1 SE (68%) and ±1.96 SE (95%).
- **NGPI/NOPI.** `net_increase(P)` on price levels, internally 100·log —
  everywhere, no exceptions.
- **Skew-t.** `fit_skewt_quantiles` matches the 7 quantiles
  (5/10/25/50/75/90/95) by L-BFGS-B over (ξ, ω, α), profiling over integer
  ν ∈ 2…20; crossed quantiles are rearranged and counted (`n_crossed`).
- **Caches.** `data/cache/match_<SPEC>_<DATASET>_H<h>.rds` — unique per
  spec/dataset/horizon; the MATLAB filename-collision bug cannot recur.

## Status / testing

Written and reviewed offline (R could not be installed in the sandbox used to
prepare this port — no root, CRAN blocked). The algorithms were validated
numerically against the project CSVs in Python during the code review;
`validate_core.R` re-runs those exact checks in R (correlations, NGPI
construction, aggregate Q05 NGPI coefficients at h = 1/12, Germany energy Q95
at h = 12, point-vs-average equivalence at h = 1). **Run `validate_core.R`
first**; if all checks pass, the numerical core is confirmed on your machine.
Expect to touch small things on first `run_all.R` (typos, a ggplot warning) —
the structure is sound but it has not executed end-to-end.
