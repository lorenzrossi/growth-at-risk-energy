# growth-at-risk-energy

Replication code for **"The effects of gas and oil prices on European
industrial production: a Growth-at-Risk approach"** (Lorenzo Rossi,
University of Milan).

The paper applies the Growth-at-Risk framework of Adrian, Boyarchenko &
Giannone (2019, *AER*) to European industrial production, conditioning on
gas and oil price dynamics (price changes and Hamilton-style net price
increases, NGPI/NOPI), at the EU27 aggregate, sector (MIG energy), and
country level (Germany, Italy, France), over 2004–2023 monthly data.

## Repository layout

```
r_pipeline/        R implementation      (quantreg + sn + ggplot2)
python_pipeline/   Python implementation (statsmodels + scipy + matplotlib)
```

The two pipelines are functionally equivalent and share the same input CSVs
(bundled in each `input/` folder) and the same conventions:

- one canonical NGPI/NOPI construction (100·log scale, `net_increase`);
- direct-projection quantile regressions at horizons 1–12, τ = 0.05…0.95,
  with either the single-month target (coefficient figures) or the
  Adrian-et-al. h-month-average target (density/GaR results);
- moving-block bootstrap inference (Künsch 1989), block length max(12, h);
- skew-t (Azzalini–Capitanio) quantile matching profiling over integer ν,
  with quantile-crossing rearrangement;
- exceedance correlations (Ang & Chen 2002) with the Hong–Tu–Zhou test,
  tail sample sizes, and ex-2020 robustness;
- fixed seed (123) end to end.

## Quick start

Python (validated end-to-end):

```sh
cd python_pipeline
pip install -r requirements.txt
python validate_core.py     # must print 8/8 PASS
python run_all.py
```

R:

```sh
cd r_pipeline
Rscript -e 'install.packages(c("quantreg","sn","ggplot2"))'
Rscript validate_core.R
Rscript run_all.R
```

Each pipeline README documents inputs, configuration, and the one missing
optional input (`EU27_MIG_energy.csv`, the EU-level MIG-energy sector series)
needed for the EU energy-sector results.

## Data

Input CSVs are derived from Eurostat (industrial production, sts_inpr_m),
TTF natural gas and Brent oil price series. If you make this repository
public, check the redistribution terms of the underlying data providers.
