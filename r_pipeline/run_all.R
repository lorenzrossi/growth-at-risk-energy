# run_all.R -- master script for the GaR project (R port).
# Runs the full pipeline from raw CSVs to figures with a fixed seed.
#
# Usage:  Rscript run_all.R          (from the R_port folder)
#    or:  source("run_all.R")        (from an R session)

## ---- configuration ---------------------------------------------------------
cfg <- list(
  root = normalizePath("."),

  # --- input files: EDIT THESE PATHS -----------------------------------------
  files = list(
    aggregate = "input/DataIP_GaR.csv",      # EU27: Time, IP, GAS_PRICE, ...
    brent     = "input/DataOIL.csv",         # Time, BRENT_PRICE, ...
    panel     = "input/country_sector.csv",  # Time, nace_r2, geo, IP, GAS_PRICE, BRENT_PRICE
    eu_sectors = list(                       # optional EU27 sector CSVs
      energy = "input/EU27_MIG_energy.csv"   # Time, IP, GAS_PRICE
    )
  ),

  # --- estimation settings ---------------------------------------------------
  horizons = 1:12,          # coefficient figures
  horizons_density = c(1, 3, 6, 12),  # skew-t fitting (slow) -- extend to 1:12 for final run
  taus = seq(0.05, 0.95, by = 0.05),
  n_boot = 1000,            # block-bootstrap replications
  use_cache = TRUE,         # cache skew-t fits (spec-specific names, no collisions)
  seed = 123
)

## ---- dependencies ----------------------------------------------------------
required <- c("quantreg", "sn", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0)
  stop("Install missing packages first: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")

## ---- pipeline --------------------------------------------------------------
set.seed(cfg$seed)
message("=== 00: data construction ===")
source(file.path(cfg$root, "scripts", "00_read_data.R"), local = FALSE)

set.seed(cfg$seed)
message("=== 01: exceedance correlations (motivation) ===")
source(file.path(cfg$root, "scripts", "01_exceedance.R"), local = FALSE)

set.seed(cfg$seed)
message("=== 02: exploratory quantile regressions ===")
source(file.path(cfg$root, "scripts", "02_exploratory_qr.R"), local = FALSE)

set.seed(cfg$seed)
message("=== 03: skew-t densities / GaR ===")
source(file.path(cfg$root, "scripts", "03_density_gar.R"), local = FALSE)

set.seed(cfg$seed)
message("=== 04: subsample comparison ===")
source(file.path(cfg$root, "scripts", "04_subsample.R"), local = FALSE)

set.seed(cfg$seed)
message("=== 05: country analysis ===")
source(file.path(cfg$root, "scripts", "05_countries.R"), local = FALSE)

message("=== run_all complete: see output/ ===")
