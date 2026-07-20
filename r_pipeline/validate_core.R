# validate_core.R -- quick numerical validation of the R port against
# independently verified reference values (computed from the same CSVs with
# Python/statsmodels during the code review; those in turn matched the
# magnitudes reported in the MATLAB-based draft).
#
# Run AFTER 00_read_data.R:   Rscript validate_core.R
# Every check prints the R value, the reference, and PASS/FAIL (tolerances
# are loose for the bootstrap-free point estimates; QR solvers may differ
# in the last decimals when the solution is non-unique).

cfg <- list(root = normalizePath("."),
            files = list(
              aggregate = "input/DataIP_GaR.csv",
              brent     = "input/DataOIL.csv",
              panel     = "input/country_sector.csv",
              eu_sectors = list()),
            horizons = c(1, 12),
            taus = seq(0.05, 0.95, by = 0.05),
            n_boot = 0, use_cache = FALSE, seed = 123)

set.seed(cfg$seed)
source(file.path(cfg$root, "scripts", "00_read_data.R"))
source(file.path(cfg$root, "R", "quantile_reg.R"))

dl <- readRDS(file.path(cfg$root, "data", "prepared.rds"))

check <- function(name, value, ref, tol) {
  ok <- is.finite(value) && abs(value - ref) <= tol
  cat(sprintf("%-55s R = %+.4f | ref = %+.4f | %s\n",
              name, value, ref, if (ok) "PASS" else "FAIL"))
  ok
}

results <- c()

## 1. Unconditional correlations --------------------------------------------
d <- dl$aggregate
ok <- stats::complete.cases(d$Diff_IP, d$Diff_GP, d$Diff_BP)
results <- c(results,
  check("corr(Diff_IP, Diff_GP), 2004-2023",
        stats::cor(d$Diff_IP[ok], d$Diff_GP[ok]), 0.077, 0.02),
  check("corr(Diff_IP, Diff_BP), 2004-2023",
        stats::cor(d$Diff_IP[ok], d$Diff_BP[ok]), 0.499, 0.02))

## 2. NGPI construction -------------------------------------------------------
results <- c(results,
  check("max(NGPI) aggregate (100*log scale)",
        max(d$NGPI, na.rm = TRUE), 37.49, 0.5))

## 3. Aggregate Q5 NGPI coefficient at h = 1 and h = 12 ----------------------
okn <- stats::complete.cases(d$Diff_IP, d$NGPI)
y <- d$Diff_IP[okn]; x <- d$NGPI[okn]
res <- qr_direct(y, matrix(x, dimnames = list(NULL, "NGPI")),
                 horizons = c(1, 12), taus = cfg$taus,
                 include_lag_y = TRUE, target = "point",
                 var_names = c("Intercept", "NGPI", "Diff_IP_lag"))
b_h1  <- res[["1"]]$coef["NGPI", which.min(abs(cfg$taus - 0.05))]
b_h12 <- res[["12"]]$coef["NGPI", which.min(abs(cfg$taus - 0.05))]
results <- c(results,
  check("Aggregate Q05 NGPI coefficient, h = 1",  b_h1,  -0.073, 0.03),
  check("Aggregate Q05 NGPI coefficient, h = 12", b_h12, -0.196, 0.04))

## 4. Germany energy sector, Q95 NGPI at h = 12 ------------------------------
gkey <- grep("^Germany\\..*energy", names(dl), ignore.case = TRUE, value = TRUE)
if (length(gkey) == 1) {
  g <- dl[[gkey]]
  okg <- stats::complete.cases(g$Diff_IP, g$NGPI)
  resg <- qr_direct(g$Diff_IP[okg],
                    matrix(g$NGPI[okg], dimnames = list(NULL, "NGPI")),
                    horizons = 12, taus = cfg$taus,
                    include_lag_y = TRUE, target = "point",
                    var_names = c("Intercept", "NGPI", "Diff_IP_lag"))
  bg <- resg[["12"]]$coef["NGPI", which.min(abs(cfg$taus - 0.95))]
  results <- c(results,
    check("Germany energy Q95 NGPI coefficient, h = 12", bg, -0.145, 0.04))
} else {
  cat("Germany energy dataset not found -- skipped check 4\n")
}

## 5. Averaged-target sanity: h = 1 equals point target ----------------------
res_avg <- qr_direct(y, matrix(x, dimnames = list(NULL, "NGPI")),
                     horizons = 1, taus = cfg$taus,
                     include_lag_y = TRUE, target = "average",
                     var_names = c("Intercept", "NGPI", "Diff_IP_lag"))
results <- c(results,
  check("Averaged target == point target at h = 1 (Q05 NGPI)",
        res_avg[["1"]]$coef["NGPI", which.min(abs(cfg$taus - 0.05))], b_h1, 1e-6))

cat(sprintf("\n%d of %d checks passed.\n", sum(results), length(results)))
if (!all(results)) stop("Validation failed -- inspect the FAIL lines above.")
