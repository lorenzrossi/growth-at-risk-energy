# 03_density_gar.R -- skew-t density fitting, GaR, Expected Shortfall /
# Longrise, and inverse-CDF comparisons (energy-augmented vs IP-only vs raw).
# Replaces: MainIPGas.m, MainIPNGPI.m, MainIPBrent.m, MainIPNOPI.m,
#           Step2match_improved.m, QuantilesInterpolation_improved.m
#
# Adrian-style averaged target: ybar(t+1..t+h) regressed on X(t).
# Note: results are cached per SPEC and HORIZON in data/cache with unique
# names (the MATLAB filename-collision bug cannot recur by construction).

source(file.path(cfg$root, "R", "utils.R"))
source(file.path(cfg$root, "R", "quantile_reg.R"))
source(file.path(cfg$root, "R", "skewt_fit.R"))
source(file.path(cfg$root, "R", "plots.R"))

dl <- readRDS(file.path(cfg$root, "data", "prepared.rds"))
fig_dir <- file.path(cfg$root, "output", "fig_density")
cache_dir <- file.path(cfg$root, "data", "cache")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

specs <- list(
  GAS   = list(var = "Diff_GP", label = "IP and Nat. Gas Price Change"),
  NGPI  = list(var = "NGPI",    label = "IP and NGPI"),
  BRENT = list(var = "Diff_BP", label = "IP and Brent Price Change"),
  NOPI  = list(var = "NOPI",    label = "IP and NOPI")
)
datasets <- c(AGG = "aggregate")
if ("EU27.energy" %in% names(dl)) datasets <- c(datasets, ENERGY = "EU27.energy")

# density-plot dates (COVID onset, Russian invasion)
focus_dates <- as.Date(c("2020-04-01", "2022-03-01"))

for (ds in names(datasets)) {
  d <- dl[[datasets[ds]]]
  for (sp_name in names(specs)) {
    sp <- specs[[sp_name]]
    if (!sp$var %in% names(d)) next
    ok <- stats::complete.cases(d$Diff_IP, d[[sp$var]])
    y <- d$Diff_IP[ok]; x <- d[[sp$var]][ok]; tm <- d$Time[ok]

    res_main <- qr_direct(y, matrix(x, dimnames = list(NULL, sp$var)),
                          horizons = cfg$horizons_density, taus = cfg$taus,
                          include_lag_y = TRUE, target = "average",
                          var_names = c("Intercept", sp$var, "Diff_IP_lag"))
    res_ip   <- qr_direct(y, NULL, horizons = cfg$horizons_density,
                          taus = cfg$taus, include_lag_y = TRUE,
                          target = "average")

    for (h in cfg$horizons_density) {
      cache <- file.path(cache_dir, sprintf("match_%s_%s_H%d.rds", sp_name, ds, h))
      if (file.exists(cache) && cfg$use_cache) {
        m <- readRDS(cache)
      } else {
        r  <- res_main[[as.character(h)]]
        ri <- res_ip[[as.character(h)]]
        yq_unc <- stats::quantile(trailing_mean(y, h), cfg$taus, na.rm = TRUE)
        m <- list(
          full   = step2_match(r$yq,  yq_unc, cfg$taus, verbose = FALSE),
          iponly = step2_match(ri$yq, yq_unc, cfg$taus, verbose = FALSE)
        )
        saveRDS(m, cache)
      }
      r <- res_main[[as.character(h)]]

      ## GaR / ES / Longrise over time
      save_pdf(plot_es_longrise(m$full, tm,
                 sprintf("ES and Longrise, %s (%s, h = %d)", sp_name, ds, h)),
               file.path(fig_dir, sprintf("ES_%s_%s_H%d.pdf", sp_name, ds, h)))

      ## fan chart
      save_pdf(plot_fan(r, tm, trailing_mean(y, h),
                 sprintf("Predicted distribution, %s (%s, h = %d)", sp_name, ds, h)),
               file.path(fig_dir, sprintf("Fan_%s_%s_H%d.pdf", sp_name, ds, h)))

      ## inverse CDF comparison at focus dates
      for (fi in seq_along(focus_dates)) {
        fd <- focus_dates[fi]         # keep Date class (for() would strip it)
        jt <- which(tm == fd)
        if (length(jt) != 1 || jt + h > length(tm)) next
        save_pdf(plot_inverse_cdf(m$full, m$iponly, r, jt, h,
                   labels = c(sp$label, "IP only", "Raw QR")),
                 file.path(fig_dir, sprintf("InvCDF_%s_%s_H%d_%s.pdf",
                           sp_name, ds, h, format(fd, "%Y_%m"))))
      }

      ## GaR series to CSV
      utils::write.csv(
        data.frame(Time = tm, GaR05 = m$full$gar,
                   ES = m$full$shortfall, Longrise = m$full$longrise,
                   skew = m$full$moments[, "skewness"]),
        file.path(fig_dir, sprintf("GaR_%s_%s_H%d.csv", sp_name, ds, h)),
        row.names = FALSE)

      message(sprintf("03: %s / %s h=%d done (crossed: %d)",
                      ds, sp_name, h, m$full$n_crossed))
    }
  }
}

message("03_density_gar complete.")
