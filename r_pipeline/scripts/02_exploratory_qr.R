# 02_exploratory_qr.R -- coefficient term-structure figures and QR scatters
# for all four energy variables, aggregate + EU energy sector.
# Replaces: ExploratoryAndFit.m, ExploratoryAndFit_NGPI.m,
#           ExploratoryAndFitBP.m, ExploratoryAndFitNOPI.m, run_static_qr.m
#
# Single-month target (y(t+h)); moving-block bootstrap, block = max(12, h).

source(file.path(cfg$root, "R", "utils.R"))
source(file.path(cfg$root, "R", "quantile_reg.R"))
source(file.path(cfg$root, "R", "plots.R"))

dl <- readRDS(file.path(cfg$root, "data", "prepared.rds"))
fig_dir <- file.path(cfg$root, "output", "fig_exploratory")
tab_dir <- file.path(cfg$root, "output", "tab_exploratory")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

specs <- list(
  list(var = "Diff_GP", label = "Nat. Gas Price Change"),
  list(var = "NGPI",    label = "NGPI"),
  list(var = "Diff_BP", label = "Brent Price Change"),
  list(var = "NOPI",    label = "NOPI")
)

datasets <- c(aggregate = "aggregate")
if ("EU27.energy" %in% names(dl)) datasets <- c(datasets, energy = "EU27.energy")

for (ds in names(datasets)) {
  d <- dl[[datasets[ds]]]
  for (sp in specs) {
    if (!sp$var %in% names(d)) { message("skip ", ds, "/", sp$var); next }
    ok <- stats::complete.cases(d$Diff_IP, d[[sp$var]])
    y <- d$Diff_IP[ok]; x <- d[[sp$var]][ok]

    res <- qr_direct(y, X = matrix(x, dimnames = list(NULL, sp$var)),
                     horizons = cfg$horizons, taus = cfg$taus,
                     include_lag_y = TRUE, target = "point",
                     n_boot = cfg$n_boot,
                     var_names = c("Intercept", sp$var, "Diff_IP_lag"))

    tab <- add_stars(coef_table(res, sp$var))
    utils::write.csv(tab, file.path(tab_dir,
                     sprintf("%s_%s_coefs.csv", ds, sp$var)), row.names = FALSE)

    for (tau in c(0.05, 0.95)) {
      save_pdf(plot_coef_by_horizon(res, sp$var, tau, sp$label),
               file.path(fig_dir, sprintf("%s_%s_Q%02d_vs_horizon.pdf",
                                          ds, sp$var, round(100 * tau))))
    }
    for (h in intersect(c(1, 3, 6, 12), cfg$horizons)) {
      save_pdf(plot_coef_by_quantile(res, sp$var, h, sp$label),
               file.path(fig_dir, sprintf("%s_%s_coef_by_quantile_H%d.pdf",
                                          ds, sp$var, h)))
      save_pdf(plot_qr_scatter(res, y, x, h, x_label = sp$label),
               file.path(fig_dir, sprintf("%s_%s_scatter_H%d.pdf",
                                          ds, sp$var, h)))
    }
    message(sprintf("02: %s / %-8s done (n=%d)", ds, sp$var, res[["1"]]$n))
  }
}

message("02_exploratory_qr complete.")
