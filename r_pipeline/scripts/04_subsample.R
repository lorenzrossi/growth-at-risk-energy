# 04_subsample.R -- subsample comparison: pre-2020 vs full vs ex-2020.
# Replaces: SUB_analysis.m, SUB_NGPI.m, SUB_BP.m (and adds the ex-2020
# robustness recommended in the review: full sample minus calendar 2020).

source(file.path(cfg$root, "R", "utils.R"))
source(file.path(cfg$root, "R", "quantile_reg.R"))
source(file.path(cfg$root, "R", "plots.R"))

dl <- readRDS(file.path(cfg$root, "data", "prepared.rds"))
fig_dir <- file.path(cfg$root, "output", "fig_subsample")
tab_dir <- file.path(cfg$root, "output", "tab_subsample")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

specs <- list(
  list(var = "Diff_GP", label = "Nat. Gas Price Change", tau = 0.05),
  list(var = "NGPI",    label = "NGPI",                  tau = 0.05),
  list(var = "Diff_BP", label = "Brent Price Change",    tau = 0.05)
)
datasets <- c(aggregate = "aggregate")
if ("EU27.energy" %in% names(dl)) datasets <- c(datasets, energy = "EU27.energy")

fit_sample <- function(d, var, mask) {
  ok <- stats::complete.cases(d$Diff_IP, d[[var]]) & mask
  y <- d$Diff_IP[ok]; x <- d[[var]][ok]
  qr_direct(y, matrix(x, dimnames = list(NULL, var)),
            horizons = cfg$horizons, taus = cfg$taus,
            include_lag_y = TRUE, target = "point", n_boot = cfg$n_boot,
            var_names = c("Intercept", var, "Diff_IP_lag"))
}

for (ds in names(datasets)) {
  d <- dl[[datasets[ds]]]
  yr <- as.integer(format(d$Time, "%Y"))
  masks <- list(full = rep(TRUE, nrow(d)),
                pre2020 = yr <= 2019,
                ex2020 = yr != 2020)
  # for the energy sector the paper focuses on the upper tail
  taus_show <- if (ds == "energy") c(0.95) else c(0.05)

  for (sp in specs) {
    if (!sp$var %in% names(d)) next
    tabs <- list()
    for (mn in names(masks)) {
      res <- fit_sample(d, sp$var, masks[[mn]])
      tt <- coef_table(res, sp$var); tt$sample <- mn
      tabs[[mn]] <- tt
    }
    tab <- do.call(rbind, tabs)
    utils::write.csv(tab, file.path(tab_dir,
                     sprintf("%s_%s_subsamples.csv", ds, sp$var)), row.names = FALSE)

    for (tau in taus_show) {
      df <- tab[abs(tab$tau - tau) < 1e-9, ]
      p <- ggplot(df, aes(h, coef, colour = sample)) +
        geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se,
                        fill = sample), alpha = .15, colour = NA) +
        geom_line(linewidth = 1) +
        geom_hline(yintercept = 0, linetype = 2) +
        scale_x_continuous(breaks = unique(df$h)) +
        labs(x = "Horizon (months)", y = paste(sp$label, "coefficient"),
             title = sprintf("%s, Q %.2f: full vs pre-2020 vs ex-2020 (%s)",
                             sp$label, tau, ds)) +
        gar_theme()
      save_pdf(p, file.path(fig_dir, sprintf("%s_%s_Q%02d_subsamples.pdf",
                                             ds, sp$var, round(100 * tau))))
    }
    message(sprintf("04: %s / %s done", ds, sp$var))
  }
}

message("04_subsample complete.")
