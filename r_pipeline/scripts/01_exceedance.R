# 01_exceedance.R -- Figure 4 (exceedance correlations) with n-per-tail
# and ex-2020 robustness, aggregate + countries.
# Replaces: ExceedenceCorrelation_Analysis.m + ExceedanceRobustness.m

source(file.path(cfg$root, "R", "utils.R"))
source(file.path(cfg$root, "R", "exceedance.R"))
source(file.path(cfg$root, "R", "plots.R"))

dl <- readRDS(file.path(cfg$root, "data", "prepared.rds"))
fig_dir <- file.path(cfg$root, "output", "fig_exceedance")
tab_dir <- file.path(cfg$root, "output", "tab_exceedance")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

run_pair <- function(d, ycol, name, label) {
  ok <- stats::complete.cases(d$Diff_IP, d[[ycol]])
  X <- d$Diff_IP[ok]; Y <- d[[ycol]][ok]; tm <- d$Time[ok]
  ex <- format(tm, "%Y") != "2020"

  full <- exceedance_correl(X, Y)
  ex20 <- exceedance_correl(X[ex], Y[ex])

  utils::write.csv(full$table, file.path(tab_dir, paste0(name, "_full.csv")),
                   row.names = FALSE)
  utils::write.csv(ex20$table, file.path(tab_dir, paste0(name, "_ex2020.csv")),
                   row.names = FALSE)

  p <- plot_exceedance(full$table, ex20$table,
        title = sprintf("%s  (HTZ asymmetry test p = %.3f)", label, full$pval))
  save_pdf(p, file.path(fig_dir, paste0(name, ".pdf")))

  j90 <- which.min(abs(full$table$threshold - 0.90))  # fp-safe lookup
  message(sprintf("%-28s HTZ p=%.3f | n at q=0.90: lower=%d upper=%d",
                  name, full$pval,
                  full$table$n_lower[j90], full$table$n_upper[j90]))
}

## Aggregate
run_pair(dl$aggregate, "Diff_GP", "EU27_IP_vs_Gas", "EU27 IP vs Gas price changes")
run_pair(dl$aggregate, "Diff_BP", "EU27_IP_vs_Oil", "EU27 IP vs Brent price changes")

## Countries (manufacturing series as in the MATLAB code)
for (geo in c("Germany", "France", "Italy")) {
  key <- paste0(geo, ".Manufacturing")
  if (!key %in% names(dl)) next
  run_pair(dl[[key]], "Diff_GP", paste0(geo, "_IP_vs_Gas"),
           paste(geo, "IP vs Gas price changes"))
  run_pair(dl[[key]], "Diff_BP", paste0(geo, "_IP_vs_Oil"),
           paste(geo, "IP vs Brent price changes"))
}

message("01_exceedance complete.")
