# 05_countries.R -- Germany, Italy, France: manufacturing + energy sector.
# Replaces: Germany.m, Italy.m, France.m

source(file.path(cfg$root, "R", "utils.R"))
source(file.path(cfg$root, "R", "quantile_reg.R"))
source(file.path(cfg$root, "R", "plots.R"))

dl <- readRDS(file.path(cfg$root, "data", "prepared.rds"))
fig_dir <- file.path(cfg$root, "output", "fig_countries")
tab_dir <- file.path(cfg$root, "output", "tab_countries")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

countries <- c("Germany", "Italy", "France")
sectors <- c(Manufacturing = "Manufacturing",
             Energy = "MIG___energy__except_section_E_")

specs <- list(
  list(var = "Diff_GP", label = "Nat. Gas Price Change"),
  list(var = "NGPI",    label = "NGPI"),
  list(var = "Diff_BP", label = "Brent Price Change"),
  list(var = "NOPI",    label = "NOPI")
)

for (geo in countries) {
  for (sec in names(sectors)) {
    # find the dataset key produced by 00_read_data (tolerant match)
    key <- grep(paste0("^", geo, "\\."), names(dl), value = TRUE)
    key <- key[grepl(if (sec == "Energy") "energy" else "Manufacturing",
                     key, ignore.case = TRUE)]
    if (length(key) == 0) { message("skip ", geo, "/", sec); next }
    d <- dl[[key[1]]]

    for (sp in specs) {
      if (!sp$var %in% names(d)) next
      ok <- stats::complete.cases(d$Diff_IP, d[[sp$var]])
      y <- d$Diff_IP[ok]; x <- d[[sp$var]][ok]

      res <- qr_direct(y, matrix(x, dimnames = list(NULL, sp$var)),
                       horizons = cfg$horizons, taus = cfg$taus,
                       include_lag_y = TRUE, target = "point",
                       n_boot = cfg$n_boot,
                       var_names = c("Intercept", sp$var, "Diff_IP_lag"))

      tab <- add_stars(coef_table(res, sp$var))
      utils::write.csv(tab, file.path(tab_dir,
        sprintf("%s_%s_%s_coefs.csv", geo, sec, sp$var)), row.names = FALSE)

      for (tau in c(0.05, 0.95)) {
        save_pdf(plot_coef_by_horizon(res, sp$var, tau,
                   paste(geo, sec, "-", sp$label)),
                 file.path(fig_dir, sprintf("%s_%s_%s_Q%02d.pdf",
                   geo, sec, sp$var, round(100 * tau))))
      }
      message(sprintf("05: %s / %s / %-8s done", geo, sec, sp$var))
    }
  }
}

message("05_countries complete.")
