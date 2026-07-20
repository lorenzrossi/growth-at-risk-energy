# plots.R -- ggplot2 figure builders.
# Replaces: PlotQRbands.m, PlotPredictiveTS.m and the per-script plotting code
# (coefficient term structures with 68/95% bands, QR scatter plots, fan
# charts, ES/Longrise, densities, inverse CDFs, exceedance profiles).

library(ggplot2)

gar_theme <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.border = element_rect(colour = "grey30", fill = NA),
          legend.position = "bottom", plot.title = element_text(size = 12))
}

#' Coefficient across quantiles at a fixed horizon, with 68/95% bands.
plot_coef_by_quantile <- function(res, var, h, var_label = var) {
  r <- res[[as.character(h)]]
  df <- data.frame(tau = r$taus, coef = r$coef[var, ], se = r$se[var, ])
  ggplot(df, aes(tau, coef)) +
    geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se),
                fill = "grey75", alpha = .5) +
    geom_ribbon(aes(ymin = coef - se, ymax = coef + se),
                fill = "grey55", alpha = .5) +
    geom_line(colour = "red", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = expression(tau), y = paste(var_label, "coefficient"),
         title = sprintf("%s coefficients across quantiles (h = %d)", var_label, h)) +
    gar_theme()
}

#' Coefficient across horizons at a fixed quantile, with 68/95% bands.
plot_coef_by_horizon <- function(res, var, tau, var_label = var) {
  tab <- coef_table(res, var)
  df <- tab[abs(tab$tau - tau) < 1e-9, ]
  ggplot(df, aes(h, coef)) +
    geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se),
                fill = "grey75", alpha = .5) +
    geom_ribbon(aes(ymin = coef - se, ymax = coef + se),
                fill = "grey55", alpha = .5) +
    geom_line(colour = "red", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = unique(df$h)) +
    labs(x = "Horizon (months)", y = paste(var_label, "coefficient"),
         title = sprintf("%s coefficients (Q %.2f), 68%% and 95%% block-bootstrap bands",
                         var_label, tau)) +
    gar_theme()
}

#' Scatter of target vs one regressor with fitted quantile lines + OLS.
#' Lines connect the fitted values at the endpoints of the sorted regressor,
#' exactly as in the MATLAB scripts (other regressors' contributions are
#' embedded in the fitted values).
plot_qr_scatter <- function(res, y, x, h, x_label = "X", y_label = "IP growth",
                            taus_show = c(0.05, 0.5, 0.95)) {
  r <- res[[as.character(h)]]
  n <- length(y)
  yv <- if (r$target == "point") y else trailing_mean(y, h)
  df <- data.frame(x = x[1:(n - h)],
                   y = yv[(h + 1):n],
                   t = (h + 1):n)
  df <- df[stats::complete.cases(df$x, df$y), ]

  i_lo <- which.min(df$x); i_hi <- which.max(df$x)
  lines <- do.call(rbind, lapply(taus_show, function(tq) {
    j <- which.min(abs(r$taus - tq))
    yf <- r$yq[df$t, j]
    data.frame(x = df$x[c(i_lo, i_hi)], yend = yf[c(i_lo, i_hi)],
               tau = sprintf("Q%02d", round(100 * tq)))
  }))

  ggplot(df, aes(x, y)) +
    geom_point(colour = "steelblue", alpha = .6) +
    geom_line(data = lines, aes(x, yend, colour = tau), linewidth = .9) +
    labs(x = x_label,
         y = sprintf("%s %d month(s) ahead", y_label, h),
         colour = NULL,
         title = sprintf("Quantile regressions (h = %d)", h)) +
    gar_theme()
}

#' Predictive quantile fan chart with realised series.
plot_fan <- function(r, time, y_real, title = "Predicted distribution") {
  taus <- r$taus
  bands <- data.frame(
    Time = time,
    q05 = r$yq[, which.min(abs(taus - .05))],
    q25 = r$yq[, which.min(abs(taus - .25))],
    q50 = r$yq[, which.min(abs(taus - .50))],
    q75 = r$yq[, which.min(abs(taus - .75))],
    q95 = r$yq[, which.min(abs(taus - .95))],
    real = y_real
  )
  ggplot(bands, aes(Time)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = "grey80") +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey60") +
    geom_line(aes(y = q50), colour = "black") +
    geom_line(aes(y = real), colour = "red", linewidth = .4) +
    labs(y = "IPI growth", title = title) + gar_theme()
}

#' Expected shortfall / longrise over time.
plot_es_longrise <- function(match, time, title = "Expected Shortfall and Longrise") {
  df <- rbind(
    data.frame(Time = time, value = match$shortfall, measure = "Shortfall"),
    data.frame(Time = time, value = match$longrise,  measure = "Longrise")
  )
  ggplot(df, aes(Time, value, linetype = measure)) +
    geom_line(linewidth = .7, na.rm = TRUE) +
    labs(y = NULL, linetype = NULL, title = title) + gar_theme()
}

#' Inverse CDF (quantile function) comparison for one date.
plot_inverse_cdf <- function(match_full, match_iponly, r_full, date_index, h,
                             labels = c("Energy-augmented", "IP only", "Raw QR")) {
  taus <- match_full$taus
  jt <- date_index + h
  df <- rbind(
    data.frame(tau = taus, q = match_full$QST[jt, ],   model = labels[1]),
    data.frame(tau = taus, q = match_iponly$QST[jt, ], model = labels[2]),
    data.frame(tau = taus, q = r_full$yq[jt, ],        model = labels[3])
  )
  ggplot(df, aes(tau, q, colour = model, linetype = model)) +
    geom_line(linewidth = .9) +
    labs(x = expression(tau), y = "IPI growth",
         title = sprintf("Fitted inverse CDF (h = %d)", h)) +
    gar_theme()
}

#' Exceedance-correlation profile (optionally two samples overlaid).
plot_exceedance <- function(tab_full, tab_ex = NULL, title = "Exceedance correlations") {
  df <- rbind(
    data.frame(threshold = tab_full$threshold, corr = tab_full$lower,
               tail = "Lower", sample = "Full"),
    data.frame(threshold = tab_full$threshold, corr = tab_full$upper,
               tail = "Upper", sample = "Full")
  )
  if (!is.null(tab_ex)) {
    df <- rbind(df,
      data.frame(threshold = tab_ex$threshold, corr = tab_ex$lower,
                 tail = "Lower", sample = "Ex-2020"),
      data.frame(threshold = tab_ex$threshold, corr = tab_ex$upper,
                 tail = "Upper", sample = "Ex-2020"))
  }
  ggplot(df, aes(threshold, corr, colour = tail, linetype = sample)) +
    geom_line(linewidth = .9) + geom_point() +
    labs(x = "Quantile threshold", y = "Exceedance correlation", title = title) +
    gar_theme()
}
