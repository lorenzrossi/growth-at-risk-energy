# quantile_reg.R -- direct-projection quantile regressions with block-bootstrap
# inference. Replaces: QRboot.m, run_static_qr.m, rq.m (via quantreg::rq,
# Koenker's own Frisch-Newton implementation, method = "fn").

#' Direct-projection quantile (and OLS) regressions across horizons.
#'
#' For each horizon h, regresses the target at t+h on conditioning
#' variables at t. Two target conventions:
#'   target = "point"   : y(t+h)                       (exploratory/coefficient figures)
#'   target = "average" : mean(y(t+1), ..., y(t+h))    (Adrian et al. 2019; density/GaR)
#'
#' @param y        numeric vector, dependent variable (monthly growth)
#' @param X        NULL, vector, or matrix of conditioning variables (aligned with y)
#' @param horizons integer vector of horizons
#' @param taus     quantile levels
#' @param include_lag_y include y(t) as a regressor (default TRUE)
#' @param target   "point" or "average"
#' @param n_boot   bootstrap replications (0 = none)
#' @param block_length NULL -> max(12, h) per horizon
#' @param var_names names for the X columns (for output labelling)
#' @return list with one element per horizon: coef (k x ntau), se (k x ntau),
#'   ols (k), ols_se (k), yq (T x ntau fitted quantiles, aligned so row t is
#'   the fitted quantile of the target dated t; first h rows NA), boot
#'   (n_boot x k x ntau array), n, taus, var_names
qr_direct <- function(y, X = NULL, horizons = 1:12,
                      taus = seq(0.05, 0.95, by = 0.05),
                      include_lag_y = TRUE,
                      target = c("point", "average"),
                      n_boot = 0, block_length = NULL,
                      var_names = NULL) {
  target <- match.arg(target)
  stopifnot(requireNamespace("quantreg", quietly = TRUE))
  n <- length(y)

  if (!is.null(X)) X <- as.matrix(X)
  Z <- cbind(Intercept = rep(1, n))
  if (!is.null(X)) Z <- cbind(Z, X)
  if (include_lag_y) Z <- cbind(Z, y_lag = y)

  if (is.null(var_names)) {
    xn <- if (is.null(X)) character(0) else
      if (!is.null(colnames(X))) colnames(X) else paste0("X", seq_len(ncol(X)))
    var_names <- c("Intercept", xn, if (include_lag_y) "y_lag")
  }
  colnames(Z) <- var_names
  k <- ncol(Z)

  make_target <- function(yv, h) {
    if (target == "point") {
      yv                                   # will be shifted below
    } else {
      trailing_mean(yv, h)                 # ybar[t] = mean(y[(t-h+1):t])
    }
  }

  fit_one <- function(yh, Zh, tau) {
    # Koenker's Frisch-Newton interior point ("fn"), same algorithm as rq.m
    fit <- quantreg::rq.fit.fnb(Zh, yh, tau = tau)
    fit$coefficients
  }

  out <- list()
  for (h in horizons) {
    yt <- make_target(y, h)
    # align: target dated t+h regressed on Z at t
    yh <- yt[(h + 1):n]
    Zh <- Z[1:(n - h), , drop = FALSE]
    ok <- stats::complete.cases(yh, Zh)
    yh <- yh[ok]; Zh <- Zh[ok, , drop = FALSE]
    nh <- length(yh)

    coef <- matrix(NA_real_, k, length(taus),
                   dimnames = list(var_names, paste0("tau", taus)))
    for (j in seq_along(taus)) coef[, j] <- fit_one(yh, Zh, taus[j])

    ols <- stats::lm.fit(Zh, yh)$coefficients

    # fitted quantiles aligned to target date (row t = fitted for date t)
    yq <- matrix(NA_real_, n, length(taus))
    target_rows <- ((h + 1):n)[ok]
    yq[target_rows, ] <- Zh %*% coef

    # moving-block bootstrap
    se <- matrix(NA_real_, k, length(taus), dimnames = dimnames(coef))
    ols_se <- rep(NA_real_, k)
    boot <- NULL
    if (n_boot > 0) {
      bl <- if (is.null(block_length)) max(12L, h) else block_length
      boot <- array(NA_real_, c(n_boot, k, length(taus)))
      boot_ols <- matrix(NA_real_, n_boot, k)
      for (b in seq_len(n_boot)) {
        idx <- block_boot_idx(nh, bl)
        Zb <- Zh[idx, , drop = FALSE]; yb <- yh[idx]
        for (j in seq_along(taus)) {
          cb <- tryCatch(fit_one(yb, Zb, taus[j]), error = function(e) rep(NA_real_, k))
          boot[b, , j] <- cb
        }
        boot_ols[b, ] <- stats::lm.fit(Zb, yb)$coefficients
      }
      se[] <- apply(boot, c(2, 3), stats::sd, na.rm = TRUE)
      ols_se <- apply(boot_ols, 2, stats::sd, na.rm = TRUE)
    }

    out[[as.character(h)]] <- list(
      h = h, coef = coef, se = se, ols = ols, ols_se = ols_se,
      yq = yq, boot = boot, n = nh, taus = taus, var_names = var_names,
      target = target
    )
  }
  class(out) <- "qr_direct"
  out
}

#' Tidy coefficient table across horizons/quantiles for one regressor.
#' @param res result of qr_direct
#' @param var regressor name (e.g. "NGPI")
coef_table <- function(res, var) {
  do.call(rbind, lapply(res, function(r) {
    data.frame(
      h = r$h, tau = r$taus,
      coef = r$coef[var, ],
      se = r$se[var, ],
      row.names = NULL
    )
  }))
}

#' Significance stars from a normal approximation to the bootstrap
add_stars <- function(tab) {
  z <- abs(tab$coef / tab$se)
  tab$pval <- 2 * stats::pnorm(-z)
  tab$sig <- cut(tab$pval, c(-Inf, .01, .05, .10, Inf),
                 labels = c("***", "**", "*", ""))
  tab
}
