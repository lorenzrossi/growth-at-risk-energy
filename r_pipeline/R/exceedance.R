# exceedance.R -- exceedance correlations (Ang & Chen 2002) with the
# Hong-Tu-Zhou (2003) asymmetry test.
# Replaces: exceedence_correl.m (fixed version), corrcoef12.m, newey_west.m,
# plot_exceedence_correlation.m.

#' Newey-West long-run covariance of a T x k matrix of (demeaned) scores.
newey_west_cov <- function(Xi, lags = NULL) {
  Tn <- nrow(Xi)
  if (is.null(lags)) lags <- floor(4 * (Tn / 100)^(2 / 9))
  Xi <- scale(Xi, center = TRUE, scale = FALSE)
  om <- crossprod(Xi) / Tn
  if (lags > 0) {
    for (l in seq_len(lags)) {
      w <- 1 - l / (lags + 1)
      G <- crossprod(Xi[(l + 1):Tn, , drop = FALSE], Xi[1:(Tn - l), , drop = FALSE]) / Tn
      om <- om + w * (G + t(G))
    }
  }
  om
}

#' Exceedance correlations across quantile thresholds, with tail sample
#' sizes and the HTZ symmetry test.
#'
#' Lower tail: corr(X, Y | X <= Qx(1-q), Y <= Qy(1-q))
#' Upper tail: corr(X, Y | X >= Qx(q),   Y >= Qy(q)),  q in [0.5, 0.95]
#'
#' @param X,Y numeric vectors
#' @param qc  vector of upper-half thresholds (>= 0.5)
#' @param min_n minimum tail observations to compute a correlation
#' @return list(table = data.frame(threshold, lower, upper, n_lower, n_upper),
#'              teststat, pval)
exceedance_correl <- function(X, Y, qc = seq(0.5, 0.95, by = 0.05), min_n = 3L) {
  ok <- stats::complete.cases(X, Y)
  X <- X[ok]; Y <- Y[ok]
  Tn <- length(X)

  res <- data.frame(threshold = qc, lower = NA_real_, upper = NA_real_,
                    n_lower = NA_integer_, n_upper = NA_integer_)
  Xi <- NULL; rhodiff <- c(); used <- c()

  for (j in seq_along(qc)) {
    q <- qc[j]
    lo <- which(X <= stats::quantile(X, 1 - q) & Y <= stats::quantile(Y, 1 - q))
    up <- which(X >= stats::quantile(X, q)     & Y >= stats::quantile(Y, q))
    res$n_lower[j] <- length(lo); res$n_upper[j] <- length(up)
    if (length(lo) < min_n || length(up) < min_n) next

    r_lo <- stats::cor(X[lo], Y[lo]); r_up <- stats::cor(X[up], Y[up])
    res$lower[j] <- r_lo; res$upper[j] <- r_up

    # HTZ influence scores
    zx_up <- (X - mean(X[up])) / stats::sd(X[up])
    zy_up <- (Y - mean(Y[up])) / stats::sd(Y[up])
    zx_lo <- (X - mean(X[lo])) / stats::sd(X[lo])
    zy_lo <- (Y - mean(Y[lo])) / stats::sd(Y[lo])
    ind_up <- as.numeric(X >= stats::quantile(X, q)     & Y >= stats::quantile(Y, q))
    ind_lo <- as.numeric(X <= stats::quantile(X, 1 - q) & Y <= stats::quantile(Y, 1 - q))
    xi_j <- Tn / length(up) * (zx_up * zy_up - r_up) * ind_up -
            Tn / length(lo) * (zx_lo * zy_lo - r_lo) * ind_lo
    Xi <- cbind(Xi, xi_j)
    rhodiff <- c(rhodiff, r_up - r_lo)
    used <- c(used, q)
  }

  teststat <- NA_real_; pval <- NA_real_
  if (length(rhodiff) > 0) {
    om <- newey_west_cov(Xi)
    sol <- try(solve(om, rhodiff), silent = TRUE)
    if (!inherits(sol, "try-error")) {
      teststat <- Tn * sum(rhodiff * sol)
      pval <- 1 - stats::pchisq(teststat, df = length(rhodiff))
    }
  }
  list(table = res, teststat = teststat, pval = pval, thresholds_used = used)
}
