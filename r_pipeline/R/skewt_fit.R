# skewt_fit.R -- fit Azzalini skew-t distributions to estimated conditional
# quantiles (step 2 of Adrian-Boyarchenko-Giannone) and compute GaR,
# Expected Shortfall / Longrise, entropy, and moments.
# Replaces: Step2match_improved.m, QuantilesInterpolation_improved.m and the
# entire hand-ported azzalini/ suite (dskt/pskt/qskt/rskt) via the sn package
# (Azzalini's own reference implementation: sn::dst, sn::pst, sn::qst).

#' Fit a skew-t to a set of target quantiles by nonlinear least squares,
#' profiling over integer degrees of freedom (as in the AER replication).
#'
#' @param q_targ  vector of target quantiles (same length as taus)
#' @param taus    quantile levels corresponding to q_targ
#' @param q_match subset of levels used for matching
#' @param nu_grid integer grid for degrees of freedom
#' @return list(xi, omega, alpha, nu, ssq, crossed)
fit_skewt_quantiles <- function(q_targ, taus,
                                q_match = c(.05, .10, .25, .50, .75, .90, .95),
                                nu_grid = 2:20) {
  stopifnot(requireNamespace("sn", quietly = TRUE))

  crossed <- any(diff(q_targ) < 0)
  if (crossed) q_targ <- sort(q_targ)   # rearrangement (CFG 2010)

  sel <- vapply(q_match, function(q) which.min(abs(taus - q)), integer(1))
  p_sel <- taus[sel]; q_sel <- q_targ[sel]

  # initial conditions: median and normal-scaled IQR
  xi0 <- stats::median(q_targ)
  om0 <- max((stats::quantile(q_targ, .75) - stats::quantile(q_targ, .25)) /
               (stats::qnorm(.75) - stats::qnorm(.25)), 0.1)

  obj <- function(par, nu) {
    qs <- try(sn::qst(p_sel, xi = par[1], omega = par[2], alpha = par[3], nu = nu),
              silent = TRUE)
    if (inherits(qs, "try-error") || any(!is.finite(qs))) return(1e12)
    sum((q_sel - qs)^2)
  }

  best <- list(ssq = Inf)
  for (nu in nu_grid) {
    fit <- try(stats::optim(c(xi0, om0, 0), obj, nu = nu, method = "L-BFGS-B",
                            lower = c(-100, 1e-4, -100), upper = c(100, 200, 100)),
               silent = TRUE)
    if (inherits(fit, "try-error")) next
    if (fit$value < best$ssq) {
      best <- list(xi = fit$par[1], omega = fit$par[2], alpha = fit$par[3],
                   nu = nu, ssq = fit$value)
    }
  }
  best$crossed <- crossed
  best
}

#' Step 2: fit skew-t distributions to each row of a matrix of conditional
#' quantiles; compute densities, moments, entropy, GaR and ES/Longrise.
#'
#' @param YQ      T x ntau matrix of fitted conditional quantiles (rows with NA skipped)
#' @param yq_unc  vector of unconditional quantiles (same taus)
#' @param taus    quantile levels
#' @param yy      evaluation grid for the density
#' @param alpha_es tail probability for ES / Longrise / GaR (default 0.05)
#' @return list with matrices/vectors mirroring Step2match's Res struct
step2_match <- function(YQ, yq_unc, taus,
                        yy = seq(-30, 30, by = 0.02), alpha_es = 0.05,
                        q_match = c(.05, .10, .25, .50, .75, .90, .95),
                        nu_grid = 2:20, verbose = TRUE) {
  stopifnot(requireNamespace("sn", quietly = TRUE))
  Tn <- nrow(YQ); ntau <- length(taus)
  delta <- yy[2] - yy[1]
  jq50 <- which.min(abs(taus - 0.5))

  # unconditional fit
  un <- fit_skewt_quantiles(yq_unc, taus, q_match, nu_grid)
  p_unc <- sn::dst(yy, xi = un$xi, omega = un$omega, alpha = un$alpha, nu = un$nu)

  par_mat <- matrix(NA_real_, Tn, 4,
                    dimnames = list(NULL, c("location", "scale", "shape", "df")))
  PST <- matrix(NA_real_, Tn, length(yy))
  QST <- matrix(NA_real_, Tn, ntau)
  mom <- matrix(NA_real_, Tn, 4,
                dimnames = list(NULL, c("mean", "variance", "skewness", "kurtosis")))
  ent <- matrix(NA_real_, Tn, 2, dimnames = list(NULL, c("left", "right")))
  es  <- rep(NA_real_, Tn)   # expected shortfall
  el  <- rep(NA_real_, Tn)   # expected longrise
  gar <- rep(NA_real_, Tn)   # Growth-at-Risk (alpha_es quantile of fitted dist.)
  n_crossed <- 0L

  p_grid <- seq(0.01, alpha_es, by = 0.01)   # for ES integral (delta1 = 0.01)

  for (t in seq_len(Tn)) {
    q_targ <- YQ[t, ]
    if (any(is.na(q_targ))) next
    ft <- fit_skewt_quantiles(q_targ, taus, q_match, nu_grid)
    if (!is.finite(ft$ssq)) next
    if (isTRUE(ft$crossed)) n_crossed <- n_crossed + 1L
    par_mat[t, ] <- c(ft$xi, ft$omega, ft$alpha, ft$nu)

    d <- sn::dst(yy, xi = ft$xi, omega = ft$omega, alpha = ft$alpha, nu = ft$nu)
    PST[t, ] <- d
    QST[t, ] <- sn::qst(taus, xi = ft$xi, omega = ft$omega, alpha = ft$alpha, nu = ft$nu)

    m1 <- sum(yy * d * delta)
    v  <- sum((yy - m1)^2 * d * delta)
    mom[t, ] <- c(m1, v,
                  sum((yy - m1)^3 * d * delta) / v^1.5,
                  sum((yy - m1)^4 * d * delta) / v^2)

    med <- QST[t, jq50]
    lg <- log(d) - log(p_unc)
    lg[!is.finite(lg)] <- 0
    ent[t, "left"]  <- sum(lg * d * (yy < med) * delta)
    ent[t, "right"] <- sum(lg * d * (yy > med) * delta)

    qq_lo <- sn::qst(p_grid,      xi = ft$xi, omega = ft$omega, alpha = ft$alpha, nu = ft$nu)
    qq_hi <- sn::qst(1 - p_grid,  xi = ft$xi, omega = ft$omega, alpha = ft$alpha, nu = ft$nu)
    es[t]  <- mean(qq_lo)
    el[t]  <- mean(qq_hi)
    gar[t] <- sn::qst(alpha_es, xi = ft$xi, omega = ft$omega, alpha = ft$alpha, nu = ft$nu)

    if (verbose && t %% 24 == 0)
      message(sprintf("  skew-t fit: %d / %d", t, Tn))
  }

  if (n_crossed > 0)
    warning(sprintf("Quantile crossing rearranged in %d of %d observations.", n_crossed, Tn))

  list(par = par_mat, PST = PST, QST = QST, yy = yy, moments = mom,
       entropy = ent, shortfall = es, longrise = el, gar = gar,
       unc = un, p_unc = p_unc, n_crossed = n_crossed, taus = taus)
}
