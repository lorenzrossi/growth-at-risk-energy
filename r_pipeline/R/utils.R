# utils.R -- shared utilities for the GaR project (R port)
# Replaces: netIncrease.m, blockBootIdx.m, printpdf.m

#' Hamilton (2003) net price increase on the 100*log scale.
#' Single canonical construction for both NGPI and NOPI.
#'
#' @param P numeric vector of PRICE LEVELS (not logs, not differences)
#' @param n_months trailing window (default 12)
#' @return numeric vector, zero for the first n_months observations
net_increase <- function(P, n_months = 12L) {
  if (any(P <= 0, na.rm = TRUE)) stop("P must contain strictly positive price levels.")
  logP <- 100 * log(P)
  n <- length(logP)
  NI <- numeric(n)
  for (t in (n_months + 1L):n) {
    NI[t] <- max(0, logP[t] - max(logP[(t - n_months):(t - 1L)]))
  }
  NI
}

#' Moving-block bootstrap indices (Kunsch, 1989).
#' Use instead of iid resampling whenever residuals are serially dependent
#' (direct multi-horizon projections in particular).
#'
#' @param n sample size
#' @param block_length block size; recommended max(12, h) for monthly data
block_boot_idx <- function(n, block_length = 12L) {
  block_length <- min(max(1L, as.integer(round(block_length))), n)
  n_blocks <- ceiling(n / block_length)
  starts <- sample.int(n - block_length + 1L, n_blocks, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) s:(s + block_length - 1L)))
  idx[seq_len(n)]
}

#' 100*log first difference (monthly log growth in percent)
log_diff <- function(x) c(NA_real_, 100 * diff(log(x)))

#' Trailing h-month mean, aligned like MATLAB filter(ones(1,h)/h, 1, y):
#' out[t] = mean(y[(t-h+1):t]), NA for t < h.
trailing_mean <- function(y, h) {
  n <- length(y)
  out <- rep(NA_real_, n)
  if (h == 1) return(y)
  cs <- cumsum(y)
  out[h:n] <- (cs[h:n] - c(0, cs[seq_len(n - h)])) / h
  out
}

#' Save a ggplot as PDF (replaces printpdf.m)
save_pdf <- function(plot, filename, width = 8, height = 5.5) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, device = "pdf")
  invisible(filename)
}

#' Restrict a data.frame with a Time column to Jan 2004 -- Dec 2023
subsample_window <- function(df, from = "2004-01-01", to = "2023-12-31") {
  df[df$Time >= as.Date(from) & df$Time <= as.Date(to), , drop = FALSE]
}
