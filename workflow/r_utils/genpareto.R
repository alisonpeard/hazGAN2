"
Distribution definition script.

TODO: 

Must include functions:
- cdf(q:vector, params:list) -> vector
- threshold_selector(var:vector) -> list(params:list, p.value:float, pk:float)
"
library(eva)
library(goftest)
library(evd)

cdf <- function(q, params) {
  shape <- params$shape
  scale <- params$scale
  return(eva::pgpd(q, scale = scale, shape = shape, loc = 0))
}

threshold_selector <- function(var, nthresholds = 28, nsim = 5, alpha = 0.05) {
  gc(full = TRUE)
  thresholds <- quantile(var, probs = seq(0.7, 0.98, length.out = nthresholds))
  fits <- gpdSeqTestsWithFallback(var, thresholds, method = "ad", nsim = nsim)
  valid_pk <- fits$ForwardStop

  k    <- min(which(valid_pk > alpha)) # lowest index being "accepted"
  if (!is.finite(k)) {
    stop("All thresholds rejected under H0: X ~ GPD with α = 0.05")
    k <- 1
  }
  gc(full = TRUE)
  return(list(
    params   = list(
      thresh = fits$threshold[k],
      scale  = fits$est.scale[k],
      shape  = fits$est.shape[k]
    ),
    p.value  = fits$p.values[k],
    pk       = fits$ForwardStop[k]
  ))
}

gpdBackup <- function(var, threshold) {
  ad_test <- function(x, shape, scale, eps = 0.05) {
    cdf <- function(x) eva::pgpd(x, loc = 0, shape = shape, scale = scale)
    result <- goftest::ad.test(x, cdf)
    return(list(p.value = result$p.value))
  }

  # extract exceedances
  exceedances <- var[var > threshold] - threshold
  exceedances <- sort(exceedances)
  numabove <- length(exceedances)

  # TODO: alternative to POT::fitgpd --> evd::fpot(x, threshold)
  fit <- evd::fpot(var, threshold, std.err = FALSE)
  scale <- fit$param[1]
  shape <- fit$param[2]

  # goodness-of-fit test
  gof  <- ad_test(exceedances, shape, scale)

  # results
  return(list(
              thresh = threshold,
              shape = shape,
              scale = scale,
              p.value = gof$p.value,
              num.above = numabove))
}

gpdBackupSeqTests <- function(var, thresholds) {
  nthresh <- length(thresholds)
  shapes <- vector(length = nthresh)
  scales <- vector(length = nthresh)
  p.values <- vector(length = nthresh)
  num.above <- vector(length = nthresh)

  for (k in seq_along(thresholds)) {
    thresh        <- thresholds[k]
    fit           <- gpdBackup(var, thresh)
    shapes[k]     <- fit$shape
    scales[k]     <- fit$scale
    p.values[k]   <- fit$p.value
    num.above[k]  <- fit$num.above
  }

  ForwardStop <- rev(eva:::pSeqStop(rev(p.values))$ForwardStop)

  out <- list(threshold = thresholds, num.above = num.above,
              p.value = p.values,
              ForwardStop = ForwardStop, est.scale = scales,
              est.shape = shapes)
  return(as.data.frame(out))
}

gpdSeqTestsWithFallback <- function(var, thresholds, method, nsim) {
  gc(full = TRUE)
  fits <- tryCatch({
    fits <- eva::gpdSeqTests(var, thresholds = thresholds,
                             method = method, nsim = nsim)
    fits
  },
  error = function(e) {
    fits <- gpdBackupSeqTests(var, thresholds)
    fits
  })
  gc(full = TRUE)
  return(fits)
}