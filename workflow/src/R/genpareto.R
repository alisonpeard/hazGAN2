"
Distribution definition script.

TODO: 

Must include functions: 
- cdf(q:vector, params:list) -> vector
- threshold_selector(var:vector) -> list(params:list, p.value:float, pk:float)
"
suppressPackageStartupMessages({
  library(eva, quietly = TRUE)
  library(goftest, quietly = TRUE)
  library(evd, quietly = TRUE)
})

cdf <- function(q, params) {
  shape <- params$shape
  scale <- params$scale
  return(eva::pgpd(q, scale = scale, shape = shape, loc = 0))
}

threshold_selector <- function(var, nthresholds = 28, nsim = 1, alpha = 0.05) {
  
  thresholds <- quantile(var, probs = seq(0.7, 0.98, length.out = nthresholds))

  log_debug(paste0(
    "genpareto.R::threshold_selector - ",
    "Threshold selector called with: ",
    "\nVariable length: ", length(var),
    "\nVariable median: ", median(var),
    "\nVariable 90th percentile: ", quantile(var, 0.9),
    "\nVariable 95th percentile: ", quantile(var, 0.95),
    "\nVariable 98th percentile: ", quantile(var, 0.98),
    "\nVariable head: ", paste(head(var), collapse = ", "),
    "\nThresholds (tail): ", paste(tail(thresholds), collapse = ", ")

  ))

  fits <- gpdSeqTestsWithFallback(var, thresholds, method = "ad", nsim = nsim)
  valid_pk <- fits$ForwardStop

  k <- min(which(valid_pk > alpha)) # lowest index being "accepted"
  if (!is.finite(k)) {
    error_msg <- paste0(
      "genpareto.R::threshold_selector - ",
      "All thresholds rejected under H0: X ~ GPD with α = ", alpha,
      ".\nForward stop values: ", valid_pk, "\n\n"
    )
    stop(error_msg)
    k <- 1
  }
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

  if (numabove < 10) {
    return(NULL)
  }

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
  nthresh   <- length(thresholds)
  shapes    <- vector(length = nthresh)
  scales    <- vector(length = nthresh)
  p.values  <- vector(length = nthresh)
  num.above <- vector(length = nthresh)

  for (k in seq_along(thresholds)) {
    thresh        <- thresholds[k]
    fit           <- gpdBackup(var, thresh)
    if (!is.null(fit)) {
      shapes[k]     <- fit$shape
      scales[k]     <- fit$scale
      p.values[k]   <- fit$p.value
      num.above[k]  <- fit$num.above
    } else {
      shapes[k]     <- NA
      scales[k]     <- NA
      p.values[k]   <- NA
      num.above[k]  <- NA
    }
  }

  p.values <- p.values[!is.na(p.values)]

  # give up parametric fitting if all p-values are NA
  if (length(p.values) == 0) {
    msg <- paste0(
        "genpareto.R::gpdBackupSeqTests - ",
        "No valid p-values found.\n\n",
        "Variable statistics: ",
        "\n===========================",
        "\nThresholds checked: ", paste0(thresholds, collapse = ", "),
        "\nLength of thresholds: ", length(thresholds),
        "\nVariable length: ", length(var),
        '\nvar median: ', median(var),
        "\nvar 90th percentile: ", quantile(var, 0.9),
        "\nvar 95th percentile: ", quantile(var, 0.95),
        "\nvar 98th percentile: ", quantile(var, 0.98),
        "\nVariable head: ", paste0(head(var), collapse = ", "),
        "\n\n"
        )
    stop(msg)
  }

  ForwardStop <- rev(eva:::pSeqStop(rev(p.values))$ForwardStop)

  out <- list(threshold   = thresholds,
              num.above   = num.above,
              p.value     = p.values,
              ForwardStop = ForwardStop,
              est.scale   = scales,
              est.shape   = shapes)

  return(as.data.frame(out))
}

gpdSeqTestsWithFallback <- function(var, thresholds, method, nsim) {
  # try to run with eva first
  fits <- tryCatch({
    fits <- eva::gpdSeqTests(var, thresholds = thresholds,
                             method = method, nsim = nsim)
    fits
  }, # use evd::fpot as backup
  error = function(e) {
    # 2nd-lowest level error handling
    log_warn(paste0("genpareto.R::gpdSeqTestsWithFallback::eva::gpdSeqTests - ", e))
    fits <- gpdBackupSeqTests(var, thresholds)
    fits
  })
  return(fits)
}