suppressPackageStartupMessages({
  library(eva, quietly = TRUE)
  library(goftest, quietly = TRUE)
  library(evd, quietly = TRUE)
  library(logger, quietly = TRUE)
})

cdf <- function(q, params) {
  shape <- params$shape
  scale <- params$scale
  return(eva::pgpd(q, scale = scale, shape = shape, loc = 0))
}

ad_test <- function(x, shape, scale, eps = 0.05) {
  cdf <- function(x) eva::pgpd(x, loc = 0, shape = shape, scale = scale)
  result <- goftest::ad.test(x, cdf)
  return(list(p.value = result$p.value))
}

threshold_selector <- function(
  var, id, nthresholds = 28, nsim = 1, alpha = 0.05
) {
  thresholds <- quantile(var, probs = seq(0.7, 0.98, length.out = nthresholds))
  fits <- fits_with_fallback(var, thresholds, nsim = nsim)

  valid_pk <- fits$ForwardStop
  
  k <- min(which(valid_pk > alpha)) # lowest valid threshold
  if (!is.finite(k))  {
    error_msg <- paste0(
      "All thresholds for ", id, " rejected under H0: X ~ GPD."
    )
    stop(error_msg)
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


fits_with_fallback <- function(var, thresholds, nsim) {
  fits <- tryCatch({
    eva::gpdSeqTests(
      var, thresholds = thresholds, method = "ad", nsim = nsim
    )
  },
  error = function(e) {
    log_warn(paste0(
      "\neva::gpdSeqTests fit failed.\n...", e,
      "...Falling back to eva::fpot"
    ))
    fallback_fits(var, thresholds)
  })
  return(fits)
}


fallback_fits <- function(var, thresholds) {
  nthresh   <- length(thresholds)
  shapes    <- vector(length = nthresh)
  scales    <- vector(length = nthresh)
  p_values  <- vector(length = nthresh)
  numabove  <- vector(length = nthresh)

  for (k in seq_along(thresholds)) {

    thresh <- thresholds[k]
    fit <- fallback_genpareto(var, thresh)

    if (!is.null(fit)) {
      shapes[k]   <- fit$shape
      scales[k]   <- fit$scale
      p_values[k] <- fit$p.value
      numabove[k] <- fit$num.above
    } else {
      shapes[k]   <- NA
      scales[k]   <- NA
      p_values[k] <- NA
      numabove[k] <- NA
    }
  }

  if (all(is.na(p_values))) {
    stop("All fallback fits failed.")
  }

  valid_mask <- which(!is.na(p_values))
  pk_tmp <- rev(eva:::pSeqStop(rev(p_values[valid_mask]))$ForwardStop)
  pk_values <- rep(NA, nthresh)
  pk_values[valid_mask] <- pk_tmp

  out <- list(
    threshold   = thresholds,
    num.above   = numabove,
    p.value     = p_values,
    ForwardStop = pk_values,
    est.scale   = scales,
    est.shape   = shapes
  )

  return(as.data.frame(out))
}


fallback_genpareto <- function(var, threshold) {
  exceedances <- var[var > threshold] - threshold
  exceedances <- sort(exceedances)
  numabove    <- length(exceedances)

  if (numabove < 10) {
    log_warn(paste0(
      "Only ", numabove, " exceedances",
      " for threshold ", threshold,
      ". Skipping..."
    ))
    return(NULL)
  }

  fit <- tryCatch({
    evd::fpot(var, threshold, std.err = FALSE)
  },
  error = function(e) {
    log_warn(paste0("evd::fpot fit failed: ", e))
    return(NULL)
  })

  scale <- fit$param[1]
  shape <- fit$param[2]

  gof  <- ad_test(exceedances, shape, scale)

  return(list(
    thresh = threshold,
    shape = shape,
    scale = scale,
    p.value = gof$p.value,
    num.above = numabove
  ))
}
