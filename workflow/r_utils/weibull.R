"
Distribution definition script.

Must include functions:
- cdf(q:vector, params:list) -> vector
- threshold_selector(var:vector) -> list(params:list, p.value:float, pk:float)
"
library(goftest)

cdf <- function(q, params) {
  shape <- params$shape
  scale <- params$scale
  return(pweibull(q, shape = shape, scale = scale))
}

threshold_selector <- function(var, alpha = 0.05, nthresholds = 50) {
  loglikelihood <- function(params, data) {
    # Calculate negative Weibull log-likelihood.
    shape <- params[1]
    scale <- params[2]
    -sum(dweibull(data, shape = shape, scale = scale, log = TRUE))
  }
  ad_test <- function(x, shape, scale, eps=0.05){
    cdf <- function(x) pweibull(x, shape = shape, scale = scale)
    result <- goftest::ad.test(x, cdf)
    return(list(p.value = result$p.value))
  }

  thresholds <- quantile(var, probs = seq(0.7, 0.98, length.out = nthresholds))

  shapes    <- vector(length = length(thresholds))
  scales    <- vector(length = length(thresholds))
  num_above <- vector(length = length(thresholds))
  pvals     <- vector(length = length(thresholds))

  for (i in seq_along(thresholds)) {
    q <- thresholds[i]
    exceedances <- var[var > q] - q

    # initial estimates
    mean_exc <- mean(exceedances)
    init_shape <- 2  # like Rayleigh distribution for winds
    init_scale <- mean_exc

    # fit MLE
    fit <- optim(c(init_shape, init_scale),
                 loglikelihood,
                 data = exceedances,
                 method = "L-BFGS-B",
                 lower = c(0.1, 0.1),
                 upper = c(2, 10))

    shapes[i]      <- fit$par[1]
    scales[i]      <- fit$par[2]
    num_above[i]   <- length(exceedances)
    pvals[i]       <- ad_test(exceedances, shapes[i], scales[i])$p.value
  }

  # https://doi.org/10.1214/17-AOAS1092
  pk <- rev(eva:::pSeqStop(rev(pvals))$ForwardStop)
  # m   <- length(pvals)
  # int <- seq(1, m, 1)
  # pk  <- cumsum(-log(1 - pvals[int]))/int

  k   <- min(which(pk > alpha)) # lowest index being "accepted"
  if (!is.finite(k)) {
    stop("All thresholds rejected under H0:X~Weibull with α=0.05")
    k <- 1
  }

  thresh <- thresholds[k]
  shape  <- shapes[k]
  scale  <- scales[k]
  pval   <- pvals[k]
  pk     <- pk[k]
  exceedances <- var[var > thresh] - thresh
  num_above   <- length(exceedances)

  return(list(
    params = list(
      thresh = thresh,
      scale = scale,
      shape = shape
    ),
    p.value = pval,
    pk = pk,
    n_exceed = num_above
  ))
}
