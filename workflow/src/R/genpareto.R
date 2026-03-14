# Threshold selection EQD method
# GPD nll + penalty terms => strict CDF(x) < 1
# doi:10.1080/00401706.2024.2421744

suppressPackageStartupMessages({
  library(eva, quietly = TRUE)
  library(goftest, quietly = TRUE)
  library(evd, quietly = TRUE)
  library(logger, quietly = TRUE)
})


cdf <- function(q, params) {
  shape <- params$shape
  scale <- params$scale
  eva::pgpd(q, scale = scale, shape = shape, loc = 0)
}


nll_genpareto <- function(params, x, max_x) {
  # genpareto with strict CDF(x) < 1
  scale <- params[1]
  shape <- params[2]

  # penalty terms
  if (scale <= 0) return(1e10)
  if (shape < 0) {
    if (max_x >= (-scale / shape)) return(1e10)
  }
  
  ll <- sum(eva::dgpd(x, scale = scale, shape = shape, log = TRUE))

  # stability check
  if (is.na(ll) || is.infinite(ll)) return(1e10)
  -ll
}


fit_genpareto <- function(x) {
  init_scale <- mean(x)
  init_shape <- 0.1

  fit <- tryCatch({
    optim(
      c(init_scale, init_shape),
      nll_genpareto,
      x = x,
      method = "L-BFGS-B",
      lower = c(1e-6, -0.9),
      upper = c(Inf, 0.9),
      max_x = max(x)
    )
  },
  error = function(e) {
    log_warn(paste0("optim fit failed: ", e))
    return(NULL)
  })

  if (is.null(fit)) {
    return(NULL)
  }

  scale <- fit$par[1]
  shape <- fit$par[2]

  list(param = c(scale, shape))
}


eqd_genpareto <- function(x, thresh, nboot = 100, m = 100) {
  # doi: 10.1080/00401706.2024.2421744
  meandists <- rep(Inf, length(thresh))
  sddists   <- rep(NA, length(thresh))
  scales    <- rep(NA, length(thresh))
  shapes    <- rep(NA, length(thresh))
  numaboves <- rep(NA, length(thresh))
  
  p_grid    <- (1:m) / (m + 1)
  
  for (i in seq_along(thresh)) {
    u <- thresh[i]
    excess <- x[x > u] - u
    numaboves[i] <- length(excess)
    
    if (numaboves[i] > 20) {
      fit0 <- fit_genpareto(excess)
      if(is.null(fit0)) next
      scales[i] <- fit0$param[1]
      shapes[i] <- fit0$param[2]
      
      par_init <- if (shapes[i] < 0) c(mean(excess), 0.1) else c(scales[i], shapes[i])
      
      # bootstrap
      dists <- vapply(seq_len(nboot), function(b) {
        xb <- sample(excess, numaboves[i], replace = TRUE)
        max_xb <- max(xb)
      
        fit_b <- optim(
          par_init,
          nll_genpareto,
          x = xb,
          method = "L-BFGS-B",
          lower = c(1e-6, -0.9),
          upper = c(Inf, 0.9),
          max_x = max_xb
        )
        
        if (fit_b$convergence != 0) return(NA)
        
        q_emp <- quantile(xb, probs = p_grid, names = FALSE)
        q_gpd <- eva::qgpd(p_grid, scale = fit_b$par[1], shape = fit_b$par[2])
        mean(abs(q_emp - q_gpd))
      }, FUN.VALUE = numeric(1))
      
      meandists[i] <- mean(dists, na.rm = TRUE)
      sddists[i]   <- sd(dists, na.rm = TRUE)
    }
  }
  i_best <- which.min(meandists)
  if (length(i_best) == 0) return(NULL)
  if (is.infinite(meandists[i_best])) return(NULL)
  
  thresh <- thresh[i_best]
  scale  <- scales[i_best]
  shape  <- shapes[i_best]
  numabove <- numaboves[i_best]

  list(
    thresh     = unname(thresh),
    scale      = scale,
    shape      = shape,
    numabove   = numabove
    # diagnostic = list(
    #   meandist = meandists[i_best],
    #   sddist   = sddists[i_best]
    # )
  )
}


ad_test <- function(x, scale, shape, eps = 0.05) {
  cdf <- function(x) eva::pgpd(x, loc = 0, scale = scale, shape = shape)
  result <- goftest::ad.test(x, cdf)
  result$p.value
}


threshold_selector <- function(
  x, id, nthresholds = 15, nsim = 30, alpha = 0.05
) {
  thresholds <- quantile(
    x, probs = seq(0.7, 0.98, length.out = nthresholds)
  )
  fit <- eqd_genpareto(x, thresholds, nboot = nsim)

  if (is.null(fit) || is.null(fit$thresh)) {
    return(list(
      params = list(thresh = NA, scale = NA, shape = NA),
      p.value = NA,
      status = "Failure: No valid threshold found with enough data"
    ))
  }
  
  p <- tryCatch(
    ad_test(
      x = x[x > fit$thresh] - fit$thresh,
      scale = fit$scale,
      shape = fit$shape
    ),
    error = function(e) NA
  )

  return(list(
    params   = list(
      thresh = fit$thresh,
      scale  = fit$scale,
      shape  = fit$shape
      # diagnostic = fit$diagnostic
    ),
    p.value = p,
    pk = NA,
    status = "Success"
  ))
}


if (FALSE) {
  # run diagnostic tests
  library(eva)
  print("My case")
  x <- rnorm(1249)
  threshold_selector(x, "1249 points");
  
  print("too few points case:")
  x <- rnorm(10)
  threshold_selector(x, "too-few points");
  
  print("zero-variance tail case:")
  x <- c(rnorm(100), rep(5, 50))
  threshold_selector(x, "zero-variance tail");
  
  print("extremely heavy (Cauchy) tail (ξ<-0.5)")
  x <- abs(rcauchy(200))
  threshold_selector(x, "heavy tail");
}