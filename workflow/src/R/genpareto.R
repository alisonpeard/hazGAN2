# Threshold selection EQD method
# GPD nll + penalty terms => strict CDF(x) < 1
# NB: Hardcoded to avoid negative shape params
# EQD method: 10.1080/00401706.2024.2421744

suppressPackageStartupMessages({
  library(eva, quietly = TRUE)
  library(goftest, quietly = TRUE)
  library(evd, quietly = TRUE)
  library(logger, quietly = TRUE)
})


cdf <- function(q, params) {
  shape <- params$shape
  scale <- params$scale
  loc   <- params$loc
  eva::pgpd(q, scale = scale, shape = shape, loc = loc)
}


nll_genpareto <- function(params, x, xmax) {
  # genpareto with strict CDF(x) < 1
  scale <- params[1]
  shape <- params[2]
  loc   <- params[3]

  if (scale <= 0) return(1e10)
  if (any(x < loc)) return(1e10)

  ll <- sum(eva::dgpd(x, scale = scale, shape = shape, loc = loc, log = TRUE))
  if (is.na(ll) || is.infinite(ll)) return(1e10)

  # penalty blows up as upper bound approaches xmax
  penalty <- 0
  if (shape < 0) {
    upper <- loc - scale / shape
    margin <- (upper - xmax) / upper
    if (margin <= 0) return(1e10)
    penalty <- 1 / margin

  -ll + penalty
}
}


fit_genpareto <- function(
  x, lower = c(1e-6, -0.9, -1e-6), upper = c(Inf, 0.9, 1e-6)
) {
  init_scale <- mean(x)
  init_shape <- 0.1
  init_loc   <- 0

  fit <- tryCatch({
    optim(
      c(init_scale, init_shape, init_loc),
      nll_genpareto,
      x = x,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      xmax = max(x)
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
  loc   <- fit$par[3]

  list(param = c(scale, shape, loc))
}


eqd_genpareto <- function(x, thresh, nboot = 100, m = 100) {
  # doi: 10.1080/00401706.2024.2421744
  meandists <- rep(Inf, length(thresh))
  sddists   <- rep(NA, length(thresh))
  scales    <- rep(NA, length(thresh))
  shapes    <- rep(NA, length(thresh))
  locs      <- rep(NA, length(thresh))
  numaboves <- rep(NA, length(thresh))

  p_grid    <- (1:m) / (m + 1)

  for (i in seq_along(thresh)) {
    u <- thresh[i]
    excess <- x[x > u] - u
    numaboves[i] <- length(excess)
    xmax <- max(excess)

    if (numaboves[i] > 20) {
      fit0 <- fit_genpareto(excess, lower = c(1e-6, -0.9, -1e-6))
      if (is.null(fit0)) next
      scales[i] <- fit0$param[1]
      shapes[i] <- fit0$param[2]
      locs[i]   <- fit0$param[3]

      par_init <- if (shapes[i] < 0) {
        c(mean(excess), 0.1, 0)
      } else {
        c(scales[i], shapes[i], 0)
      }

      # bootstrap
      dists <- vapply(seq_len(nboot), function(b) {
        xb <- sample(excess, numaboves[i], replace = TRUE)

        fit_b <- optim(
          par_init,
          nll_genpareto,
          x = xb,
          method = "L-BFGS-B",
          lower = c(1e-6, -0.9, -1e-6),
          upper = c(Inf, 0.9, 1e-6),
          xmax = xmax
        )

        if (fit_b$convergence != 0) return(NA)

        q_emp <- quantile(xb, probs = p_grid, names = FALSE)
        q_gpd <- eva::qgpd(p_grid, scale = fit_b$par[1], shape = fit_b$par[2], loc = fit_b$par[3])
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
  loc    <- locs[i_best]
  numabove <- numaboves[i_best]

  list(
    thresh     = unname(thresh),
    scale      = scale,
    shape      = shape,
    loc        = loc,
    numabove   = numabove
  )
}


ad_test <- function(x, scale, shape, loc = 0, eps = 0.05) {
  cdf <- function(x) eva::pgpd(x, loc = loc, scale = scale, shape = shape)
  result <- goftest::ad.test(x, cdf)
  result$p.value
}


# full params:
# x, id, nthresholds = 15, nsim = 10, alpha = 0.05
# quick run params:
# x, id, nthresholds = 8, nsim = 1, alpha = 0.05
threshold_selector <- function(
  x, id, nthresholds = 8, nsim = 1, alpha = 0.05
) {
  thresholds <- quantile(
    x, probs = seq(0.7, 0.98, length.out = nthresholds)
  )
  fit <- eqd_genpareto(x, thresholds, nboot = nsim)

  if (is.null(fit) || is.null(fit$thresh)) {
    return(list(
      params = list(thresh = NA, scale = NA, shape = NA, loc = NA),
      p.value = NA,
      status = "Failure: No valid threshold found with enough data"
    ))
  }

  # manual AD test now that not using
  # Bader method anymore
  p <- tryCatch(
    ad_test(
      x = x[x > fit$thresh] - fit$thresh,
      scale = fit$scale,
      shape = fit$shape,
      loc   = fit$loc
    ),
    error = function(e) NA
  )

  return(list(
    params   = list(
      thresh = fit$thresh,
      scale  = fit$scale,
      shape  = fit$shape,
      loc    = fit$loc
      # diagnostic = fit$diagnostic
    ),
    p.value = p,
    pk = p,
    status = "Success"
  ))
}