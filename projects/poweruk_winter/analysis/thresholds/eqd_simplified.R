# Negative log-likelihood functions
nll <- function(par, x, distn="genpareto") {
  if (distn == "genpareto") {
    nll_genpareto(par, x)
  } else if (distn == "weibull") {
    nll_weibull(par, x)
  }
}

nll_genpareto <- function(params, x) {
  scale <- params[1]
  shape <- params[2]
  
  # Check scale is positive
  if (scale <= 0) return(1e10)
  
  # Check GPD support: when shape < 0, need x < -scale/shape
  if (shape < 0 && any(x >= -scale/shape)) {
    return(1e10)
  }
  
  # Calculate log-likelihood with error handling
  ll <- tryCatch({
    sum(eva::dgpd(x, scale = scale, shape = shape, log = TRUE))
  }, error = function(e) {
    return(-Inf)
  }, warning = function(w) {
    return(-Inf)
  })
  
  # Check for invalid values
  if (!is.finite(ll) || is.na(ll)) return(1e10)
  
  # Return NEGATIVE log-likelihood (for minimization)
  return(-ll)
}

nll_weibull <- function(params, x) {
  scale <- params[1]
  shape <- params[2]
  
  if (scale <= 0 || shape <= 0) return(1e10)
  
  ll <- tryCatch({
    sum(dweibull(x, shape = shape, scale = scale, log = TRUE))
  }, error = function(e) {
    return(-Inf)
  }, warning = function(w) {
    return(-Inf)
  })
  
  if (!is.finite(ll) || is.na(ll)) return(1e10)
  
  return(-ll)
}

# Quantile function
qdistn <- function(x, params, distn="genpareto") {
  if (distn == "genpareto") {
    scale <- params[1]
    shape <- params[2]
    eva::qgpd(x, scale = scale, shape = shape)
  } else if (distn == "weibull") {
    scale <- params[1]
    shape <- params[2]
    qweibull(x, shape = shape, scale = scale)  # FIXED
  }
}

# Main function
eqd <- function(data, thresh, B = 100, m = 500, distn = "genpareto") {
  
  # Check inputs are valid
  if (!is.numeric(data)) stop("Data must be a vector")
  if (!is.numeric(thresh)) stop("u to be tested needs to be a vector")
  if (B <= 0 | B %% 1 != 0) stop("Number of bootstrapped samples must be a positive integer")
  if (m <= 0 | m %% 1 != 0) stop("Number of equally spaced probabilities must be a positive integer")
  
  meandistances <- xis <- sigmas <- num_excess <- numeric(length(thresh))
  
  if (distn == "genpareto") {
    lower <- c(1e-6, -Inf)
  } else if (distn == "weibull") {
    lower <- c(1e-6, 1e-6)
  }
  
  upper <- c(Inf, Inf)
  
  for (i in 1:length(thresh)) {
    u <- thresh[i]
    excess <- data[data > u] - u
    num_excess[i] <- length(excess)
    
    if (num_excess[i] > 10) {
      
      scale0 <- mean(excess)
      if (distn == "genpareto") {
        par_init <- c(scale0, 0.1) 
      } else if (distn == "weibull") {
        par_init <- c(scale0, 1)
      }
      
      # Initial fit
      init.fit <- optim(par = par_init, fn = nll, x = excess, distn = distn,
                        method = "L-BFGS-B", lower = lower, upper = upper)
      
      sigmas[i] <- init.fit$par[1]
      xis[i] <- init.fit$par[2]
      distances <- numeric(B)
      
      # Bootstrap loop
      for (b in 1:B) {
        X <- sample(excess, num_excess[i], replace = TRUE)
        scale0 <- mean(X)
        
        # FIXED: Proper if-else instead of ifelse
        if (distn == "genpareto") {
          if (xis[i] < 0) {
            pars_init <- c(scale0, 0.1)
          } else {
            pars_init <- c(sigmas[i], xis[i])
          }
        } else if (distn == "weibull") {
          pars_init <- c(sigmas[i], xis[i])
        }
        
        # Bootstrap fit
        fit <- optim(par = pars_init, fn = nll, x = X, distn = distn,
                     method = "L-BFGS-B", lower = lower, upper = upper)
        
        quants <- qdistn((1:m)/(m + 1), params = fit$par, distn = distn)
        distances[b] <- (1/m) * sum(abs(quantile(X, probs = (1:m)/(m+1)) - quants))
      }
      meandistances[i] <- mean(distances)
    } else {
      meandistances[i] <- NA
    }
  }
  
  chosen_index <- which.min(meandistances)
  chosen_threshold <- thresh[chosen_index]
  xi <- xis[chosen_index]
  sigma <- sigmas[chosen_index]
  len <- num_excess[chosen_index]
  
  result <- list(thresh = chosen_threshold, par = c(sigma, xi), 
                 num_excess = len, dists = meandistances)
  return(result)
}

# Test
set.seed(123)
#x <- eva::rgpd(n = 500, scale = 0.5, shape = 0.1)
x <- rweibull(n = 500, scale = 0.5, shape = 2)
thresholds <- quantile(x, probs = seq(0.0, 0.98, length.out = 20))
fits <- eqd(x, thresholds, B = 50, distn = "weibull")  # Reduced B for testing
fits
