suppressPackageStartupMessages({
  library(dplyr, quietly = TRUE)
  library(future, quietly = TRUE)
  library(furrr, quietly = TRUE)
  library(lubridate, quietly = TRUE)
  library(logger, quietly = TRUE)
})

`%ni%` <- Negate(`%in%`)

ecdf <- function(x) {
  # modified to Weibull plotting positions
  x <- sort(x)
  n <- length(x)
  if (n < 1) {
    stop("x must have ≥1 non-missing values")
  }
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals))) / (n + 1),
                    method = "constant",
                    rule = 2,
                    f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}


scdf <- function(train, params, cdf) {
  # Using excesses and setting loc=0
  # This is for flexibility with cdf choice
  params_upper <- params[
    names(params) %in% c("thresh_upper", "scale_upper", "shape_upper")
  ]
  params_lower <- params[
    names(params) %in% c("thresh_lower", "scale_lower", "shape_lower")
  ]

  # rename without "_upper" or "_lower"
  params_upper <- setNames(
    params_upper, c("thresh", "scale", "shape")
  )
  params_lower <- setNames(
    params_lower, c("thresh", "scale", "shape")
  )

  loc_upper <- params$thresh
  loc_lower <- params$thresh

  calculator <- function(x) {
    u <- ecdf(train)(x)

    mask_upper <- (x > loc_upper) & !is.na(loc_upper)
    mask_lower <- (x <= loc_lower) & !is.na(loc_lower)

    # fit the upper if there are valid parameters
    if (any(mask_upper)) {
      x_upper <- x[mask_upper]
      loc_upper <- loc_upper[mask_upper]
      params_upper <- params_upper[mask_upper]

      exceedances_upper <- x_upper - loc_upper
      pthresh_upper <- ecdf(train)(loc_upper)
      u_upper <- 1 - (1 - pthresh_upper)

      u_upper <- u_upper * (1 - cdf(exceedances_upper, params_upper))
      u[mask_upper] <- u_upper
    }

    # fit the lower if there are valid parameters
    if (any(mask_lower)) {
      x_lower <- x[mask_lower]
      loc_lower <- loc_lower[mask_lower]
      params_lower <- params_lower[mask_lower]

      exceedances_lower <- loc_lower - x_lower
      pthresh_lower <- ecdf(train)(loc_lower)
      u_lower <- pthresh_lower * (1 - cdf(exceedances_lower, params_lower))
      u[mask_lower] <- u_lower
    }
    return(u)
  }
  return(calculator)
}


ljung_box <- function(variable, threshold, grid_i, tail) {
  # test independence of exceedences
  excesses <- variable[variable > threshold]
  nexcesses <- length(excesses)

  if (nexcesses < 30) {
    log_warn(paste0(
      "stats::ljung_box - ",
      "Only ", nexcesses,
      " exceedances in gridcell ",
      grid_i, " (", tail, ")", ". ",
      "Skipping Ljung-Box test."
    ))
    p_box <- 0
  } else {
    # log_debug(paste(excesses, collapse = ", "))
    p_box <- Box.test(excesses)[["p.value"]]

    if (is.na(p_box)) {
      log_warn(paste0(
        "stats::fit_gridcell - ",
        "Ljung-Box test returned NA for gridcell ",
        grid_i, " (", tail, ")", ". Value: ", round(p_box, 4),
        ", Num exceedances: ", nexcesses,
        ", Mean exceedance: ", round(mean(excesses), 4)
      ))
    } else if (p_box < 0.1) {
      log_warn(paste0(
        "stats::fit_gridcell - ",
        "p-value ≤ 10% for H0: independent exceedences in ",
        grid_i, " (", tail, ")", ". Value: ", round(p_box, 4),
        ", Num exceedances: ", nexcesses,
        ", Mean exceedance: ", round(mean(excesses), 4)
      ))
    } else {
      log_success(paste0(
        "stats::fit_gridcell - ",
        "p-value > 10% for H0: independent exceedences in ",
        grid_i, " (", tail, ")", ". Value: ", round(p_box, 4),
        ", Num exceedances: ", nexcesses,
        ", Mean exceedance: ", round(mean(excesses), 4)
      ))
    }
  }
  return(p_box)
}