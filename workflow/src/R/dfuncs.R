"Functions related to creating and declustering an event intensity
time series. Coles (2001) runs declustering."
suppressPackageStartupMessages({
  library(extRemes)
})
source("workflow/src/R/rfuncs.R")

`%ni%` <- Negate(`%in%`)

gridsearch <- function(series, var, qmin = 60, qmax = 99, rmin = 1, rmax = 14) {
  "Unit tests for this?"
  qvec <- c(qmin:qmax) / 100
  rvec <- c(rmin:rmax)

  nclusters <- matrix(nrow = length(rvec), ncol = length(qvec))
  ext_ind   <- matrix(nrow = length(rvec), ncol = length(qvec))
  pvals     <- matrix(nrow = length(rvec), ncol = length(qvec))

  series_var <- series[[var]]
  thresholds  <- quantile(series_var, qvec)

  print("Testing combinations:")
  for (i in seq_along(rvec)){
    for (j in seq_along(qvec)){
      thresh <- thresholds[j]

      d <- decluster(series_var, thresh = thresh,
                     r = rvec[i], method = "runs")

      # NOTE: theta = 1 a lot, double-check?
      e <- extremalindex(c(d), thresh, r = rvec[i],
                         method = "runs") # Coles (2001) §5.3.2

      p <- Box.test(c(d)[c(d) > thresh], type = "Ljung")

      nclusters[i, j] <- e[["number.of.clusters"]]
      ext_ind[i, j]   <- e[["extremal.index"]]
      pvals[i, j]     <- p$p.value
    }
  }
  print("Before filtering:")
  print(table(is.finite(nclusters)))

  print("After extremal index filter:")
  nclusters[ext_ind < 0.8] <- -Inf # theta < 1 => extremal dependence
  print(table(is.finite(nclusters)))

  print("After p-value filter:")
  nclusters[pvals < 0.1]  <- -Inf  # H0: independent exceedances
  print(table(is.finite(nclusters)))

  ind <- which(nclusters == max(nclusters), arr.ind = TRUE)
  r <- rvec[ind[1]]
  q <- qvec[ind[2]]
  p <- pvals[ind[1], ind[2]]
  return(list(r = r, q = q, p = p))
}

identify_events <- function(daily, rfunc) {
  args    <- rfunc$args
  rfunc   <- match.fun(rfunc$func)
  series  <- aggregate(. ~ time, daily[, c("time", args)], rfunc)

  # gridsearch run lengths and thresholds
  result <- gridsearch(series, args)
  r <- result$r
  q <- result$q
  p <- result$p

  thresh <- quantile(series[, args], q)

  # final declustering
  declustering <- decluster(series[, args], thresh, r = r)
  events   <- attr(declustering, "clusters")
  times    <- series$time[series[, args] > thresh]
  variable <- series[, args][series[, args] > thresh]
  metadata <- data.frame(time = times, event = events, variable = variable)

  # event stats
  events <- metadata |>
    group_by(event) |>
    mutate(
      event.size = n(),
      var_max = max(variable)
    ) |>
    filter(variable == var_max) |>
    group_by(event) |>
    summarise(
      variable = first(variable), # only for assigning rps
      #time = first(time),
      event.size = first(event.size)
    )

  # Ljung-box again
  p <- Box.test(c(events$variable), type = "Ljung")$p.value
  cat(paste0("Final Ljung-Box p-value: ", round(p, 4), '\n'))

  # event frequency
  m <- nrow(events)
  nyears <- length(unique(year(daily$time)))
  lambda <- m / nyears
  metadata$lambda <- lambda
  cat(paste0("Number of events: ", m, '\n'))

  # assign return periods
  survival_prob <- 1 - (
    rank(events$variable, ties.method = "average") / (m + 1)
  )
  rp <- 1 / (lambda * survival_prob)
  events$event.rp <- rp

  # remaining metadata
  metadata <- left_join(
    metadata,
    events[c("event", "event.rp", "event.size")],
    by = c("event")
  )
  metadata <- metadata |> rename_with(~ args, variable)

  # gridsearch fitted parameters
  metadata$run_length    <- rep(r, nrow(metadata))
  metadata$p_independent <- rep(p, nrow(metadata))
  metadata$event_thresh  <- rep(thresh, nrow(metadata))

  return(metadata)
}