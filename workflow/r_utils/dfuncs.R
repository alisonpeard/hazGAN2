"Functions related to creating and declustering an event intensity
time series. Coles (2001) runs declustering."
library(extRemes)
library(methods) # for S4 methods / snakemake
# source("workflow/r_utils/utils.R")

`%ni%` <- Negate(`%in%`)

gridsearch <- function(series, var, qmin = 60, qmax = 99, rmin = 1, rmax = 14) {
  "Unit tests for this?"
  qvec <- c(qmin:qmax) / 100
  rvec <- c(rmin:rmax)
  
  print("Initial data summary:")
  print(summary(series[[var]]))
  
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
      
      # print(sprintf("r=%d, q=%.2f: ext_ind=%.3f",
      #               rvec[i], qvec[j], e[["extremal.index"]]))
      
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
  series <- aggregate(. ~ time, daily[, c("time", args)], rfunc)
  
  # gridsearch run lengths and thresholds
  result <- gridsearch(series, var)
  r <- result$r
  q <- result$q
  p <- result$p
  
  thresh <- quantile(series[[var]], q)
  cat(paste0(
    "Final selection from gridsearch:\n",
    "Run length: ", r, "\n",
    "Quantile: : ", q, "\n",
    "Threshold: ", round(thresh, 4), "\n",
    "P-value (H0:independent): ", round(p, 4), "\n"
  ))
  
  # final declustering
  declustering <- decluster(series[[var]], thresh, r = r)
  storms <- attr(declustering, "clusters")
  times <- series$time[series[[var]] > thresh]
  variable <- series[[var]][series[[var]] > thresh]
  metadata <- data.frame(time = times, storm = storms, variable = variable)
  
  # storm stats
  storms <- metadata |>
    group_by(storm) |>
    mutate(storm.size = n()) |>
    slice(which.max(variable)) |>
    summarise(
      variable = max(variable),
      time = time,
      storm.size = storm.size
    )
  
  # Ljung-box again
  p <- Box.test(c(storms$variable), type = "Ljung")$p.value
  cat(paste0("Final Ljung-Box p-value: ", round(p, 4), '\n'))
  
  # storm frequency
  m <- nrow(storms)
  nyears <- length(unique(year(daily$time)))
  lambda <- m / nyears
  metadata$lambda <- lambda
  cat(paste0("Number of storms: ", m, '\n'))
  
  # assign return periods
  survival_prob <- 1 - (
    rank(storms$variable, ties.method = "average") / (m + 1)
  )
  rp <- 1 / (lambda * survival_prob)
  storms$storm.rp <- rp
  
  # remaining metadata
  metadata <- left_join(metadata,
                        storms[c("storm", "storm.rp", "storm.size")],
                        by = c("storm"))
  metadata <- metadata |> rename_with(~ var, variable)
  
  return(metadata)
}