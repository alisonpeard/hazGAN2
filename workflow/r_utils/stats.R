library(dplyr)
library(future)
library(furrr)
library(lubridate)
library(logger)

`%ni%` <- Negate(`%in%`)

ecdf <- function(x) {
  # modified to Weibull plotting positions
  x <- sort(x)
  n <- length(x)
  if (n < 1) {
    stop("'x' must have ≥1 non-missing values")
  }
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals))) / (n + 1),
                    method = "constant", #yleft = 0, yright = 1,
                    rule = 2, # take values at extremes
                    f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}

scdf <- function(train, params, cdf){
  # Using excesses and setting loc=0
  # This is for flexibility with cdf choice
  loc <- params$thresh
  calculator <- function(x){
    u <- ecdf(train)(x)
    pthresh <- ecdf(train)(loc)
    tail_mask <- x > loc
    x_tail <- x[tail_mask]
    exceedances <- x_tail - loc
    u_tail <- 1 - (1 - pthresh) * (1 - cdf(exceedances, params))
    u[tail_mask] <- u_tail
    return(u)
  }
  return(calculator)
}

load_distn <- function(distn) {
  # find appropriate function definition file
  module <- paste0("workflow/r_utils/", distn, ".R")

  # check it exists
  if (!file.exists(module)) {
    stop(paste("Module for distribution not found:", module))
  }

  env <- new.env()
  source(module, local = env)

  return(list(
    cdf = env$cdf,
    threshold_selector = env$threshold_selector
  ))
}

marginal_transformer <- function(df, metadata, var, q,
                                 distn = "genpareto",
                                 chunksize = 128,
                                 log_file = tempfile(fileext = ".log")) {
  # load functions for supplied distribution
  distn <- load_distn(distn)

  # only take days when an event is occurring
  df <- df[df$time %in% metadata$time, ]

  # chunk data for memory efficiency
  gridcells <- unique(df$grid)
  gridchunks <- split(
    gridcells, ceiling(seq_along(gridcells) / chunksize)
  )
  gridchunks <- unname(gridchunks)
  nchunks <- length(gridchunks)

  # save each chunk to RDS for worker access
  tmps <- vector("list", length = nchunks)
  for (i in seq_along(gridchunks)) {
    tmps[[i]] <- tempfile(fileext = ".rds")
    saveRDS(df[df$grid %in% gridchunks[[i]], ], tmps[[i]])
  }
  rm(df)

  # setup multiprocessing, test with `plan(sequential)`
  ncores <- max(1, min(availableCores(), nchunks))
  log_info(paste0("Using ", ncores, " cores"))
  plan(multisession, workers = ncores)

  # fit marginals using POT methods
  fit_gridcell <- function(grid_i, df) {
    # extract marginal
    gridcell <- df[df$grid == grid_i, ]
    gridcell <- left_join(
      gridcell,
      metadata[, c("time", "event", "event.rp")],
      by = c("time" = "time")
    )

    # extract maximum for each event
    maxima <- gridcell |>
      group_by(event) |>
      slice(which.max(get(var))) |>
      summarise(
        variable = max(get(var)),
        time = time,
        event.rp = event.rp,
        grid = grid
      )

    # omit test years from fitting functions
    if (exists("TEST.YEARS")) {
      train <- maxima[year(maxima$time) %ni% TEST.YEARS,]
    } else {
      train <- maxima
    }

    # main fitting functions happen here
    maxima <- tryCatch({
      # choose a threshold and fit parameters
      fit    <- distn$threshold_selector(train$variable)
      params <- fit$params
      pval   <- fit$p.value
      pk     <- fit$pk

      # add parameters and p-values to dataframe
      maxima    <- maxima |> mutate(!!!params)
      maxima$p  <- pval
      maxima$pk <- pk

      # transform variable to uniform
      maxima$scdf <- scdf(
        train$variable, params, cdf = distn$cdf
      )(maxima$variable)
      maxima$ecdf <- ecdf(train$variable)(maxima$variable)
      maxima
    }, error = function(e) {
      # log error
      calls <- sys.calls()
      calls_str <- paste(sapply(calls, deparse), collapse = "\n")
      msg <- paste0(
        "MLE failed for grid cell ", grid_i, "\n",
        "Error message: ", conditionMessage(e), "\n",
        "Call: ", deparse(conditionCall(e)), "\n",
        "Traceback:\n", calls_str, "\n"
      )

      log_error(skip_formatter(msg))

      # fallback to empirical fits
      maxima$thresh <- NA
      maxima$scale  <- NA
      maxima$shape  <- NA
      maxima$p      <- 0
      maxima$pk     <- 0
      maxima$ecdf <- ecdf(train$variable)(maxima$variable)
      maxima$scdf <- maxima$ecdf
      maxima
    })

    # validate exceedences are independent
    excesses <- maxima$variable[maxima$variable > maxima$thresh]
    nexcesses <- length(excesses)
    if (nexcesses < 30) {
      log_warn(paste0(
        "Only ", nexcesses,
        " exceedences in gridcell ",
        grid_i, ". ",
        "Skipping Ljung-Box test."
      ))
      p_box <- 0
    } else {
      # Box test H0: independent exceedences
      p_box <- Box.test(excesses)[["p.value"]]

      if (is.na(p_box)) {
        log_warn(paste0(
          "Ljung-Box test returned NA for gridcell ",
          grid_i, ". Value: ", round(p_box, 4)
        ))
      } else if (p_box < 0.1) {
        log_warn(paste0(
          "p-value ≤ 10% for H0: independent exceedences in ",
          grid_i, ". Value: ", round(p_box, 4)
        ))
      } else {
        log_info(paste0(
          "p-value > 10% for H0: independent exceedences in ",
          grid_i, ". Value: ", round(p_box, 4)
        ))
      }
    }

    maxima$box.test <- p_box
    return(maxima)
  }

  # wrapper for fit_gridcell()
  fit_gridchunk <- function(i) {
    # Reconfigure logging in workers
    library(logger)

    log_file <- as.character(log_file)
    log_appender(appender_file(log_file, append = TRUE))
    log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
    # log_layout(layout_format_generator(format = "%Y-%m-%d %H:%M:%S - [%l] - %m"))
    log_threshold(INFO)

    log_info(paste0("Fitting gridchunk ", i))

    df <- readRDS(tmps[[i]])
    gridchunk <- gridchunks[[i]]
    maxima <- lapply(gridchunk, function(grid_i) {
      fit_gridcell(grid_i, df)
    })
    bind_rows(maxima)
  }

  # fit gridcells with multiprocessing
  transformed <- future_map_dfr(
    .x = seq_along(tmps),
    .f = fit_gridchunk,
    .options = furrr_options(
      seed = TRUE,
      scheduling = 1
    )
  )

  unlink(tmps)

  # return transformed variable
  fields <- c("event", "variable", "time", "event.rp",
              "grid", "thresh", "scale", "shape", "p",
              "pk", "ecdf", "scdf", "box.test")

  transformed <- transformed[, fields]
  return(transformed)
}
