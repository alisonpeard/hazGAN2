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
    stop("ecdf - 'x' must have ≥1 non-missing values")
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


scdf <- function(train, params, cdf, two_tailed = FALSE) {
    # wrapper for scdf_one_tail and scdf_two_tailed
    if (two_tailed) {
        return(scdf_two_tailed(train, params, cdf))
    } else {
        return(scdf_one_tail(train, params, cdf))
    }
}


scdf_one_tail <- function(train, params, cdf){
  # Using excesses and setting loc=0
  # This is for flexibility with cdf choice
  params_upper <- params[names(params) %in% c("thresh_upper", "scale_upper", "shape_upper")]
  loc <- params$thresh_upper
  calculator <- function(x) {
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


scdf_two_tailed <- function(train, params, cdf) {
  # Using excesses and setting loc=0
  # This is for flexibility with cdf choice
  params_upper <- params[names(params) %in% c("thresh_upper", "scale_upper", "shape_upper")]
  params_lower <- params[names(params) %in% c("thresh_lower", "scale_lower", "shape_lower")]

  loc_upper <- params$thresh_upper
  loc_lower <- params$thresh_lower

  calculator <- function(x) {
    u <- ecdf(train)(x)

    pthresh_upper <- ecdf(train)(loc_upper)
    mask_upper <- x > loc_upper
    x_upper <- x[mask_upper]
    exceedances_upper <- x_upper - mask_upper
    u_upper <- 1 - (1 - pthresh) * (1 - cdf(exceedances_upper, params_upper))
    u[mask_upper] <- u_upper

    pthresh_lower <- ecdf(train)(loc_lower)
    mask_lower <- x <= loc_lower
    x_lower <- x[mask_lower]
    exceedances_lower <- loc_lower - x_lower
    u_lower <- pthresh_lower * (1 - cdf(exceedances_lower, params_lower))
    u[mask_lower] <- u_lower

    return(u)
  }
  return(calculator)
}


load_distn <- function(distn) {
  # find appropriate function definition file
  module <- paste0("workflow/src/R/", distn, ".R")

  # check it exists
  if (!file.exists(module)) {
    stop(paste(
      "stats.R::load_distn - Module for distribution not found:", module
    ))
  }

  env <- new.env()
  source(module, local = env)

  return(list(
    cdf = env$cdf,
    threshold_selector = env$threshold_selector
  ))
}


format_stack_trace <- function(calls, max_depth = 10) {
  calls <- tail(calls, max_depth)
  formatted <- character(length(calls))
  for (i in seq_along(calls)) {
    formatted[i] <- paste0(i, ": ", deparse1(calls[[i]]))
  }
  paste(formatted, collapse = "\n  ")
}


log_fit_gridcell_error <- function(grid_i, error, max_frames = 10) {
  msg <- conditionMessage(error)
  stack_formatted <- sapply(stack, function(call) deparse1(call))

  error_report <- paste0(
    "stats.R::fit_gridcell - ",
    "ERROR IN GRID CELL ", grid_i, "\n",
    "===========================\n",
    "WARN: ", msg, "\n"
  )

  log_error(skip_formatter(error_report))
}


fit_gridcell <- function(
  grid_i, df, metadata,
  distn, two_tailed,
  hfunc, hfunc_args,
  log_file, log_level
  ) {
  # fit marginals using POT methods
  library(logger)
  log_appender(appender_file(log_file, append = TRUE))
  # log_appender(appender_tee(
  #   appender_file(log_file, append = TRUE),
  #   appender_console(lines = 1)
  # ))
  log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
  log_threshold(log_level)

  # extract marginal
  gridcell <- df[df$grid == grid_i, ]
  gridcell <- left_join(
    gridcell,
    metadata[, c("time", "event", "event.rp")],
    by = c("time" = "time")
  )

  footprint <- hfunc(gridcell, hfunc_args)
  footprint$variable <- footprint[[hfunc_args[1]]]

  # Check for no footprints early
  if (nrow(footprint) == 0) {
    log_warn(paste0(
      "stats.R::fit_gridcell - Empty footprint dataframe for grid cell ",
      grid_i
    ))
    # Create a 1-row dummy frame with all required columns
    return(data.frame(
      event = NA, variable = NA, time = NA, event.rp = NA,
      grid = grid_i, thresh = NA, scale = NA, shape = NA,
      p = 0, pk = 0, ecdf = NA, scdf = NA, box.test = 0
    ))
  }

  # omit test years from fitting functions
  if (exists("TEST_YEARS")) {
    train <- footprint[year(footprint$time) %ni% TEST_YEARS, ]
  } else {
    train <- footprint
  }

  # main fitting functions happen here
  footprint <- tryCatch({
    # choose a threshold and fit parameters
    fit    <- distn$threshold_selector(train$variable)
    params <- fit$params
    pval   <- fit$p.value
    pk     <- fit$pk

    # append "_upper" suffix to params
      params <- setNames(
          params, paste0(names(params), "_upper")
      )

    # add parameters and p-values to dataframe
    footprint          <- footprint |> mutate(!!!params)
    footprint$p_upper  <- pval
    footprint$pk_upper <- pk

    if (two_tailed) {
      # fit lower tail parameters
      fit_lower <- distn$threshold_selector(-train$variable)
      params_lower <- fit_lower$params
      pval_lower   <- fit_lower$p.value
      pk_lower     <- fit_lower$pk

      # append "_lower" suffix to params
      params_lower <- setNames(
        params_lower, paste0(names(params_lower), "_lower")
      )

      # add lower tail parameters and p-values to dataframe
      footprint <- footprint |> mutate(!!!params_lower)
      footprint$p_lower  <- pval_lower
      footprint$pk_lower <- pk_lower
    } else {
      # if not two-tailed, set lower tail params to NA
      footprint$thresh_lower <- NA
      footprint$scale_lower  <- NA
      footprint$shape_lower  <- NA
      footprint$p_lower      <- NA
      footprint$pk_lower     <- NA
    }

    # transform variable to uniform
    footprint$scdf <- scdf(
      train$variable, params, cdf = distn$cdf, two_tailed = two_tailed
    )(footprint$variable)

    footprint$ecdf <- ecdf(train$variable)(footprint$variable)
    footprint
  }, error = function(e) {
    # highest level error handling
    log_fit_gridcell_error(grid_i, e)

    # fallback to empirical fits
    footprint$thresh <- Inf
    footprint$scale  <- NA
    footprint$shape  <- NA
    footprint$p      <- 0
    footprint$pk     <- 0
    footprint$ecdf <- ecdf(train$variable)(footprint$variable)
    footprint$scdf <- footprint$ecdf
    footprint
  })

  # validate exceedences are independent
  excesses <- footprint$variable[footprint$variable > footprint$thresh]
  nexcesses <- length(excesses)
  if (nexcesses < 30) {
    log_warn(paste0(
      "stats.R::fit_gridcell - ",
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
        "stats.R::fit_gridcell - ",
        "Ljung-Box test returned NA for gridcell ",
        grid_i, ". Value: ", round(p_box, 4)
      ))
    } else if (p_box < 0.1) {
      log_warn(paste0(
        "stats.R::fit_gridcell - ",
        "p-value ≤ 10% for H0: independent exceedences in ",
        grid_i, ". Value: ", round(p_box, 4)
      ))
    } else {
      log_info(paste0(
        "stats.R::fit_gridcell - ",
        "p-value > 10% for H0: independent exceedences in ",
        grid_i, ". Value: ", round(p_box, 4)
      ))
    }
  }

  footprint$box.test <- p_box
  return(footprint)
}


marginal_transformer <- function(df, metadata, var, q,
                                 hfunc = "max",
                                 hfunc_args = NULL,
                                 distn = "genpareto",
                                 two_tailed = FALSE,
                                 chunksize = 128,
                                 log_file = tempfile(fileext = ".log"),
                                 log_level = INFO) {

  # reconfigure logging
  log_file <- as.character(log_file)
  log_appender(appender_file(log_file, append = TRUE))
  log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
  log_threshold(log_level)

  # load functions for specified extremal distribution
  distn <- load_distn(distn)

  # only take days when an event is occurring
  df <- df[df$time %in% metadata$time, ]

  # make a grid for easy chunking
  log_info("Creating grid for chunking")
  df$grid <- paste(df$lat, df$lon, sep = "_")
  df$grid <- as.integer(factor(df$grid))
  coords  <- df[df$time == min(df$time), ]
  coords  <- coords[, c("lat", "lon", "grid")]
  coords  <- coords[!duplicated(coords), ]

  # chunk data for memory efficiency
  log_info("Chunking data for memory efficiency")
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

  # Save metadata to file
  metadata <- metadata[, c("time", "event", "event.rp")]
  metadata_file <- tempfile(fileext = ".rds")
  saveRDS(metadata, metadata_file)
  rm(metadata)

  # setup multiprocessing
  ncores <- max(1, min(availableCores(), nchunks))
  plan(multisession, workers = ncores)
  log_debug(paste0("stats.R::marginal_transformer - Available cores: ", availableCores()))
  log_debug(paste0("stats.R::marginal_transformer - Number of chunks: ", nchunks))
  log_debug(paste0("stats.R::marginal_transformer - Using ", ncores, " cores"))

  fit_gridchunk <- function(i) {
    library(dplyr)
    library(lubridate)

    df <- readRDS(tmps[[i]])
    metadata <- readRDS(metadata_file)

    gridchunk <- gridchunks[[i]]
    maxima <- lapply(gridchunk, function(grid_i) {
      fit_gridcell(
        grid_i, df, metadata,
        distn, two_tailed,
        hfunc, hfunc_args,
        log_file, log_level
        )
    })
    bind_rows(maxima)
  }

  # fit gridcells with multiprocessing
  tryCatch({
    transformed <- future_map_dfr(
        .x = seq_along(tmps),
        .f = fit_gridchunk,
        .options = furrr_options(
          seed = TRUE,
          scheduling = 1,
          globals = list(
            tmps = tmps,
            metadata_file = metadata_file,
            gridchunks = gridchunks,
            fit_gridcell = fit_gridcell,
            distn = distn,
            two_tailed = two_tailed,
            hfunc = hfunc,
            hfunc_args = hfunc_args,
            log_file = log_file,
            log_level = log_level,
            log_fit_gridcell_error = log_fit_gridcell_error,
            ecdf = ecdf,
            scdf = scdf
          )
        )
      )
  }, error = function(e) {
    log_error(paste("Error in parallel processing:", conditionMessage(e)))
    stop(e)
  }, finally = {
    log_info("Cleaning up temporary files")
    unlink(tmps)
    unlink(metadata_file)
  })

  # map gridcell back to lat/lon
  log_info("Mapping gridcell back to lat/lon")
  transformed <- left_join(
    transformed,
    coords,
    by = c("grid" = "grid")
  )

  # return transformed variable
  fields <- c("event", "variable", "time", "event.rp",
              "lat", "lon",
              "thresh_upper", "scale_upper", "shape_upper", "p_upper", "pk_upper",
              "thresh_lower", "scale_lower", "shape_lower", "p_lower", "pk_lower",
              "ecdf", "scdf", "box.test")

  transformed <- transformed[, fields]
  return(transformed)
}
