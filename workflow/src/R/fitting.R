suppressPackageStartupMessages({
  library(dplyr, quietly = TRUE)
  library(future, quietly = TRUE)
  library(furrr, quietly = TRUE)
  library(lubridate, quietly = TRUE)
  library(logger, quietly = TRUE)
})

source("workflow/src/R/stats.R")

`%ni%` <- Negate(`%in%`)

setup_logger <- function(log_file, log_level) {
  # setup logger in each function
  if (!is.null(log_file) && !is.null(log_level)) {
    log_appender(appender_file(log_file, append = TRUE))
    log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
    log_threshold(log_level)
  }
}

load_distn <- function(distn) {

  module <- paste0("workflow/src/R/", distn, ".R")

  if (!file.exists(module)) {
    stop(paste(
      "Module for distribution not found:", module
    ))
  }

  env <- new.env()
  source(module, local = env)

  return(list(
    cdf = env$cdf,
    threshold_selector = env$threshold_selector
  ))
}


fit_tails <- function(
  footprint, train, distn,
  two_tailed = FALSE, grid_i = NULL,
  log_file = NULL, log_level = NULL
) {

  setup_logger(log_file, log_level)

  upper_result <- tryCatch({

    fit <- distn$threshold_selector(train$variable, grid_i)
    params <- setNames(fit$params, paste0(names(fit$params), "_upper"))

    list(
      params = params,
      p = fit$p.value,
      pk = fit$pk,
      success = TRUE
    )
  }
  , error = function(e) {
    log_error(paste0(
      "Upper tail fitting failed for ", grid_i, ":\n", e$message
    ))
    list(
      params = list(
        thresh_upper = NA,
        scale_upper = NA,
        shape_upper = NA
      ),
      p = 0,
      pk = 0,
      success = FALSE
    )
  })

  # add upper tail results to footprint
  footprint <- footprint |> mutate(!!!upper_result$params)
  footprint$p_upper <- upper_result$p
  footprint$pk_upper <- upper_result$pk

  if (two_tailed) {
    lower_result <- tryCatch({
      fit <- distn$threshold_selector(-train$variable, grid_i)
      log_info(paste0(
        "Lower tail parameters for ", grid_i, ": ",
        paste(fit$params, collapse = ", ")
      ))

      params_lower <- setNames(
        fit$params, paste0(names(fit$params), "_lower")
      )

      list(
        params = params_lower,
        p = fit$p.value,
        pk = fit$pk,
        success = TRUE
      )
    }, error = function(e) {
      log_error(
        paste0(
          "Lower tail fitting failed for grid cell ",
          grid_i, ": ", e$message
        )
      )
      list(
        params = list(thresh_lower = NA, scale_lower = NA, shape_lower = NA),
        p = 0,
        pk = 0,
        success = FALSE
      )
    })
    # add lower tail results to footprint
    footprint <- footprint |> mutate(!!!lower_result$params)
    footprint$p_lower <- lower_result$p
    footprint$pk_lower <- lower_result$pk

  } else {
    # Set lower tail params to NA if not two-tailed
    footprint$thresh_lower <- NA
    footprint$scale_lower <- NA
    footprint$shape_lower <- NA
    footprint$p_lower <- NA
    footprint$pk_lower <- NA
  }

  # semiparametric CDF
  # if (upper_result$success || lower_result$success) {
  upper_params <- c("thresh_upper", "scale_upper", "shape_upper")
  lower_params <- c("thresh_lower", "scale_lower", "shape_lower")

  params <- c(
    footprint[upper_params],
    footprint[lower_params]
  )

  footprint$scdf <- scdf(
    train$variable, params, cdf = distn$cdf
  )(footprint$variable)

  # empirical CDF
  footprint$ecdf <- ecdf(train$variable)(footprint$variable)

  return(footprint)
}


fit_gridcell <- function(
  grid_i, df, metadata,
  distn, two_tailed,
  hfunc, hfunc_args,
  log_file = NULL, log_level = NULL
) {
  # fit marginals using POT methods
  setup_logger(log_file, log_level)

  # extract marginal
  gridcell <- df[df$grid == grid_i, ]
  gridcell <- left_join(
    gridcell,
    metadata[, c("time", "event", "event.rp")],
    by = c("time" = "time")
  )

  footprint <- hfunc(gridcell, hfunc_args)

  # Check for no footprints early
  if (nrow(footprint) == 0) {
    log_warn(paste0(
      "stats::fit_gridcell - Empty footprint dataframe for grid cell ",
      grid_i
    ))
    # Create a dummy frame
    unique_events <- gridcell |>
      select(event, event.rp, lat, lon) |>
      distinct()
    
    return(data.frame(
      event = unique_events$event,
      variable = NA,  # No meaningful variable value since hfunc failed
      time = NA,      # No meaningful time since hfunc failed  
      event.rp = unique_events$event.rp,
      lat = unique_events$lat,
      lon = unique_events$lon,
      grid = grid_i,
      thresh_upper = NA, scale_upper = NA, shape_upper = NA,
      p_upper = 0, pk_upper = 0, box.test.upper = 0,
      thresh_lower = NA, scale_lower = NA, shape_lower = NA,
      p_lower = 0, pk_lower = 0, box.test.lower = 0,
      ecdf = NA, scdf = NA
    ))
  }

  train <- footprint # legacy from previous train/test split

  # main fitting functions happen here
  footprint <- fit_tails(
    footprint, train, distn, two_tailed, grid_i, log_file, log_level
  )

  # validate exceedences are independent
  pbox_upper <- ljung_box(
    footprint$variable, footprint$thresh_upper, grid_i, "upper"
  )
  footprint$box.test.upper <- pbox_upper

  if (two_tailed) {
    pbox_lower <- ljung_box(
      -footprint$variable, -footprint$thresh_lower, grid_i, "lower"
    )
    footprint$box.test.lower <- pbox_lower
  } else {
    footprint$box.test.lower <- NA
  }

  return(footprint)
}


marginal_transformer <- function(df, metadata, var,
                                 hfunc = "max",
                                 hfunc_args = NULL,
                                 distn = "genpareto",
                                 two_tailed = FALSE,
                                 chunksize = 128,
                                 log_file = tempfile(fileext = ".log"),
                                 log_level = INFO) {

  log_info(paste0("Starting marginals transformation for ", var))

  # load functions for specified extremal distribution
  distn <- load_distn(distn)

  # only take days when an event is occurring
  df <- df[df$time %in% metadata$time, ]
  log_debug(paste0("Extracted ", nrow(df), " rows for variable ", var))

  # make a grid for easy chunking
  log_info("Creating grid for chunking")
  df$grid <- paste(df$lat, df$lon, sep = "_")
  df$grid <- as.integer(factor(df$grid))
  coords  <- df[df$time == min(df$time), ]
  coords  <- coords[, c("lat", "lon", "grid")]
  coords  <- coords[!duplicated(coords), ]

  # chunk data for memory efficiency
  log_info("Chunking data")
  gridcells <- unique(df$grid)
  gridchunks <- split(
    gridcells, ceiling(seq_along(gridcells) / chunksize)
  )
  gridchunks <- unname(gridchunks)
  nchunks <- length(gridchunks)
  log_debug(paste0("Chunksize: ", chunksize))

  # check no coordinates are lost
  gridchunk_sum <- 0
  for (chunk in gridchunks) {
    gridchunk_sum <- gridchunk_sum + length(chunk)
  }
  if (gridchunk_sum != length(gridcells)) {
    stop(paste0(
      "stats::marginal_transformer - ",
      "Chunking lost some grid cells. ",
      "Expected: ", length(gridcells), ", got: ", gridchunk_sum
    ))
  }
  rm(gridchunk_sum)

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
  log_debug(paste0("Available cores: ", availableCores()))
  log_debug(paste0("Number of chunks: ", nchunks, " and cores: ", ncores))

  fit_gridchunk <- function(i) {

    suppressPackageStartupMessages({
      library(dplyr)
      library(lubridate)
      library(logger)
    })

    log_appender(appender_file(log_file, append = TRUE))
    log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
    log_threshold(log_level)

    df <- readRDS(tmps[[i]]) # load subset of df for grid chunk
    meta <- readRDS(metadata_file)

    gridchunk <- gridchunks[[i]]

    footprints <- lapply(gridchunk, function(grid_i) {
      fit_gridcell(
        grid_i, df, meta,
        distn, two_tailed,
        hfunc, hfunc_args,
        log_file, log_level
      )
    })
    bind_rows(footprints)
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
          fit_tails = fit_tails,
          distn = distn,
          two_tailed = two_tailed,
          hfunc = hfunc,
          hfunc_args = hfunc_args,
          setup_logger = setup_logger,
          log_file = log_file,
          log_level = log_level,
          ecdf = ecdf,
          scdf = scdf,
          ljung_box = ljung_box
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
  log_debug(paste0("536: Mapped transformed has ", nrow(transformed), " rows."))

  # return transformed variable
  fields <- c("event", "variable", "time", "event.rp",
              "lat", "lon",
              "thresh_upper", "scale_upper", "shape_upper",
              "p_upper", "pk_upper", "box.test.upper",
              "thresh_lower", "scale_lower", "shape_lower",
              "p_lower", "pk_lower", "box.test.lower",
              "ecdf", "scdf")

  transformed <- transformed[, fields]

  return(transformed)
}