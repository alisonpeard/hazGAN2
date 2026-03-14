suppressPackageStartupMessages({
  library(dplyr, quietly = TRUE)
  library(future, quietly = TRUE)
  library(furrr, quietly = TRUE)
  library(lubridate, quietly = TRUE)
  library(logger, quietly = TRUE)
})

source("workflow/src/R/stats.R") # provides: ecdf_wb, scdf_wb, ljung_box


setup_logger <- function(logfile, loglevel) {
  # setup logger in each function
  if (!is.null(logfile) && !is.null(loglevel)) {
    log_appender(appender_file(logfile, append = TRUE))
    log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
    log_threshold(loglevel)
  }
}

load_distn <- function(distn) {
  module <- paste0("workflow/src/R/", distn, ".R")
  if (!file.exists(module)) {
    stop(paste0(
      "No R module found for distribution: '", module
    ))
  }
  env <- new.env()
  source(module, local = env)
  return(list(
    cdf = env$cdf,
    threshold_selector = env$threshold_selector
  ))
}


fit_tail <- function(x, distn, grid_label, label) {
  tryCatch({
    fit <- distn$threshold_selector(x, grid_label)
    list(
      params = setNames(fit$params, paste0(names(fit$params), label)),
      p = fit$p.value,
      pk = fit$pk,
      success = TRUE
    )
  }, error = function(e) {
    log_error(paste0(
      "Tail fit failed for ", grid_label, label, ":\n", e$message
    ))
    list(
      params = setNames(
        list(NA, NA, NA),
        paste0(c("thresh", "scale", "shape"), label)
      ),
      p = 0,
      pk = 0,
      success = FALSE
    )
  })
}


fit_margin_tails <- function(
  margin, distn,
  two_tailed = FALSE,
  grid_i = NULL
) {

  # fit semiparametric distribution to upper tails
  upper_result <- fit_tail(margin$variable, distn, grid_i, "_upper")
  margin <- margin |> mutate(!!!upper_result$params)
  margin$p_upper <- upper_result$p
  margin$pk_upper <- upper_result$pk

  # (optional) fit lower tails of margin too
  if (two_tailed) {
    lower_result <- fit_tail(-margin$variable, distn, grid_i, "_lower")
    lower_result$params$thresh_lower <- -lower_result$params$thresh_lower
    margin <- margin |> mutate(!!!lower_result$params)
    margin$p_lower <- lower_result$p
    margin$pk_lower <- lower_result$pk

  } else {
    # Set lower tail params to NA if not two-tailed
    lower_result <- list(
      params = list(
        thresh_lower = NA,
        scale_lower = NA,
        shape_lower = NA
      )
    )
    margin$thresh_lower <- NA
    margin$scale_lower  <- NA
    margin$shape_lower  <- NA
    margin$p_lower      <- NA
    margin$pk_lower     <- NA
  }

  # apply empirical and semiparametric distribution transforms
  params <- c(upper_result$params, lower_result$params)
  margin$scdf <- scdf_wb(margin$variable, params, cdf = distn$cdf)(margin$variable)
  margin$ecdf <- ecdf_wb(margin$variable)(margin$variable)

  return(margin)
}


fit_margin <- function(
  grid_i, df, metadata,
  distn, two_tailed,
  hfunc, hfunc_args
) {
  # extract marginal
  gridcell <- df[df$grid == grid_i, ]
  gridcell <- left_join(
    gridcell,
    metadata[, c("time", "event", "event.rp")],
    by = c("time" = "time")
  )

  margin <- hfunc(gridcell, hfunc_args)

  # Return empty if no data
  if (nrow(margin) == 0) {
    log_warn(paste0(
      "Empty footprint dataframe for grid cell ", grid_i
    ))
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

  # main fitting functions happen here
  margin <- fit_margin_tails(margin, distn, two_tailed, grid_i)

  # validate exceedences are independent
  pbox_upper <- ljung_box(
    margin$variable, margin$thresh_upper, grid_i, "upper"
  )
  margin$box.test.upper <- pbox_upper

  if (two_tailed) {
    pbox_lower <- ljung_box(
      -margin$variable, -margin$thresh_lower, grid_i, "lower"
    )
    margin$box.test.lower <- pbox_lower
  } else {
    margin$box.test.lower <- NA
  }

  return(margin)
}


fit_chunk <- function(
  i, tmps, metadata_file, gridchunks,
  distn, two_tailed, hfunc, hfunc_args,
  logfile, loglevel
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(lubridate)
    library(logger)
  })

  setup_logger(logfile, loglevel)

  df <- readRDS(tmps[[i]]) # load subset of df for grid chunk
  meta <- readRDS(metadata_file)

  gridchunk <- gridchunks[[i]]

  footprints <- lapply(gridchunk, function(grid_i) {
    fit_margin(
      grid_i, df, meta,
      distn, two_tailed,
      hfunc, hfunc_args
    )
  })
  bind_rows(footprints)
}


marginal_transformer <- function(df, metadata, var,
                                 hfunc = "max",
                                 hfunc_args = NULL,
                                 distn = "genpareto",
                                 two_tailed = FALSE,
                                 chunksize = 128,
                                 logfile = tempfile(fileext = ".log"),
                                 loglevel = INFO) {

  log_info(paste0("Starting marginals transformation for ", var))

  # load functions for specified extremal distribution
  distn <- load_distn(distn)
  hfunc <- if (is.character(hfunc)) match.fun(hfunc) else hfunc
  force(hfunc)

  # only take days when an event is occurring
  df <- df[df$time %in% metadata$time, ]
  log_debug(paste0("Extracted ", nrow(df), " rows for variable ", var))

  # make a grid for easy chunking
  df$grid <- paste(df$lat, df$lon, sep = "_")
  df$grid <- as.integer(factor(df$grid))
  coords <- unique(df[, c("lat", "lon", "grid")])

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
  stopifnot(sum(lengths(gridchunks)) == length(gridcells))

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

  # fit gridcells with multiprocessing
  transformed <- tryCatch({
    future_map_dfr(
      .x = seq_along(tmps),
      .f = fit_chunk,
      tmps, metadata_file, gridchunks,
      distn, two_tailed, hfunc, hfunc_args,
      logfile, loglevel,
      .options = furrr_options(
        seed = TRUE,
        scheduling = 1,
        globals = list(
          fit_margin = fit_margin,
          fit_margin_tails = fit_margin_tails,
          fit_tail = fit_tail,
          setup_logger = setup_logger,
          ecdf_wb = ecdf_wb,    # Weibull plotting positions (stats.R)
          scdf_wb = scdf_wb,
          ljung_box = ljung_box,
          hfunc = hfunc
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
  transformed <- left_join(transformed, coords, by = "grid")
  log_debug(paste0("Mapped transformed has ", nrow(transformed), " rows."))

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