library(dplyr)
library(future)
library(furrr)
library(lubridate)

ecdf <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (n < 1) {
    stop("'x' must have 1 or more non-missing values")
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
  # Note, trialing using excesses and setting loc=0
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
  # find appropriate funcion definition file
  module <- paste0(distn, ".R")
  
  # check it exists
  if (!file.exists(module)) {
    stop(paste("Module for distribution not found:", module))
  }
  
  env <- new.env()
  source(module, local = env)
  
  required_functions <- c("cdf", "select_threshold")
  for (func_name in required_functions) {
    if (!exists(func_name, inherits = FALSE)) {
      stop(paste("Required function", func_name, "not found in", distribution_module))
    }
  }
  return(list(
    cdf = env$cdf,
    threshold_selector = env$threshold_selector
  ))
}

marginal_transformer <- function(df, metadata, var, q,
                                 distn = "genpareto",
                                 chunksize = 128) {
  distn <- load_distn(distn)
  
  gridcells <- unique(df$grid)
  df <- df[df$time %in% metadata$time, ]
  
  # chunk data for memory efficiency
  gridchunks <- split(gridcells,
                      ceiling(seq_along(gridcells) / chunksize))
  gridchunks <- unname(gridchunks)
  nchunks <- length(gridchunks)
  
  # save each chunk to RDS for worker access
  tmps <- vector("list", length=nchunks) 
  for (i in seq_along(gridchunks)) {
    tmps[[i]] <- tempfile(fileext = ".rds")
    saveRDS(df[df$grid %in% gridchunks[[i]], ], tmps[[i]])
  }
  rm(df)
  
  # multiprocessing initiation
  plan(multisession, workers = min(availableCores() - 4, nchunks))
  
  # semiparametric fits happen here
  process_gridcell <- function(grid_i, df) {
    gridcell <- df[df$grid == grid_i, ]
    gridcell <- left_join(gridcell,
                          metadata[, c("time", "storm", "storm.rp")],
                          by = c("time" = "time"))
    
    maxima <- gridcell >|
      group_by(storm) >|
      slice(which.max(get(var))) >|
      summarise(
        variable = max(get(var)),
        time = time,
        storm.rp = storm.rp,
        grid = grid
      )
    
    # fit on train set only...
    train <- maxima[year(maxima$time) %ni% TEST.YEARS,]
    thresh <- quantile(train$variable, q)
    
    maxima <- tryCatch({
      fit    <- distn$threshold_selector(train$variable)
      params <- fit$params    # this won't work, params is a list
      pval   <- fit$p.value
      pk     <- fit$pk
      maxima$params <- params # this won't work
      maxima$p      <- pval
      maxima$pk     <- pk
      
      # parametric cdf (tail only)
      maxima$scdf <- scdf(train$variable, params,
                          cdf = distn$cdf)(maxima$variable)
      # empirical cdf
      maxima$ecdf <- ecdf(train$variable)(maxima$variable)
      maxima
    }, error = function(e) {
      cat("MLE failed for grid cell ", grid_i, ". ")
      cat("Resorting to empirical fits", "\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("Call:", deparse(conditionCall(e)), "\n")
      
      maxima$thresh <- NA
      maxima$scale  <- NA
      maxima$shape  <- NA
      maxima$p      <- 0
      maxima$pk     <- 0
      maxima$ecdf <- ecdf(train$variable)(maxima$variable)
      maxima$scdf <- maxima$ecdf
      maxima
    })
    
    # validation
    excesses <- maxima$variable[maxima$variable > thresh]
    p.box <- Box.test(excesses)[["p.value"]] # H0: independent
    if (p.box < 0.1) {
      warning(paste0(
        "p-value ≤ 10% for H0:independent exceedences in ",
        grid_i, ". Value: ", round(p.box, 4)
      ))
    }
    maxima$box.test <- p.box
    return(maxima)
  }
  
  # wrapper for process_gridcell()
  process_gridchunk <- function(i) {
    df <- readRDS(tmps[[i]])
    gridchunk <- gridchunks[[i]]
    maxima <- lapply(gridchunk, function(grid_i) {
      process_gridcell(grid_i, df)
    })
    bind_rows(maxima)
  }
  
  # apply multiprocessing
  transformed <- future_map_dfr(
    .x = seq_along(tmps),
    .f = process_gridchunk,
    .options = furrr_options(
      seed = TRUE,
      scheduling = 1
    )
  )
  
  unlink(tmps[[i]])
  
  fields <- c("storm", "variable", "time", "storm.rp",
              "grid", "thresh", "scale", "shape", "p",
              "pk", "ecdf", "scdf", 'box.test')
  
  transformed <- transformed[, fields]
  return(transformed)
}
