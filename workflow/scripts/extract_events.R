"Identify independent events using runs declustering."
if (TRUE) { # all set up in this block
  # --- MANUAL DEBUGGING BLOCK ---
  if (!exists("snakemake")) {
    library(methods)
    # set wd to repository root to access src modules
    wd <- dirname(rstudioapi::getSourceEditorContext()$path)
    wd <- file.path(wd, "..", "..")
    setwd(wd)
    message("working directory set to: ", getwd())
    
    # create snakemake object
    setClass("Snakemake", slots = c(input="list", output="list", params="list", log="list"))
    snakemake <- new("Snakemake",
                     input  = list(netcdf = "/Users/alison/Documents/dphil/data/hazGAN2/projects/poweruk_winter/results/processing/resampled_all.parquet"),
                     output = list(medians = "~/Desktop/test_climatology.parquet", 
                                   metadata = "~/Desktop/test_event_metadata.parquet", 
                                   daily = "~/Desktop/event_cubes.parquet"),
                     params = list(
                       fields=list(
                         u10_gust=list(
                           init=list(
                             func="scale_to_gust",
                             args=c("u10", "v10", "i10fg")
                           ),
                           hfunc=list(
                             func="l2norm_argmax",
                             args=c("u10_gust", "v10_gust")
                           ),
                           obj="max"
                        ),
                         v10_gust=list(
                           init=list(
                             func="scale_to_gust",
                             args=c("u10", "v10", "i10fg")
                           ),
                           hfunc=list(
                             func="l2norm_argmax",
                             args=c("v10_gust", "u10_gust")
                           ),
                           obj="max"
                         ),
                        r30=list(
                          init=list(
                            func="sum_30_days",
                            args=c("tp")
                          ),
                          hfunc=list(
                            func="first",
                            args=c("r30")
                          ),
                          obj="max"
                        )),
                       rfunc="your_rfunc",
                       sfunc="winter",
                       R_funcs="/Users/alison/Documents/dphil/data/hazGAN2/projects/poweruk_winter/src/funcs.R"
                       ),
                     log= list(file = "~/Desktop/debug_extract_events.log")
    )
  }
  # ------------------------------
  
  suppressPackageStartupMessages({
    library(logger)
    library(arrow)
    library(lubridate)
    library(dplyr)
    require(ggplot2)
    library(tidync)
    library(purrr)
    library(duckdb)
  })
  
  source("workflow/src/R/dfuncs.R")
  source("workflow/src/R/sfuncs.R")
  source("workflow/src/R/rfuncs.R")
  
  if (!is.null(snakemake@params[["R_funcs"]])) {
    # overwrite with any custom functions defined in project/src/funcs.R
    rfunc_file <- snakemake@params[["R_funcs"]]
    source(rfunc_file)
  }
  
  # additional function for logging
  log_size <- function(obj) {
    log_info(paste0(
      "loaded object size: ",
      format(object.size(obj), units = "GB")
    ))
  }
  
  # configure logging
  log_appender(appender_file(snakemake@log[["file"]]))
  log_appender(appender_stdout(lines = 1))
  log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
  log_threshold(INFO)
  
  # load snakemake config
  log_info("Loading snakemake config...")
  INPUT        <- snakemake@input[["netcdf"]]
  MEDIANS_OUT  <- snakemake@output[["medians"]]
  METADATA_OUT <- snakemake@output[["metadata"]]
  DAILY_OUT    <- snakemake@output[["daily"]]
  # RESX         <- snakemake@params[["resx"]]
  # RESY         <- snakemake@params[["resy"]]
  # XMIN         <- snakemake@params[["xmin"]]
  # XMAX         <- snakemake@params[["xmax"]]
  # YMIN         <- snakemake@params[["ymin"]]
  # YMAX         <- snakemake@params[["ymax"]]
  FIELDS       <- snakemake@params[["fields"]]
  RFUNC        <- snakemake@params[["rfunc"]]
  SFUNC        <- snakemake@params[["sfunc"]]
  FIELD_NAMES  <- names(FIELDS)
} # keep all set up in here (collapse to focus on analysis)

# ====== 1.LOAD AND DESEASONALIZE DATA =========================================
con <- dbConnect(duckdb())
params_list <- list()
deseas <- NULL

month_filter <- switch(SFUNC,
  "winter" = c("'December','January','February'"),
  NULL
)

month_query <- if(!is.null(month_filter)) {
  paste0("WHERE monthname(time) IN (", paste(month_filter), ")")
} else ""

# loop through fields to remove seasonality
for (k in seq_along(FIELD_NAMES)) {
  field <- FIELD_NAMES[k]
  log_info(paste0("processing field: ", field))
  
  query <- paste0("SELECT lat, lon, time, ", field,
  " FROM '", INPUT, "' ", month_query)
  daily_k <- setDT(dbGetQuery(con, query)) # ~ 0.9 GB
  daily_k[, time := as.Date(time)]
  log_size(daily_k)

  if (FIELDS[[k]]$obj == "min") {
    daily_k[, (field) := get(field) * (-1)]
  }

  log_info(paste0("deseasonalising ...")) # from hfuncs.py or custom

  deseas_k <- deseasonalise(daily_k, field, method = SFUNC)

  log_info("finished deseasonalising, assigning parameters...")
  params_k <- deseas_k$params
  setnames(params_k, "m_val", field)
  params_list[[k]] <- params_k

  log_info("left-joining with final df...")
  setDT(deseas_k$df) # remind R that it's a data.table
  if (is.null(deseas)) {
    deseas <- copy(deseas_k$df)
  } else {
    setDT(deseas)
    deseas[deseas_k$df, on = .(lat, lon, time),
           (field) := get(paste0("i.", field))]
  }
  rm(daily_k, deseas_k)
  gc()
}

params <- Reduce(
  function(x, y) merge(x, y, by = c("month", "lat", "lon"), all = TRUE), 
  params_list
)

log_info(paste0(
  "Finished deseasonalising. Saving ",
  DAILY_OUT, " ..."
))

write_parquet(deseas, DAILY_OUT)
write_parquet(params, MEDIANS_OUT)

dbDisconnect(con)
stop("not finished debugging")

# ====== 2.EXTRACT EVENTS ======================================================
log_info("Extracting events")
metadata <- identify_events(deseas, RFUNC)

log_info(paste0(
  "Finished event extraction. Saving ",
  METADATA_OUT, " ..."
))

write_parquet(metadata, METADATA_OUT)
