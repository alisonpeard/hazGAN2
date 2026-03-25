"Identify independent events using runs declustering."
suppressPackageStartupMessages({
  library(logger)
  library(arrow)
  library(lubridate)
  library(dplyr)
  require(ggplot2)
  library(tidync)
  library(purrr)
})

source("workflow/src/R/dfuncs.R")
source("workflow/src/R/sfuncs.R")
source("workflow/src/R/rfuncs.R")

if (!is.null(snakemake@params[["R_funcs"]])) {
  # overwrite with any custom functions defined in project/src/funcs.R
  rfunc_file <- snakemake@params[["R_funcs"]]
  source(rfunc_file)
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
RESX         <- snakemake@params[["resx"]]
RESY         <- snakemake@params[["resy"]]
XMIN         <- snakemake@params[["xmin"]]
XMAX         <- snakemake@params[["xmax"]]
YMIN         <- snakemake@params[["ymin"]]
YMAX         <- snakemake@params[["ymax"]]
FIELDS       <- snakemake@params[["fields"]]
RFUNC        <- snakemake@params[["rfunc"]]
SFUNC        <- snakemake@params[["sfunc"]]
FIELD_NAMES  <- names(FIELDS)

# ====== 1.LOAD AND DESEASONALIZE DATA =========================================
log_info("loading data")
src <- tidync(INPUT)
params_list <- list()
deseasonalized_all <- NULL

# deseasonalize each field
for (k in seq_along(FIELD_NAMES)) {
  field <- FIELD_NAMES[k]
  log_info(paste0("processing field: ", field))
  daily  <- src |> hyper_tibble(select_var = field, force = TRUE)
  log_info(paste0(
    "loaded object size: ",
    format(object.size(daily), units = "GB")
  ))

  daily <- daily[, c("lat", "lon", "time", field)]
  if (FIELDS[[k]]$obj == "min") {
    daily[field] <- (-1) * daily[field]
  }
  daily$time <- as.Date(daily$time)
  lats <- unique(daily$lat)
  lons <- unique(daily$lon)

  log_info(paste0("deseasonalizing field: ", field))
  deseasonalized <- deseasonalize(daily, field, method = SFUNC)

  log_info("Finished deseasonalizing, assigning parameters")
  params_list[[k]] <- deseasonalized$params[, c("month", "lat", "lon", field)]

  if (is.null(output)) {
    deseasonalized_all <- deseasonalized$df
  } else {
    deseasonalized_all <- output |>
      left_join(
        deseasonalized$df[, c("lat", "lon", "time", field)],
        by = c("lat", "lon", "time"),
        suffix = c("", paste0("_", field))
      )
  }
  rm(daily, deseasonalized); gc()
}

params <- Reduce(
  function(x, y) left_join(x, y, by = c("month", "lat", "lon")),
  params_list
)

log_info(paste0(
  "Finished deseasonalizing. Saving ",
  DAILY_OUT, " ..."
))

write_parquet(deseasonalized_all, DAILY_OUT)
write_parquet(params, MEDIANS_OUT)

stop("not finished debugging")

# ====== 2.EXTRACT EVENTS ======================================================
log_info("Extracting events")
metadata <- identify_events(daily, RFUNC)

log_info(paste0(
  "Finished event extraction. Saving ",
  METADATA_OUT, " ..."
))

write_parquet(metadata, METADATA_OUT)
