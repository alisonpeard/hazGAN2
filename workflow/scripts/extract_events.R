# Identify independent events using runs declustering.
suppressPackageStartupMessages({
  library(logger)
  library(arrow)
  library(lubridate)
  library(dplyr)
  require(ggplot2)
  library(tidync)
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
log_info("Loading data")
start <- Sys.time()
src <- tidync(INPUT)

daily  <- src |> hyper_tibble(force = TRUE)

# negate any variables where minimum is of interest
log_info("Negating any variables where minimum is of interest")
daily <- daily[, c("lat", "lon", "time", FIELD_NAMES)]
for (k in seq_along(FIELDS)) {
  if (FIELDS[[k]]$obj == "min") {
    daily[[FIELD_NAMES[k]]] <- (-1) * daily[[FIELD_NAMES[k]]]
  }
}

# additional processing
daily$time <- as.Date(daily$time)

# remove seasonality (sfuncs)
log_info("Removing seasonality")
log_info(paste0("Dataframe rows: ", nrow(daily)))
log_info(paste0("Dataframe columns: ", ncol(daily)))

# faster (hopefully)
lats <- unique(daily$lat)
lons <- unique(daily$lon)
latlon <- data.frame(expand.grid(lat = lats, lon = lons))
months <- month.name
params <- latlon[rep(seq_len(nrow(latlon)), each = length(months)), ]
params$month <- rep(months, nrow(latlon))

# deseasonalize each field
for (k in seq_along(FIELD_NAMES)) {
  field <- FIELD_NAMES[k]
  log_info(paste0("Deseasonalizing field: ", field))

  deseasonalized <- deseasonalize(daily, field, method = SFUNC)
  log_info("Finished deseasonalizing, assigning parameters")

  # daily[, field] <- deseasonalized$df
  daily <- deseasonalized$df |>
    left_join(daily, by = c("lat", "lon", "time"), suffix = c("", "_old")) |>
    rename(!!field := !!sym(field))
  
  params[, field] <- left_join(
    params[, c("month", "lat", "lon")],
    deseasonalized$params,
    by = c("month" = "month", "lat" = "lat", "lon" = "lon")
  )[field]
}

# ====== 2.EXTRACT EVENTS ======================================================
log_info("Extracting events")
metadata <- identify_events(daily, RFUNC)

# ====== 3.SAVE RESULTS ========================================================
log_info(paste0("Finished event extraction. Saving ",
                METADATA_OUT, ", and ", DAILY_OUT, " ..."))
write_parquet(params, MEDIANS_OUT)
write_parquet(daily, DAILY_OUT)
write_parquet(metadata, METADATA_OUT)

# ==============================================================================
