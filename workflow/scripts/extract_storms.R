"Identify independent storm events using runs declustering."
rm(list = ls())
library(logger)
library(arrow)
library(lubridate)
library(dplyr)
require(ggplot2)
library(CFtime)
library(tidync)

source("workflow/r_utils/dfuncs")
source("workflow/r_utils/sfuncs")

# configure logging
log_appender(appender_file(snakemake@log[[1]]))
log_layout(layout_glue("{time} - {level} - {msg}"))
log_threshold(INFO)

# load snakemake config
INPUT        <- snakemake@input[['netcdf']]
MEDIANS_OUT  <- snakemake@output[['medians']]
METADATA_OUT <- snakemake@output[['metadata']]
DAILY_OUT    <- snakemake@output[['daily']]
RESX         <- snakemake@params[['resx']]
RESY         <- snakemake@params[['resy']]
XMIN         <- snakemake@params[['xmin']]
XMAX         <- snakemake@params[['xmax']]
YMIN         <- snakemake@params[['ymin']]
YMAX         <- snakemake@params[['ymax']]
FIELDS       <- snakemake@params[['fields']]
RFUNC        <- snakemake@params[['rfunc']]
SFUNC        <- snakemake@params[['sfunc']]
FIELD_NAMES  <- names(FIELDS)

# ====== 1.LOAD AND DESEASONALIZE DATA =========================================
log_info("Loading data and removing seasonality...")
start <- Sys.time()
src <- tidync(INPUT)

daily  <- src |> hyper_tibble(force = TRUE)
coords <- src |> activate("grid") |> hyper_tibble(force = TRUE)
daily  <- left_join(daily, coords, by = c("lon", "lat"))

rm(coords)

# negate any variables where minimum is of interest
daily <- daily[, cbind(c("grid", "time"), FIELD_NAMES)]
for (k in seq_along(FIELDS)) {
  if (FIELDS[k]$obj == "min") {
    daily[[FIELD_NAMES[k]]] <- (-1) * daily[[FIELD_NAMES[k]]]
  }
}

daily$time <- as.Date(STARTDATE) + days(daily$time)
daily$grid <- as.integer(daily$grid)

# remove seasonality (sfuncs)
for (k in seq_along(FIELD_NAMES)) {
  field <- FIELD_NAMES[k]
  daily[[field]] <- deseasonalize(daily, field, method = SFUNC)
}

# ====== 2.EXTRACT EVENTS ======================================================
log_info("Extracting events") # (dfuncs)
metadata <- identify_events(daily, RFUNC)

# ====== 3.SAVE RESULTS ========================================================
log_info(paste0("Finished event extraction. Saving ",
                MEDIANS_OUT, ", ", METADATA_OUT, ", and ", DAILY_OUT, " ..."))
write.csv(medians, MEDIANS_OUT, row.names = FALSE)
write_parquet(metadata, METADATA_OUT)
write_parquet(daily, DAILY_OUT)

# ==============================================================================
