"Identify independent storm events.

Includes event metadata such as return period and frequency. Only fit ECDF and 
GPD to 1940-2022 data. Keep 2021 as a holdout test set.

Files:
------
- input:
  - data_1941_2022.nc

- output:
  - storms.parquet
  - metadata.parquet
  - medians.csv

LAST RUN: 17-02-2025
"
#%%######### START #############################################################
# Clear environment
rm(list = ls())

library(arrow)
library(lubridate)
library(dplyr)
require(ggplot2)
library(CFtime)
library(tidync)
source("r_utils/utils.R")

# configure logging
log_appender(appender_file(snakemake@log[[1]]))
log_layout(layout_glue("{time} - {level} - {msg}"))
log_threshold(INFO)

# load snakemake config
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

#%%######### LOAD AND STANDARDISE DATA #########################################
log_info("Loading and standardising data...")
start <- Sys.time()
src <- tidync(INPUT)

daily  <- src %>% hyper_tibble(force = TRUE)
coords <- src %>% activate("grid") %>% hyper_tibble(force = TRUE)
daily  <- left_join(daily, coords, by = c("lon", "lat"))

rm(coords)

# ensure maxising everything from now on
daily <- daily[, cbind(c("grid", "time"), FIELD_NAMES)]
for (k in seq_along(FIELDS)) {
  if (FIELDS[k]$obj == "min") {
    daily[[FIELD_NAMES[k]]] <- (-1) * daily[[FIELD_NAMES[k]]]
  }
}

daily$time <- as.Date(STARTDATE) + days(daily$time)
daily$grid <- as.integer(daily$grid)

for (k in seq_along(FIELD_NAMES)) {
  field <- FIELD_NAMES[k]
  daily[[field]] <- deseasonalize(daily, field, method = SFUNC)
}

#%%######## EXTRACT STORMS #####################################################
log_info("Extracting storms...")
metadata <- storm_extractor(daily, RFUNC)

########### SAVE RESULTS #######################################################
log_info(paste0("Finished storm extraction. Saving ",
                MEDIANS_OUT, ", ", METADATA_OUT, ", and ", DAILY_OUT, " ..."))
write.csv(medians, MEDIANS_OUT, row.names = FALSE)
write_parquet(metadata, METADATA_OUT)
write_parquet(daily, DAILY_OUT)

########### END ################################################################
