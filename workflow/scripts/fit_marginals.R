library(arrow)
library(logger)
library(dplyr)

source("workflow/r_utils/stats.R")

# configure logging
log_file <- snakemake@log[["file"]]
print(paste0("Log file: ", log_file))
log_appender(appender_file(log_file))
log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
log_threshold(INFO)

# load snakemake rule parameters
METADATA <- snakemake@input[["metadata"]]
DAILY    <- snakemake@input[["daily"]]
OUTPUT   <- snakemake@output[["events"]]
FIELDS   <- snakemake@params[["fields"]]
Q        <- snakemake@params[["q"]]

# load input data
daily    <- read_parquet(DAILY)
metadata <- read_parquet(METADATA)

# transform marginals (hardcoded for three fields)
log_info("Tranforming fields...")
fields  <- names(FIELDS)
distns  <- sapply(FIELDS, function(x) x$distn)
nfields <- length(fields)

events_field1 <- marginal_transformer(
  daily, metadata, fields[1], Q, distn = distns[1],
  log_file = log_file
)
log_info(paste0("Finished fitting: ", fields[1]))
events_field2 <- marginal_transformer(
  daily, metadata, fields[2], Q, distn = distns[2],
  log_file = log_file
)
log_info(paste0("Finished fitting: ", fields[2]))
events_field3 <- marginal_transformer(
  daily, metadata, fields[3], Q, distn = distns[3],
  log_file = log_file
)
log_info(paste0("Finished fitting: ", fields[3]))

# combine fields
renamer <- function(df, var) {
  df <- df |>
    rename_with(~ paste0(., ".", var),
                -c("grid", "event", "event.rp", "variable"))
  df <- df |> rename_with(~ var, "variable")
  return(df)
}

log_info("Done. Putting it all together...")
events_field1 <- renamer(events_field1, fields[1])
events_field2 <- renamer(events_field2, fields[2])
events_field3 <- renamer(events_field3, fields[3])

events <- events_field1 |>
  inner_join(events_field2, by = c("grid", "event", "event.rp")) |>
  inner_join(events_field3, by = c("grid", "event", "event.rp"))

events$thresh.q <- Q # keep track of threshold used

# save results
log_info("Saving results...")
write_parquet(events, OUTPUT)

nevents <- nrow(events)
log_info(paste0("Finished! ", nevents, " events processed."))
