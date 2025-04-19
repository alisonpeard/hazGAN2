rm(list = ls())
library(arrow)
library(logger)
library(dplyr)

source("workflow/r_utils/stats.R")

# configure logging
log_appender(appender_file(snakemake@log[[1]]))
log_layout(layout_glue("{time} - {level} - {msg}"))
log_threshold(INFO)

# load snakemake rule parameters
METADATA <- snakemake@input[["metadata"]]
DAILY    <- snakemake@input[["daily"]]
OUTPUT   <- snakemake@output[["storms"]]
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

storms_field1 <- marginal_transformer(
  daily, metadata, fields[1], Q, distn = distns[1]
)
storms_field2 <- marginal_transformer(
  daily, metadata, fields[2], Q, distn = distns[2]
)
storms_field3 <- marginal_transformer(
  daily, metadata, fields[3], Q, distn = distns[3]
)

# combine fields
renamer <- function(df, var) {
  df <- df |>
    rename_with(~ paste0(., ".", var),
                -c("grid", "event", "event.rp", "variable"))
  df <- df |> rename_with(~ var, "variable")
  return(df)
}

log_info("Done. Putting it all together...")
storms_field1 <- renamer(storms_field1, fields[1])
storms_field2 <- renamer(storms_field2, fields[2])
storms_field3 <- renamer(storms_field3, fields[3])

storms <- storms_field1 |>
  inner_join(storms_field2, by = c("grid", "event", "event.rp")) |>
  inner_join(storms_field3, by = c("grid", "event", "event.rp"))

storms$thresh.q <- Q # keep track of threshold used

# save results
log_info("Saving results...")
write_parquet(storms, OUTPUT)

nevents <- nrow(storms)
log_info(paste0("Finished! ", nevents, " events processed."))
