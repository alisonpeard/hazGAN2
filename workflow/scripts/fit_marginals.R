suppressPackageStartupMessages({
    library(arrow, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(logger, quietly = TRUE)
    source("workflow/r_utils/stats.R")
})

# configure logging
log_file <- snakemake@log[["file"]]
log_level <- snakemake@log[["level"]]
log_appender(appender_file(log_file))
log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
log_threshold(log_level)
log_debug(paste0("Log file: ", log_file))
log_debug(paste0("Logger level: ", log_level))


# load snakemake rule paramet
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

field_summary <- function(i) {
  summary_msg <- paste0(
    "fit_marginals.R - ",
    "Fitting field: ", fields[i], "\n",
    "Using distribution: ", distns[i], "\n",
    "Variable mean: ", mean(daily[[fields[i]]], na.rm = TRUE), "\n",
    "Variable sd: ", sd(daily[[fields[i]]], na.rm = TRUE), "\n",
    "Variable min: ", min(daily[[fields[i]]], na.rm = TRUE), "\n",
    "Variable max: ", max(daily[[fields[i]]], na.rm = TRUE), "\n",
    "Q70: ", quantile(daily[[fields[i]]], 0.7, na.rm = TRUE), "\n",
    "Q95: ", quantile(daily[[fields[i]]], 0.95, na.rm = TRUE), "\n",
    "Q99: ", quantile(daily[[fields[i]]], 0.99, na.rm = TRUE), "\n"
  )
}

log_debug(field_summary(1))
events_field1 <- marginal_transformer(
  daily, metadata, fields[1], Q, distn = distns[1],
  log_file = log_file, log_level = log_level
)
log_info(paste0("Finished fitting: ", fields[1]))

log_debug(field_summary(2))
events_field2 <- marginal_transformer(
  daily, metadata, fields[2], Q, distn = distns[2],
  log_file = log_file, log_level = log_level
)
log_info(paste0("Finished fitting: ", fields[2]))

log_debug(field_summary(3))
events_field3 <- marginal_transformer(
  daily, metadata, fields[3], Q, distn = distns[3],
  log_file = log_file, log_level = log_level
)
log_info(paste0("fit_marginals.R - Finished fitting: ", fields[3]))

# combine fields
renamer <- function(df, var) {
 df <- df |>
    rename_with(~ paste0(., ".", var),
                -c("grid", "event", "event.rp", "variable"))
  df <- df |> rename_with(~ var, "variable")
  return(df)
}

log_success("fit_marginals.R - Done. Putting it all together...")
events_field1 <- renamer(events_field1, fields[1])
events_field2 <- renamer(events_field2, fields[2])
events_field3 <- renamer(events_field3, fields[3])

events <- events_field1 |>
  inner_join(events_field2, by = c("grid", "event", "event.rp")) |>
  inner_join(events_field3, by = c("grid", "event", "event.rp"))

events$thresh.q <- Q # keep track of threshold used

# save results
log_info("fit_marginals.R - Saving results...")
write_parquet(events, OUTPUT)

nevents <- length(unique(events$event))
log_success(paste0("fit_marginals.R - Finished! ", nevents, " events processed."))
