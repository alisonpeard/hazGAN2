suppressPackageStartupMessages({
  library(arrow, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(logger, quietly = TRUE)
})

print("fit_marginals.R - Starting...")

source("workflow/src/R/fitting.R")
source("workflow/src/R/hfuncs.R")


if (!is.null(snakemake@params[["R_funcs"]])) {
  # overwrite with any custom functions defined in project/src/funcs.R
  rfunc_file <- snakemake@params[["R_funcs"]]
  source(rfunc_file)
}
print("Loaded custom R functions")

# configure logging
log_file <- snakemake@log[["file"]]
log_level <- snakemake@log[["level"]]
log_appender(appender_file(log_file))
log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
log_threshold(log_level)

# load snakemake rule paramet
METADATA <- snakemake@input[["metadata"]]
DAILY    <- snakemake@input[["daily"]]
OUTPUT   <- snakemake@output[["events"]]
FIELDS   <- snakemake@params[["fields"]]
Q        <- snakemake@params[["q"]]

# load input data
daily    <- read_parquet(DAILY)
metadata <- read_parquet(METADATA)

# get functions and args for resampling time
log_info("Transforming fields...")
fields  <- names(FIELDS)
distns  <- sapply(FIELDS, function(x) x$distn)
two_tailed <- sapply(FIELDS, function(x) x$two_tailed)
hfuncs  <- sapply(FIELDS, function(x) match.fun(paste0("hfunc_", x$hfunc$func)))
hfunc_args  <- sapply(FIELDS, function(x) x$hfunc$args)
nfields <- length(fields)


field_summary <- function(i) {
  paste0("Fitting field: ", fields[i], "\n")
}

# figure out expected output length: metadata
num_events <- length(unique(metadata$event))
num_x <- length(unique(daily$lon))
num_y <- length(unique(daily$lat))
num_fields <- length(fields)
expected_rows <- num_events * num_x * num_y * num_fields

# main: fit marginals and transform each field
log_debug(field_summary(1))
events_field1 <- marginal_transformer(
  daily, metadata, fields[1], Q,
  hfunc = hfuncs[[1]], hfunc_args = hfunc_args[[1]],
  distn = distns[[1]], two_tailed = two_tailed[[1]],
  log_file = log_file, log_level = log_level
)
log_info(paste0("Finished fitting: ", fields[1]))

log_debug(field_summary(2))
events_field2 <- marginal_transformer(
  daily, metadata, fields[2], Q,
  hfunc = hfuncs[[2]], hfunc_args = hfunc_args[[2]],
  distn = distns[[2]], two_tailed = two_tailed[[2]],
  log_file = log_file, log_level = log_level
)
log_info(paste0("Finished fitting: ", fields[2]))

log_debug(field_summary(3))
events_field3 <- marginal_transformer(
  daily, metadata, fields[3], Q,
  hfunc = hfuncs[[3]], hfunc_args = hfunc_args[[3]],
  distn = distns[[3]], two_tailed = two_tailed[[3]],
  log_file = log_file, log_level = log_level
)
log_info(paste0("fit_marginals.R - Finished fitting: ", fields[3]))

# combine fields
renamer <- function(df, var) {
  df <- df |>
    rename_with(~ paste0(., ".", var),
                -c("lat", "lon", "event", "event.rp", "variable"))
  df <- df |> rename_with(~ var, "variable")
  return(df)
}

log_success("fit_marginals.R - Done. Putting it all together...")

events_field1 <- renamer(events_field1, fields[1])
events_field2 <- renamer(events_field2, fields[2])
events_field3 <- renamer(events_field3, fields[3])

events <- events_field1 |>
  inner_join(events_field2, by = c("lat", "lon", "event", "event.rp")) |>
  inner_join(events_field3, by = c("lat", "lon", "event", "event.rp"))

events$thresh.q <- Q # keep track of threshold used

log_debug(paste0(
  "Combined events data has ", nrow(events), " rows."
))

# save results
log_info("fit_marginals.R - Saving results...")
write_parquet(events, OUTPUT)

nevents <- length(unique(events$event))
log_success(paste0(
  "fit_marginals.R - Finished! ", nevents, " events processed."
))
