suppressPackageStartupMessages({
  library(arrow, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(logger, quietly = TRUE)
})

source("workflow/src/R/fitting.R")
source("workflow/src/R/hfuncs.R")

if (!is.null(snakemake@params[["R_funcs"]])) {
  # overwrite with any custom functions defined in project/src/funcs.R
  rfunc_file <- snakemake@params[["R_funcs"]]
  source(rfunc_file)
  print(paste0("loaded custom R functions from ", rfunc_file))
}

# configure logging
log_file <- snakemake@log[["file"]]
log_level <- snakemake@log[["level"]]
log_appender(appender_file(log_file))
log_layout(layout_glue_generator(format = "{time} - {level} - {msg}"))
log_threshold(log_level)

# load snakemake rule parameters
METADATA <- snakemake@input[["metadata"]]
CUBES    <- snakemake@input[["cubes"]]
OUTPUT   <- snakemake@output[["events"]]
FIELDS   <- snakemake@params[["fields"]]

# load input data
log_info(paste0("memory before cubes: ", gc()[[2, 2]], " MB"))
cubes    <- read_parquet(CUBES) # 1.11 GB?
log_info(paste0("memory after cubes and before metadata: ", gc()[[2, 2]], " MB"))
metadata <- read_parquet(METADATA) # 1.1 MB?
log_info(paste0("memory after metadata: ", gc()[[2, 2]], " MB"))

log_info(paste0("cubes dims: ", nrow(cubes), " x ", ncol(cubes)))
log_info(paste0("metadata dims: ", nrow(metadata), " x ", ncol(metadata)))

stop("deliberate early exit")

# get functions and args for resampling time dimension
fields  <- names(FIELDS)
distns  <- sapply(FIELDS, function(x) x$distn)
two_tailed <- sapply(FIELDS, function(x) x$two_tailed)
hfuncs  <- sapply(FIELDS, function(x) match.fun(paste0("hfunc_", x$hfunc$func)))
hfunc_args  <- sapply(FIELDS, function(x) x$hfunc$args)
nfields <- length(fields)


field_summary_msg <- function(i) {
  paste0("Fitting field: ", fields[i], "\n")
}

# main: fit marginals and transform each field
log_info("\ntransforming fields...")

log_debug(field_summary_msg(1))
events_field1 <- marginal_transformer(
  cubes, metadata, fields[1],
  hfunc = hfuncs[[1]], hfunc_args = hfunc_args[[1]],
  distn = distns[[1]], two_tailed = two_tailed[[1]],
  logfile = log_file, loglevel = log_level
)
log_info(paste0("Finished fitting: ", fields[1]))

log_debug(field_summary_msg(2))
events_field2 <- marginal_transformer(
  cubes, metadata, fields[2],
  hfunc = hfuncs[[2]], hfunc_args = hfunc_args[[2]],
  distn = distns[[2]], two_tailed = two_tailed[[2]],
  logfile = log_file, loglevel = log_level
)
log_info(paste0("Finished fitting: ", fields[2]))

log_debug(field_summary_msg(3))
events_field3 <- marginal_transformer(
  cubes, metadata, fields[3],
  hfunc = hfuncs[[3]], hfunc_args = hfunc_args[[3]],
  distn = distns[[3]], two_tailed = two_tailed[[3]],
  logfile = log_file, loglevel = log_level
)
log_info(paste0("fit_marginals.R - Finished fitting: ", fields[3]))
log_success("fit_marginals.R - Done. Putting it all together...")

# combine fields
renamer <- function(df, var) {
  df <- df |>
    rename_with(
      ~ paste0(., ".", var), -c("lat", "lon", "event", "event.rp", "variable"
    ))
  df <- df |> rename_with(~ var, "variable")
  return(df)
}

events_field1 <- renamer(events_field1, fields[1])
events_field2 <- renamer(events_field2, fields[2])
events_field3 <- renamer(events_field3, fields[3])

events <- events_field1 |>
  inner_join(events_field2, by = c("lat", "lon", "event", "event.rp")) |>
  inner_join(events_field3, by = c("lat", "lon", "event", "event.rp"))

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
