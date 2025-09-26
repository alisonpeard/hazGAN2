# add any project-specific functions you want to source / overwrite here
# custom method for cyclone extraction from existing csv

suppressPackageStartupMessages({
  library(extRemes)
  library(arrow)
  library(dplyr)
})

`%ni%` <- Negate(`%in%`)

input_path <- "projects/poweruk_winter/resources/cyclones_midlands.csv"

identify_events <- function(daily, rfunc = NULL) {
  # read text file with dates from ..resources/cyclones_midlands.txt
  event_data <- read.csv(input_path)
  event_data <- event_data[event_data$wind > 0.0, ]
  event_data$time <- as.Date(event_data$date)
  event_data$year <- format(event_data$time, "%Y")
  event_data$event <- as.numeric(factor(event_data$cyclone_id))

  # assign ids to daily data... 
  # could also collapse daily along time here instead
  # of later...
  metadata <- inner_join(
    event_data[c("time", "wind", "event")],
    daily,
    by = "time"
  )
  if (is.null(metadata$event)) {
    stop("No events found in the data")
  }

  # clean up metadata
  metadata <- rename(metadata, variable = "wind")
  metadata <- metadata[, c("time", "event", "variable")]

  # extract event statistics
  events <- metadata |>
    group_by(time) |>
    summarise(
      event = first(event),
      variable = max(variable),
      .groups = "drop"
    ) |>
    group_by(event) |>
    summarise(
      event.size = n(),
      variable = max(variable),  # if you want to keep variable info
      .groups = "drop"
    )

  # Ljung-Box
  p <- Box.test(c(events$variable), type = "Ljung")$p.value
  cat(paste0(
    "Correlation between event max wind speeds (Ljung-Box p-value): ",
    round(p, 4),
    "\n"
  ))

  # get event frequency
  m <- nrow(events)
  nyears <- length(unique(event_data$year))
  lambda <- m / nyears
  metadata$lambda <- lambda
  cat(paste0("Number of events: ", m, "\n"))

  # assign return periods
  survival_prob <- 1 - (
    rank(events$variable, ties.method = "average") / (m + 1)
  )
  rp <- 1 / (lambda * survival_prob)
  events$event.rp <- rp

  # remaining metadata
  metadata <- left_join(
    metadata,
    events[c("event", "event.rp", "event.size")],
    by = c("event")
  )
  metadata <- metadata |> rename(wind = variable)
  metadata$p_independent <- rep(p, nrow(metadata))
  return(metadata)
}