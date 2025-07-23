# add any project-specific functions you want to source / overwrite here
# custom method for cyclone extraction from existing csv

suppressPackageStartupMessages({
  library(extRemes)
  library(arrow)
  library(dplyr)
})

`%ni%` <- Negate(`%in%`)

input_path <- "../resources/cyclones_midlands.csv"

identify_events <- function(daily, rfunc = NULL) {

  # read text file with dates from ..resources/cyclones_midlands.txt
  event_data <- read.csv(input_path)
  event_data <- event_data[event_data$wind > 0.0, ]
  event_data$time <- as.Date(event_data$date)
  event_data$year <- format(event_data$time, "%Y")
  event_data$event <- as.numeric(factor(event_data$cyclone_id))

  # assign event ids to daily data
  metadata <- left_join(
    daily,
    event_data[c("time", "wind", "event")],
    by = "time"
  )

  # clean up metadata
  metadata$variable <- metadata[, "wind"]
  metadata <- metadata[, c("time", "event", "variable")]

  # extract event statistics
  events <- metadata |>
    group_by(event) |>
    mutate(
      event.size = n(),
    ) |>
    slice(which.max(variable)) |>
    ungroup() |>
    select(event, variable, time, event.size)

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
  metadata <- left_join(metadata,
                        events[c("event", "event.rp", "event.size")],
                        by = c("event"))
  metadata <- metadata |> rename(wind = variable)

  return(metadata)
}