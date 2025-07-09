suppressPackageStartupMessages({
  library(extRemes)
  library(arrow)
  library(dplyr)
})

`%ni%` <- Negate(`%in%`)

CYCLONE_DATES_PATH <- "../resources/cyclones_midlands.csv"

identify_events <- function(daily, rfunc) {
  args    <- rfunc$args
  rfunc   <- match.fun(rfunc$func)

  # read text file with dates from ..resources/cyclones_midlands.txt
  event_data <- read.csv(CYCLONE_DATES_PATH)
  event_data <- event_data[event_data$wind > 0.0,]
  event_data$time <- as.Date(event_data$date)
  event_data$event <- as.numeric(factor(event_data$cyclone_id))

  metadata <- inner_join(event_data[c("time", "wind", "event")], daily, by = "time")
  times    <- metadata[, time]
  metadata$variable <- metadata[, args]
  variables <- metadata$variable
  metadata <- metadata[, c("time", "event", "variable")]

  # event stats
  events <- metadata |>
    group_by(event) |>
    mutate(
      event.size = n(),
      max_val = max(variable) #! max wind speed (hfunc hardcoded here)
    ) |>
    filter(variable == max_val) |>
    group_by(event) |>
    summarise(
      variable = first(variable),
      time = first(time),
      event.size = first(event.size)
    )

  # Ljung-Box again
  p <- Box.test(c(events$variable), type = "Ljung")$p.value
  cat(paste0("Final Ljung-Box p-value: ", round(p, 4), '\n'))

  # event frequency
  m <- nrow(events)
  nyears <- length(unique(event_data$year))
  lambda <- m / nyears
  metadata$lambda <- lambda
  cat(paste0("Number of events: ", m, '\n'))

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
  metadata <- metadata |> rename_with(~ args, variable)

  return(metadata)
}