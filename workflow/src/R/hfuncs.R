library(dplyr)
# Expected gridcell columns: event, event.rp, time, grid

hfunc_max <- function(gridcell, args) {
  res <- gridcell |>
    group_by(event) |>
    slice(which.max(get(args[1]))) |>
    ungroup() |>
    select(event, get(args[1]), time, event.rp, grid)
  return(res)
}


hfunc_sum <- function(gridcell, args) {
  stopifnot(c("event", "event.rp", "time", "grid") %in% names(gridcell))
  res <- gridcell |>
    group_by(event) |>
    summarise(
      variable = sum(get(args[1])),
      time = first(time),
      event.rp = max(event.rp),
      grid = first(grid),
      .groups = "drop"
    )
  return(res)
}

hfunc_mean <- function(gridcell, args) {
  stopifnot(c("event", "event.rp", "time", "grid") %in% names(gridcell))
  res <- gridcell |>
    group_by(event) |>
    summarise(
      variable = mean(get(args[1])),
      time = first(time),
      event.rp = max(event.rp),
      grid = first(grid),
      .groups = "drop"
    )
  return(res)
}


hfunc_arg2max <- function(gridcell, args) {
  stopifnot(c("event", "event.rp", "time", "grid") %in% names(gridcell))
  res <- gridcell |>
    group_by(event) |>
    slice(which.max(get(args[2]))) |>
    ungroup() |>
    select(event, get(args[1]), time, event.rp, grid)
  return(res)
}

hfunc_l2norm_argmax <- function(gridcell, args) {
  stopifnot(c("event", "event.rp", "time", "grid") %in% names(gridcell))
  
  l2norm <- function(u, v) {
    return(sqrt(u^2 + v^2))
  }
  gridcell$l2norm <- l2norm(gridcell[[args[1]]], gridcell[[args[2]]])

  res <- gridcell |>
    group_by(event) |>
    slice(which.max(l2norm)) |>
    ungroup() |>
    select(event, get(args[1]), time, event.rp, grid)
  return(res)
}