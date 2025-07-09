hfunc_max <- function(gridcell, args){
  res <- gridcell |>
    group_by(event) |>
    slice(which.max(get(args[1]))) |>
    summarise(
      variable = get(args[1]),
      time = time,
      event.rp = event.rp,
      grid = grid
    )
  return(res)
}


hfunc_sum <- function(gridcell, args){
  res <- gridcell |>
    group_by(event) |>
    summarise(
      variable = sum(get(args[1])),
      time = last(time),
      event.rp = first(event.rp),
      grid = first(grid),
      .groups = "drop"
    )
  return(res)
}


hfunc_arg2max <- function(gridcell, args){
  res <- gridcell |>
    group_by(event) |>
    slice(which.max(get(args[2]))) |>
    summarise(
      variable = get(args[1]),
      time = time,
      event.rp = event.rp,
      grid = grid
    )
  return(res)
}