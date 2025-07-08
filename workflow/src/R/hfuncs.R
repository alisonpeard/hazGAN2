hfunc_max <- function(gridcell, var, vars){
  res <- gridcell |>
    group_by(event) |>
    slice(which.max(get(var))) |>
    summarise(
      variable = get(var),
      time = time,
      event.rp = event.rp,
      grid = grid
    )
  return(res)
}


hfunc_sum <- function(gridcell, var, vars){
  res <- gridcell |>
    group_by(event) |>
    summarise(
      variable = sum(get(var)),
      time = last(time),
      event.rp = first(event.rp),
      grid = first(grid),
      .groups = "drop"
    )
  return(res)
}