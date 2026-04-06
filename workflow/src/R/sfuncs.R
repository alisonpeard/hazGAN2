# add custom named sfuncs here
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(methods) # for S4 methods / snakemake
  library(data.table) # setDT
})


filter_months <- function(df, months_to_keep = c("October")) {
  df[, month_idx := month(time)]
  target_indices <- match(months_to_keep, month.name)
  df <- df[month_idx %in% target_indices]
  return(df)
}


winter <- function(df, vars) {
  winter_months <- c("December", "January", "February")
  df <- filter_months(df, winter_months)

  # monthly_median <- df[, .(m_val = median(get(vars), na.rm = TRUE)), 
  #                      by = .(month_idx, lat, lon)]
  monthly_median <- df[, lapply(.SD, median, na.rm = TRUE),
                       by = .(month_idx, lat, lon),
                       .SDcols = vars]
  setnames(monthly_median, vars, "m_val")

  df[monthly_median, on = .(month_idx, lat, lon), m_val := i.m_val]
  df[, (vars) := get(vars) - m_val]
  df[, m_val := NULL]

  monthly_median[, month := month.name[month_idx]]
  return(list(
    df = df[, c(vars, "time", "lat", "lon"), with = FALSE],
    params = monthly_median[, c("month", "lat", "lon", "m_val"), with = FALSE]
  ))
}


monthly_medians <- function(df, var) {
  df <- df[, c(var, "time", "lat", "lon")]
  df$month <- months(df$time)
  monthly_median <- aggregate(. ~ month + lat + lon,
                              df[, c(var, "lat", "lon", "month")],
                              median)
  return(monthly_median)
}


deseasonalise <- function(df, var, method = "additive") {
  setDT(df)
  method <- match.fun(method)
  return(method(df, var))
}


# additive <- function(df, vars) {
#   df$month <- months(df$time)
#   df <- df[, c(vars, "month", "lat", "lon")]
#   monthly_median <- aggregate(. ~ month + lat + lon, df, median)
#   df$monthly_median <- left_join(
#     df[, c("month", "lat", "lon")],
#     monthly_median,
#     by = c("month" = "month", "lat" = "lat", "lon" = "lon")
#   )[vars]
#   df[vars] <- df[vars] - df$monthly_median
#   return(list(
#     df = df[c(vars, "month", "lat", "lon")],
#     params = monthly_median
#   ))
# }