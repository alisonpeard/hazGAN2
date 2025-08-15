# add custom named sfuncs here
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(methods) # for S4 methods / snakemake
})


deseasonalize <- function(df, var, method = "additive") {
  # add more methods below, var is a list of variables
  method <- match.fun(method)
  if (is.null(method)) {
    stop(paste0("Unrecognized deseasonalization method: ", method))
  }
  return(method(df, var))
}


additive <- function(df, vars) {
  df$month <- months(df$time)
  df <- df[, c(vars, "month", "lat", "lon")]
  monthly_median <- aggregate(. ~ month + lat + lon, df, median)
  df$monthly_median <- left_join(
    df[, c("month", "lat", "lon")],
    monthly_median,
    by = c("month" = "month", "lat" = "lat", "lon" = "lon")
  )[vars]
  df[vars] <- df[vars] - df$monthly_median
  return(list(
    df = df[c(vars, "month", "lat", "lon")],
    params = monthly_median
  ))
}


autumn <- function(df, vars) {
  # extract autumn months only
  autumn_months <- c("September", "October", "November")
  df <- filter_months(df, autumn_months)

  monthly_median <- df[, c(vars, "month", "lat", "lon")]
  monthly_median <- aggregate(. ~ month + lat + lon, monthly_median, median)

  df$monthly_median <- left_join(
    df[, c("month", "lat", "lon")],
    monthly_median,
    by = c("month" = "month", "lat" = "lat", "lon" = "lon")
  )[vars]
  df[vars] <- df[vars] - df$monthly_median

  df <- df[, c(vars, "time", "lat", "lon")]

  return(list(
    df = df,
    params = monthly_median
  ))
}


winter <- function(df, vars) {
  # extract winter months only
  winter_months <- c("December", "January", "February")
  df <- filter_months(df, winter_months)

  monthly_median <- df[, c(vars, "month", "lat", "lon")]
  monthly_median <- aggregate(. ~ month + lat + lon, monthly_median, median)

  df$monthly_median <- left_join(
    df[, c("month", "lat", "lon")],
    monthly_median,
    by = c("month" = "month", "lat" = "lat", "lon" = "lon")
  )[vars]
  df[vars] <- df[vars] - df$monthly_median

  df <- df[, c(vars, "time", "lat", "lon")]

  return(list(
    df = df,
    params = monthly_median
  ))
}


spring <- function(df, vars) {
  # extract spring months only
  spring_months <- c("March", "April", "May")
  df <- filter_months(df, spring_months)

  monthly_median <- df[, c(vars, "month", "lat", "lon")]
  monthly_median <- aggregate(. ~ month + lat + lon, monthly_median, median)

  df$monthly_median <- left_join(
    df[, c("month", "lat", "lon")],
    monthly_median,
    by = c("month" = "month", "lat" = "lat", "lon" = "lon")
  )[vars]
  df[vars] <- df[vars] - df$monthly_median

  df <- df[, c(vars, "time", "lat", "lon")]

  return(list(
    df = df,
    params = monthly_median
  ))
}


summer <- function(df, vars) {
  # extract summer months only
  summer_months <- c("June", "July", "August")
  df <- filter_months(df, summer_months)

  monthly_median <- df[, c(vars, "month", "lat", "lon")]
  monthly_median <- aggregate(. ~ month + lat + lon, monthly_median, median)

  df$monthly_median <- left_join(
    df[, c("month", "lat", "lon")],
    monthly_median,
    by = c("month" = "month", "lat" = "lat", "lon" = "lon")
  )[vars]
  df[vars] <- df[vars] - df$monthly_median

  df <- df[, c(vars, "time", "lat", "lon")]

  return(list(
    df = df,
    params = monthly_median
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

filter_months <- function(df, months_to_keep = c("October")) {
  # filter months to keep
  df$month <- months(df$time)
  df <- df[df$month %in% months_to_keep, ]
  return(df)
}