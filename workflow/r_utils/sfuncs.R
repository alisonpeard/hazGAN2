# add custom named sfuncs here
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(methods) # for S4 methods / snakemake
})
# source("workflow/r_utils/utils.R")

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
  df <- df[, c(vars, "month", "grid")]
  monthly_median <- aggregate(. ~ month + grid, df, median)
  df$monthly_median <- left_join(
    df[, c("month", "grid")],
    monthly_median,
    by = c("month" = "month", "grid" = "grid")
  )[vars]
  df[vars] <- df[vars] - df$monthly_median
  return(list(
    df = df[vars],
    params = monthly_median
  ))
}

autumn <- function(df, vars) {
  # extract autumn months only
  autumn_months <- c("Sep", "Oct", "Nov")
  df <- filter_months(df, autumn_months)
  return(df)
}


winter <- function(df, vars) {
  # extract winter months only
  winter_months <- c("Dec", "Jan", "Feb")
  df <- filter_months(df, winter_months)
  return(df)
}

spring <- function(df, vars) {
  # extract spring months only
  spring_months <- c("Mar", "Apr", "May")
  df <- filter_months(df, spring_months)
  return(df)
}

summer <- function(df, vars) {
  # extract summer months only
  summer_months <- c("Jun", "Jul", "Aug")
  df <- filter_months(df, summer_months)
  return(df)
}


monthly_medians <- function(df, var) {
  df <- df[, c(var, "time", "grid")]
  df$month <- months(df$time)
  monthly_median <- aggregate(. ~ month + grid,
                              df[, c(var, "grid", "month")],
                              median)
  return(monthly_median)
}

filter_months <- function(df, months_to_keep = c("Oct")) {
  # filter months to keep
  df$month <- months(df$time)
  df <- df[df$month %in% months_to_keep, ]
  return(df)
}