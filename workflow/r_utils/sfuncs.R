# add custom named sfuncs here
library(dplyr)
library(lubridate)
library(methods) # for S4 methods / snakemake
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
  return(df[vars])
}

monthly_medians <- function(df, var) {
  df <- df[, c(var, "time", "grid")]
  df$month <- months(df$time)
  monthly_median <- aggregate(. ~ month + grid,
                              df[, c(var, "grid", "month")],
                              median)
  return(monthly_median)
}