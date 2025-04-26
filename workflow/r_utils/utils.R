# remember to select source on save
library(zoo)
library(dplyr)
library(lubridate)
library(data.table)
library(progress)  # Add this
library(goftest)


`%ni%` <- Negate(`%in%`)


progress_bar <- function(n, prefix = "", suffix = "") {
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  function(i) {
    utils::setTxtProgressBar(pb, i)
    if (i == n) close(pb)
  }
}

