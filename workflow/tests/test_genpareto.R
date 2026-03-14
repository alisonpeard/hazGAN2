source("workflow/src/R/genpareto.R")
# run diagnostic tests
library(eva)
print("My case")
x <- rnorm(1249)
threshold_selector(x, "1249 points");

print("too few points case:")
x <- rnorm(10)
threshold_selector(x, "too-few points");

print("zero-variance tail case:")
x <- c(rnorm(100), rep(5, 50))
threshold_selector(x, "zero-variance tail");

print("extremely heavy (Cauchy) tail (ξ<-0.5)")
x <- abs(rcauchy(200))
threshold_selector(x, "heavy tail")