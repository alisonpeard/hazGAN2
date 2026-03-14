# %%
source("workflow/src/R/stats.R")
source("workflow/src/R/genpareto.R")

set.seed(42)
x_train <- rnorm(200)

params <- list(
  thresh_upper = 1.2,
  scale_upper  = 0.5,
  shape_upper  = 0.1,
  thresh_lower = -1.2,
  scale_lower  = 0.5,
  shape_lower  = 0.1
)

f <- scdf_wb(x_train, params, cdf = cdf)

test_vals <- c(2.0, 1.5, -2.0, -1.5, 1.2, -1.2)
cat("Results:\n")
print(f(test_vals))
# %%