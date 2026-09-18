# %%
source("workflow/src/R/stats.R")
source("workflow/src/R/genpareto.R")
source("workflow/src/R/fitting.R")

set.seed(42)

margin <- data.frame(variable = rnorm(200))

distn <- list(
  cdf = cdf,
  threshold_selector = threshold_selector
)

result <- fit_margin_tails(margin, distn, two_tailed = TRUE, grid_i = "test_cell")

cat("scdf outside (0,1):", sum(result$scdf < 0 | result$scdf > 1, na.rm=TRUE), "\n")
cat("max scdf - ecdf diff:", max(abs(result$scdf - result$ecdf)), "\n")
plot(sort(result$ecdf), sort(result$scdf)); abline(0,1)

# %%