library(withr)
library(tidync)

library(eva) # Bader method
library(threshr) # Northrop, Attalides, and Jonathon (2017) method
# library(extRemes) 
# library(ismev)
# library(evd)
library(mev)


# Threshold stability plots of Davison and Smith (1990)
source("/Users/alison/Local/github/hazGAN2/workflow/src/R/genpareto.R")
 
INPUT <- "../results/testing/pois.nc"
I_POI <- 1
I_VAR <- 1

src <- tidync(INPUT)

df  <- src |> hyper_tibble(force = TRUE)
pois <- unique(df$poi)
vars <- unique(df$field)
poi <- pois[I_POI]
var <- vars[I_VAR]

df <- df[(df$poi == poi) & (df$field == var), ]
print(paste0("Extracted df with length: ", dim(df)[1], " for (", poi, ", ", var, ")"))
x <- df$anomaly
vx <- quantile(x, 0.7)
hist(x[x>vx])

thresholds <- quantile(x, probs = seq(0.7, 0.98, length.out = 28))

"using mev"
mev::tstab.gpd(x, thresh=thresholds)


"https://cran.r-project.org/web/packages/threshr/vignettes/threshr-vignette.html#single-training-threshold"
plot(threshr::stability(x, u_vec = thresholds))
x_cv <- threshr::ithresh(x, u_vec = thresholds, n = 100)
plot(x_cv)
summary(x_cv)
print(paste0("threshr suggests quantile ", summary(x_cv)[4], " and value ", summary(x_cv)[3], "."))
plot(x_cv, which_u = "best")
plot(predict(x_cv, n_years = c(100, 1000), type = "d", npy=22))
plot(predict(x_cv, type = "d", npy=22, which_u="all"))

"my code"
fit <- threshold_selector(x, 1)
fit

"Bader method"
alpha = 0.05
fits <- gpdSeqTests(
    x, thresholds = thresholds, method = "ad", nsim = 100
  )

valid_pk <- fits$ForwardStop
k <- min(which(valid_pk > alpha)) # lowest valid threshold
k <- max(which(valid_pk <= alpha))
print(paste0("Bader method selected k=", k))
p.value <- fits$p.values[k]
pk <- fits$ForwardStop[k]
print(paste0("Bader method got p=", p.value, " and pk=", pk))
params   = list(
  thresh = fits$threshold[k],
  scale  = fits$est.scale[k],
  shape  = fits$est.shape[k]
)
print(paste0("Fitted params loc=", params$thresh, ", scale=", params$scale, ", shape=", params$shape))


"Murphy (2025)"

setwd("/Users/alison/Local/github/automated_threshold_selection/")
source("src/eqd.R")
set.seed(12345)
data_test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
thresholds1 <- seq(0.5, 2.5, by=0.1)
example1 <- eqd(data_test1, thresh = thresholds1, k=100, m=500)
example1

plot(thresholds1, example1$dists, xlab="Threshold", ylab="Metric value")
abline(v=example1$thresh, col="red", lwd=2)


example1 <- eqd(x, thresh = thresholds, k=100, m=500)
example1

plot(thresholds, example1$dists, xlab="Threshold", ylab="Metric value")
abline(v=example1$thresh, col="red", lwd=2)
