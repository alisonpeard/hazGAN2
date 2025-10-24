# Very basic XIMIS from Cook (2023)
set.seed(123)

# assume these are excesses above some threshold
R <- 250
x <- rweibull(n = R, scale = 1, shape = 2)

# fit weibull to get weibull
par <-c(0.5, 0.5)
nll <- function(x, par) {
  shape = par[1]
  scale = par[2]
  -sum(dweibull(x, shape = shape, scale = scale, log = TRUE))
}
fit <- optim(par, nll, x = x, method = "L-BFGS-B", lower = c(0.01, 0.01))

# use w to transform
omega <- fit$par[1]
z <- x^omega

# Gringorten plotting positions
gamma <- 0.57721566490153286

z <- sort(z, decreasing = TRUE)
y <- vector(length=R)
sig2 <- vector(length=R);

y[1] <- gamma + log(R)
sig2[1] <- pi^2 / 6

for (m in seq(from = 2, to = R)) {
  y[m] <- y[m-1] - 1 / m
  sig2[m] <- sig2[m-1] - 1 / (m^2)
}

# fit FT1 to y
w <- 1 / sig2

residuals <- function(par, z, y, w) {
  u <- par[1]
  d <- par[2]
  z_hat <- y * d + u
  sum(w * (z - z_hat)^2)
}

fit2 <- optim(c(1,1), residuals, z=z, y=y, w=w)
u <- fit2$par[1]
d <- fit2$par[2]

plot(y, (z - u) / d, main = "Gumbel probability plot",
     xlab = "(Z-U)/D", ylab = "y")
abline(0, 1, col = "red")

# simulate new samples
u_sample <- runif(1000, min=0.0001, max=0.9999)
y_sample <- -log(-log(u_sample))
z_sample <- d * y_sample + u
x_sample <- z_sample^(1/omega)


hist(x)
hist(x_sample)
