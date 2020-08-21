# Upload Data

data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]

inf <- dataset[,"inf"]

# Estimate model using our sample

fit <- tsreg(inf, u)

res <- residuals(fit)

# Fixed Design Wild Bootstrap

# Standard Normal Pick Distribution

set.seed(200)
sim.beta <- replicate(100, {
  inf.sim <- fit$fitted + rnorm(length(u))*fit$resids
  fit.sim <- tsreg(inf.sim, u)
  coefficients(fit.sim)[2]
})

sd(sim.beta)

# Rademacher Distribution

set.seed(200)
sim.beta <- replicate(100, {
  pick <- sample(c(-1, 1), size = length(u), replace = TRUE)
  inf.sim <- fit$fitted + pick*fit$resids
  fit.sim <- tsreg(inf.sim, u)
  coefficients(fit.sim)[2]
})

sd(sim.beta)


# Non-Pivotal

# Choosing value of beta = 0.2. Want to test that
# the SE(beta) is not sensitive to the chosen 
# values of alpha and beta.

set.seed(200)
simoutput <- replicate(1000, {
  inf.sim <- 0.2*u+rnorm(length(u))
  coefficients(tsreg(inf.sim,u))[2]
})

sd(simoutput)

set.seed(200)
simoutput <- replicate(1000, {
  inf.sim <- 2.5+0.0*u+rnorm(length(u))
  coefficients(tsreg(inf.sim,u))[2]
})

sd(simoutput)


# VAR

uinf <- ts.combine(inf, u)
inf.fit <- tsreg(inf, lags(uinf, 1))
u.fit <- tsreg(u, lags(uinf, 1))

inf.fit
u.fit

tsp(u)    # Data starts from Jan 1948 and goes up till July 2018.

last(inf)
last(u)

# August 2018 Inflation Forecast
0.09 + 0.59*0.2 + 0.005*3.9
0.2275

# Annualizing
0.2275*12

# September 2018 Inflation Forecast
# First find u for August 2018 and 
# then plug in the inflation equation.


# Multi Step Forecasting

Z <- c(0.2, 3.9, 0.09, 0.04)
F <- matrix(c(0.59, 0.03, 0, 0,
            0.005, 0.99, 0, 0,
            1, 0, 1, 0,
            0, 1, 0, 1), ncol=4)
F

# One Step Ahead Forecast 

F %*% Z

# Two Step Ahead Forecast 

F %*% F %*% Z

# Three Step Ahead Forecast 

F %*% F %*% F %*% Z


