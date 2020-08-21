library(tstools)

data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]

inf <- dataset[,"inf"]


# Inflation Equation

fit.inf <- tsreg(inf, lags(dataset, 1))
fit.inf

fit.u <- tsreg(u, ts.combine(inf, lags(dataset, 1)))
fit.u

# 12-step ahead forecast

fit.u12 <- tsreg(u, lags(dataset, 12))
fit.u12

library(vars)

# inf has contemporaneous effect on u

ds <- ts.combine(inf, u)

varfit <- VAR(ds, ic="SC", lag.max = 6)
varfit

# Plot IRF for u following supply shock

plot(irf(varfit, impulse = "inf", response= "u", boot=FALSE, n.ahead=24))

# u has contemporaneous effect on inf

ds2 <- ts.combine(u, inf)

varfit2 <- VAR(ds2, ic="SC", lag.max=6)
plot(irf(varfit2, impulse = "inf", response = "u", boot=FALSE, n.ahead=24))