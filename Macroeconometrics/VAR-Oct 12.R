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