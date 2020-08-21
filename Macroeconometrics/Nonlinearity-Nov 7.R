library(tstools)
data.raw <- read.csv("u-inf.csv", header=TRUE)
u <- ts(data.raw[,1], start=c(1948,1), frequency=12)

fit.linear <- tsreg(u, lags(u, 1))

d <- lags(u, 1) < 6

rhs.nonlinear <- ts.combine(lags(u,1), d, d*lags(u,1))
fit.nonlinear <- tsreg(u, rhs.nonlinear)

udiff <- u - lags(u,1)
rhs.linear <- lags(udiff,1)
fit.linear <- tsreg(udiff, rhs.linear)
fit.linear

rhs.nonlinear <- ts.combine(lags(udiff,1), d, d*lags(udiff,1))
fit.nonlinear <- tsreg(udiff, rhs.nonlinear)


