dcpi.raw <- read.csv("dcpi.csv", header=TRUE)
dcpi <- ts(dcpi.raw[,2], start=c(1947,2), frequency=12)
u.raw <- read.csv("unrate.csv", header=TRUE)
u <- ts(u.raw[,2], start=c(1948,1), frequency=12)
library(tstools)
dataset <- ts.intersect(u, dcpi)

fit <- tsreg(u, lags(dataset, 1))
fit

shocks.u <- fit$resids
plot(shocks.u)
plot(u)

irf <- tsreg(dcpi, lags(shocks.u, 0:12))
irf

plot(coefficients(irf)[-1], type="l", lwd=1.6,
     main="IRF of Inflation to 1 Point Unemployment Rate Rise")

fit2 <- tsreg(dcpi, ts.intersect(u, lags(dataset, 1)))
fit2

shocks.inf <- fit2$resids
plot(shocks.inf)
plot(dcpi)

irf2 <- tsreg(u, lags(shocks.inf, 0:12))
irf2

plot(coefficients(irf2)[-1], type="l", lwd=1.6,
     main="IRF of Unemployment to 1 Point Inflation Rate Rise")

irf.cumulative <- cumsum(coefficients(irf)[-1])
plot(irf.cumulative, type="l", lwd=1.6,
     main="IRF of Price Level to 1 Point Unemployment Rate Rise")

irf.cumulative <- cumsum(coefficients(irf2)[-1])
plot(irf.cumulative, type="l", lwd=1.6,
     main="IRF of Unemployment to 1 Point Inflation Rate Rise")
