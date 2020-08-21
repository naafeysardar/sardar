library(AER)

data("USMacroG")
gdp <- USMacroG[,"gdp"]

plot(gdp)

library(tstools)

fit <- tsreg(gdp, lags(gdp, 8:11))

tr <- fit$fitted
cycle <- fit$resids

plot(tr, main="Real GDP Trend")
plot(cycle, main="Real GDP Business Cycle")

plot(ts.combine(gdp, tr), plot.type = "single",
     lty=c(1,2))