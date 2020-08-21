library(vars)

data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]

inf <- dataset[,"inf"]


varfit <- VAR(dataset, p=1)
varfit

# Make Forecasts

varpred <- predict(varfit, n.ahead = 12)
varpred

# Inflation Forecast

varpred$fcst$"inf"

varpred$fcst$"inf"[, "fcst"]

library(tstools)
getVarForecasts(varfit, var="inf", n=1:12, start = c(2018, 8))

plot(getVarForecasts(varfit, var="inf", n=1:12, start = c(2018, 8)))

# Inflation forecast pushed back to the mean. 

plot(getVarForecasts(varfit, var="u", n=13:36, start = c(2019, 8)))

# U is non-stationery which is why it doesn't move back to the mean. 



varfit2 <- VAR(dataset, p=2)
getVarForecasts(varfit2, var="inf", n=1:12, start = c(2018, 8))


# Test for no. of lags

VARselect(dataset, lag.max=12)
varfit <- VAR(dataset, ic="SC", lag.max=12)

# 2-Step Ahead Forecast

inf2 <- tsreg(inf, lags(ts.combine(inf, u), 2))
inf2

last(inf)
last(u)

0.1239+0.4460*0.2+0.0052*3.9
0.23338*12


