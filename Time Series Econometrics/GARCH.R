library(urca)
data("Raotbl3")

lc <- Raotbl3[,"lc"]

lc <- ts(lc, start=c(1966,4), frequency=4)
dlc <- diff(lc)
plot(dlc)

## There is volatility based on the plot ##

library(fGarch)
fit <- garchFit(~ garch(1,1), data=dlc)
summary(fit)

fit@fitted

## Getting the fitted value of volatility (conditional variance).
fit@h.t
plot(fit@h.t)

## Forecasting the dependent variable
predict(fit, n.ahead=4)


## Swiss Pension func Index -
x = as.timeSeries(data(LPP2005REC))

## garchFit
fit = garchFit(LPP40 ~ garch(1, 1), data = 100*x, trace = FALSE)
fit

## volatility - 
# Standard Deviation:
volatility = volatility(fit, type = "sigma")
head(volatility)
class(volatility)
# Variance:
volatility = volatility(fit, type = "h")
head(volatility)
class(volatility)

## slot - 
volatility = slot(fit, "sigma.t")
head(volatility)
class(volatility)
volatility = slot(fit, "h.t")
head(volatility)
class(volatility)