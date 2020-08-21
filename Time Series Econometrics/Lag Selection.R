dcpi.raw <- read.csv("dcpi.csv", header=TRUE)
dcpi <- ts(dcpi.raw[,2], start = c(1947,2), frequency = 12)
plot (dcpi)

unrate.raw <- read.csv("unrate.csv", header=TRUE)
unrate <- ts(unrate.raw[,2], start = c(1948,1), frequency = 12)
plot (unrate)

library(tstools)
fit <- tsreg(dcpi, unrate)
summary(fit)
plot(fit$resids)
plot(residuals(fit))

fit12 <- tsreg(dcpi, lags(dcpi, 1:12))
AIC(fit12)
BIC(fit12)

arima.output <- arima(dcpi, order = c(1, 0, 2))
arima.output
AIC(arima.output)
predict(arima.output, 1)
predict(arima.output, 12)
plot(predict(arima.output, 12)$pred)

criteria <- function(lags) {
  fit <- arima(unrate, order=c(lags[1], 0, lags[2]))
  list(sic=BIC(fit), ar.lags=lags(1), ma.lags=lags(2))
  }

potential.lags <- expand.grid(0:2, 0:2)
apply(potential.lags, MARGIN = 1, criteria)