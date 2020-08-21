library(tstools)
library(AER)

data(USMacroG)

colnames(USMacroG)

dgdp <- 400*pctChange(USMacroG[,"gdp"])

plot(dgdp)

mean(dgdp)

abline(a=mean(dgdp), b=0)

tsreg(dgdp, lags(dgdp,1))

arma11 <- arima(dgdp, order = c(1,0,1))
arma11

arma12 <- arima(dgdp, order = c(1,0,2))
arma12

arma22 <- arima(dgdp, order = c(2,0,2))
arma22

AIC(arma11)
AIC(arma12)
AIC(arma22)


BIC(arma11)
BIC(arma12)
BIC(arma22)

criteria <- function(pq) {
  fit <- arima(dgdp, order=c(pq[1], 0, pq[2]))
  list(aic=AIC(fit), sic=BIC(fit), ar.lags=pq[1],
       ma.lags=pq[2])  
}

criteria(c(1,1))

potential.lags <- expand.grid(0:4, 0:4)
potential.lags

infcrit <- apply(potential.lags, MARGIN = 1, criteria)
infcrit[[1]]


aic.value <- function(crit){
  return(crit$aic)
}

aic.value((infcrit[[2]]))

aic.values <- lapply(infcrit, aic.value)
aic.values

which.min(aic.values)
infcrit[[2]]

