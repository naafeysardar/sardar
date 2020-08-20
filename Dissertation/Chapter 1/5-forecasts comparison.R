library(vars)
library(tstools)
library(tis, exclude = c("lags", "month", "quarter", "year"))

gas.shock <- 10

par(mfrow=c(2,1))

irf.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  
  rhs.nl <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(macro, 1:2), lags(gas, 0:2), lags(inter, 0:2))
  fit.nl <- tsreg(macro, rhs.nl)
  macro.nonlinear <- (coefficients(fit.nl)[1]+coefficients(fit.nl)[2]*auerbach+coefficients(fit.nl)[3]*lags(auerbach,1)
                      +coefficients(fit.nl)[4]*lags(auerbach,2)+coefficients(fit.nl)[5]*gdp+coefficients(fit.nl)[6]*lags(gdp,1)
                      +coefficients(fit.nl)[7]*lags(gdp,2)+coefficients(fit.nl)[8]*lags(pce,1)+coefficients(fit.nl)[9]*lags(pce,2)
                      +coefficients(fit.nl)[10]*gas.shock+coefficients(fit.nl)[11]*lags(gas,1)+coefficients(fit.nl)[12]*lags(gas,2)
                      +coefficients(fit.nl)[13]*(auerbach*gas.shock)+coefficients(fit.nl)[14]*lags(inter,1)+coefficients(fit.nl)[15]*lags(inter,2))
  
  rhs.l <- ts.intersect(lags(gdp, 0:2), lags(macro, 1:2), lags(gas, 0:2))
  fit.l <- tsreg(macro, rhs.l)
  macro.linear <- (coefficients(fit.l)[1]+coefficients(fit.l)[2]*gdp+coefficients(fit.l)[3]*lags(gdp,1)+coefficients(fit.l)[4]*lags(gdp,2)
                   +coefficients(fit.l)[5]*lags(pce,1)+coefficients(fit.l)[6]*lags(pce,2)+coefficients(fit.l)[7]*gas.shock
                   +coefficients(fit.l)[8]*lags(gas,1)+coefficients(fit.l)[9]*lags(gas,2))
  
  plot(macro.nonlinear, type = "l", col = "blue", xlab = "Year", ylim = c(-2, 2), ylab = "% Change", main = varname); abline(h=0)
  points(macro.linear , col = "red", type="l", lty=2)
  
  legend("topright", c("Nonlinear Forecasts","Linear Forecasts"), fill=c("blue","red"))
  
  d <- macro.nonlinear - macro.linear
  plot(d, type = "n", col = "black", xlab = "Year", ylim = c(-1, 1), ylab = "% Change", main = varname)
  nberShade()
  abline(h=0)
  lines(d)
  
  return(list(macro.nonlinear, macro.linear))
}

varinfo <- list(
  list(pce, "PCE"),
  list(nondurable, "Nondurables"),
  list(services, "Services"))
irfs <- lapply(varinfo, irf.linear)

irf.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  
  rhs.nl <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(macro, 1:2), lags(gas, 0:2), lags(inter, 0:2))
  fit.nl <- tsreg(macro, rhs.nl)
  macro.nonlinear <- (coefficients(fit.nl)[1]+coefficients(fit.nl)[2]*auerbach+coefficients(fit.nl)[3]*lags(auerbach,1)
                      +coefficients(fit.nl)[4]*lags(auerbach,2)+coefficients(fit.nl)[5]*gdp+coefficients(fit.nl)[6]*lags(gdp,1)
                      +coefficients(fit.nl)[7]*lags(gdp,2)+coefficients(fit.nl)[8]*lags(pce,1)+coefficients(fit.nl)[9]*lags(pce,2)
                      +coefficients(fit.nl)[10]*gas.shock+coefficients(fit.nl)[11]*lags(gas,1)+coefficients(fit.nl)[12]*lags(gas,2)
                      +coefficients(fit.nl)[13]*(auerbach*gas.shock)+coefficients(fit.nl)[14]*lags(inter,1)+coefficients(fit.nl)[15]*lags(inter,2))
  
  rhs.l <- ts.intersect(lags(gdp, 0:2), lags(macro, 1:2), lags(gas, 0:2))
  fit.l <- tsreg(macro, rhs.l)
  macro.linear <- (coefficients(fit.l)[1]+coefficients(fit.l)[2]*gdp+coefficients(fit.l)[3]*lags(gdp,1)+coefficients(fit.l)[4]*lags(gdp,2)
                   +coefficients(fit.l)[5]*lags(pce,1)+coefficients(fit.l)[6]*lags(pce,2)+coefficients(fit.l)[7]*gas.shock
                   +coefficients(fit.l)[8]*lags(gas,1)+coefficients(fit.l)[9]*lags(gas,2))
  
  plot(macro.nonlinear, type = "l", col = "blue", xlab = "Year", ylim = c(-10, 10), ylab = "% Change", main = varname); abline(h=0)
  points(macro.linear , col = "red", type="l", lty=2)
  
  legend("topright", c("Nonlinear Forecasts","Linear Forecasts"), fill=c("blue","red"))
  
  d <- macro.nonlinear - macro.linear
  plot(d, type = "l", col = "black", xlab = "Year", ylim = c(-3, 3), ylab = "% Change", main = varname)
  nberShade()
  abline(h=0)
  lines(d)
}

varinfo <- list(
  list(durable, "Durables"))
irfs <- lapply(varinfo, irf.linear)

