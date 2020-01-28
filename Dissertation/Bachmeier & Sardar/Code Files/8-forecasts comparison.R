library(vars)
library(tstools)

gas.shock <- 10

c.nonlinear <- (0.51-0.12*auerbach-0.27*lags(auerbach,1)+0.085*lags(auerbach,2)+0.50*gdp+0.056*lags(gdp,1)
                +0.049*lags(gdp,2)-0.117*lags(pce,1)-0.01*lags(pce,2)-0.04*gas.shock-0.018*lags(gas,1)-0.004*lags(gas,2)
                -0.02*(auerbach*gas.shock)+0.011*lags(inter,1)-0.017*lags(inter,2))

c.linear <- (0.51+0.50*gdp+0.056*lags(gdp,1)+0.049*lags(gdp,2)-0.117*lags(pce,1)-0.01*lags(pce,2)-0.04*gas.shock-
               0.018*lags(gas,1)-0.004*lags(gas,2))

plot(c.nonlinear, type = "l", col = "black", xlab = "Year", ylab = "% Change", ylim= c(-2, 2), main = "PCE"); abline(h=0)
points(c.linear , col = "black", type="l", lty=2)

d <- c.nonlinear - c.linear
plot(d, type='n', ylab="Forecast Difference", ylim= c(-2, 2))
nberShade()
lines(d)

par(mfrow=c(1,1))

irf.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  
  rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(macro, 1:2), lags(gas, 0:2), lags(inter, 0:2))
  fit <- tsreg(macro, rhs)
  macro.nonlinear <- (coefficients(fit)[1]+coefficients(fit)[2]*auerbach+coefficients(fit)[3]*lags(auerbach,1)
                      +coefficients(fit)[4]*lags(auerbach,2)+coefficients(fit)[5]*gdp+coefficients(fit)[6]*lags(gdp,1)
                      +coefficients(fit)[7]*lags(gdp,2)+coefficients(fit)[8]*lags(pce,1)+coefficients(fit)[9]*lags(pce,2)
                      +coefficients(fit)[10]*gas.shock+coefficients(fit)[11]*lags(gas,1)+coefficients(fit)[12]*lags(gas,2)
                      +coefficients(fit)[13]*(auerbach*gas.shock)+coefficients(fit)[14]*lags(inter,1)+coefficients(fit)[15]*lags(inter,2))
  
  macro.linear <- (coefficients(fit)[1]+coefficients(fit)[5]*gdp+coefficients(fit)[6]*lags(gdp,1)+coefficients(fit)[7]*lags(gdp,2)
                   +coefficients(fit)[8]*lags(pce,1)+coefficients(fit)[9]*lags(pce,2)+coefficients(fit)[10]*gas.shock
                   +coefficients(fit)[11]*lags(gas,1)+coefficients(fit)[12]*lags(gas,2))
  
  plot(macro.nonlinear, type = "l", col = "blue", xlab = "Year", ylim = c(-2, 2), ylab = "% Change", main = varname); abline(h=0)
  points(macro.linear , col = "red", type="l", lty=2)
  
  legend("topright", c("Nonlinear Forecasts","Linear Forecasts"), fill=c("blue","red"))
  
  d <- macro.nonlinear - macro.linear
  plot(d, type = "n", col = "black", xlab = "Year", ylim = c(-2, 2), ylab = "% Change", main = varname)
  nberShade()
  abline(h=0)
  lines(d)
  return(macro.linear)
}

varinfo <- list(
  list(pce, "PCE"),
  list(nondurable, "Nondurables"),
  list(services, "Services"))
irfs <- lapply(varinfo, irf.linear)

irf.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  
  rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(macro, 1:2), lags(gas, 0:2), lags(inter, 0:2))
  fit <- tsreg(macro, rhs)
  macro.nonlinear <- (coefficients(fit)[1]+coefficients(fit)[2]*auerbach+coefficients(fit)[3]*lags(auerbach,1)
                      +coefficients(fit)[4]*lags(auerbach,2)+coefficients(fit)[5]*gdp+coefficients(fit)[6]*lags(gdp,1)
                      +coefficients(fit)[7]*lags(gdp,2)+coefficients(fit)[8]*lags(pce,1)+coefficients(fit)[9]*lags(pce,2)
                      +coefficients(fit)[10]*gas.shock+coefficients(fit)[11]*lags(gas,1)+coefficients(fit)[12]*lags(gas,2)
                      +coefficients(fit)[13]*(auerbach*gas.shock)+coefficients(fit)[14]*lags(inter,1)+coefficients(fit)[15]*lags(inter,2))
  
  macro.linear <- (coefficients(fit)[1]+coefficients(fit)[5]*gdp+coefficients(fit)[6]*lags(gdp,1)+coefficients(fit)[7]*lags(gdp,2)
                   +coefficients(fit)[8]*lags(pce,1)+coefficients(fit)[9]*lags(pce,2)+coefficients(fit)[10]*gas.shock
                   +coefficients(fit)[11]*lags(gas,1)+coefficients(fit)[12]*lags(gas,2))
  
  plot(macro.nonlinear, type = "l", col = "blue", xlab = "Year", ylim = c(-10, 10), ylab = "% Change", main = varname); abline(h=0)
  points(macro.linear , col = "red", type="l", lty=2)
  
  legend("topright", c("Nonlinear Forecasts","Linear Forecasts"), fill=c("blue","red"))
  
  d <- macro.nonlinear - macro.linear
  plot(d, type = "l", col = "black", xlab = "Year", ylim = c(-10, 10), ylab = "% Change", main = varname)
  nberShade()
  abline(h=0)
  lines(d)
}

varinfo <- list(
  list(durable, "Durables"))
irfs <- lapply(varinfo, irf.linear)

