# gamma = 1, d = 0

theta <- function(y){
  den <- 1 + exp(-1*y)
  return(1/den)
  }

curve(theta, from =-5, to = 5)

# gamma = -5, d = 0

theta <- function(y){
  den <- 1 + exp(-5*y)
  return(1/den)
}

curve(theta, from = -5, to = 5)

# gamma = 0.05, d = 0

theta <- function(y){
  den <- 1 + exp(-0.05*y)
  return(1/den)
}

curve(theta, from = -5, to = 5)




gdp <- readRDS("gdp.RDS")
oil <- readRDS("oilppi.RDS")

library(tstools)
dgdp <- pctChange(gdp)
doil.raw <- pctChange(oil)
doil <- toQuarterly(doil.raw, sum)

# Define threshold

d <- (doil < 0)
inter <- lags(d, 1) * lags(doil, 1)

fit.gdp <- tsreg(dgdp, cbind(lags(dgdp, 1), lags(doil, 0:1), lags(d, 1), lags(inter, 0:1)))

fit.oil <- tsreg(doil, cbind(lags(dgdp, 1), lags(doil, 1)))

pred.oil <- function(dgdp1, doil1) {
  return(0.002 + 0.466*dgdp1 + 0.064*doil1)
}

pred.gdp <- function(dgdp1, doil0, doil1) {
  return(0.006 + 0.33*dgdp1 + 0.003*doil0 - 0.019*doil1 - 0.0005*(doil1 < 0)
         - (doil0 < 0)*doil0*0.03 - (doil1 < 0)*doil1*(-0.002))
}

first(dgdp)
first(doil)
first(fit.oil$resids)
first(fit.gdp$resids)

# Oil Shock = 1.0

`doil[T+1]` <- pred.oil(-0.00113966, 0.0875) + 1.0
`dgdp[T+1]` <- pred.gdp(-0.00113966, `doil[T+1]`, 0.0875) - 0.006312693

`doil[T+2]` <- pred.oil(`dgdp[T+1]`, `doil[T+1]`) + 0.01882849
`dgdp[T+2]` <- pred.gdp(`dgdp[T+1]`, `doil[T+2]`, `doil[T+1]`)

c(`dgdp[T+1]`, `dgdp[T+2]`)

# Write a function

gdp.fcsts <- function(dgdp1, doil1, oilshock.1,
                      gdpshock.1, oilshock.2) {
  `doil[T+1]` <- pred.oil(dgdp1, doil1) + oilshock.1
  `dgdp[T+1]` <- pred.gdp(dgdp1, `doil[T+1]`, doil1) + gdpshock.1
  
  `doil[T+2]` <- pred.oil(`dgdp[T+1]`, `doil[T+1]`) + oilshock.2
  `dgdp[T+2]` <- pred.gdp(`dgdp[T+1]`, `doil[T+2]`,
                          `doil[T+1]`)
  
  return(c(`dgdp[T+1]`, `dgdp[T+2]`))
  
}

gdp.fcsts(-0.00113966, 0.0875, 1.0, -0.006312693, 0.01882849)
gdp.fcsts(-0.00113966, 0.0875, 0.0, -0.006312693, 0.01882849)

x <- gdp.fcsts(-0.00113966, 0.0875, 1.0, -0.006312693, 0.01882849)- gdp.fcsts(-0.00113966, 0.0875, 0.0, -0.006312693, 0.01882849)

