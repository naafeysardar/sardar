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

# Define Rising Regimes

rising <- which(doil > 0)
rising

# Take the values in which oil prices are not rising (!)

falling <- which(!(doil > 0))
falling

doil.plus <- doil[rising]
doil.minus <- doil[falling]

mean(doil.plus)

# If we are in the increasing oil price regime, the increase will be 9% in that quarter.

mean(doil.minus)

# If we are in the decreasing oil price regime, the decrease will be 6% in that quarter.

# Simulation Approach to computing IRFs

gdp.irf <- function(doil1) {
  dgdp1 <- sample(dgdp, size=1)
  oilshock <- sample(fit.oil$resids, size=1)
  gdpshock <- sample(fit.gdp$resids, size=1)
  return(gdp.fcsts(dgdp1, doil1, 1.0, gdpshock, oilshock) - gdp.fcsts(dgdp1, doil1, 0.0, gdpshock, oilshock))
}

# If we start out in the increasing regime, the impact of oil price shock would be that gdp decreases by 1.8%.

gdp.irf(0.09)

# If we start out in the decreasing regime, the impact of oil price shock would be that gdp decreases by 4.8%.

gdp.irf(-0.06)

set.seed(4000)

doil1 <- 0.9    # Increasing Regime
doil1 <- -0.6  

irfs <- replicate(1000, {
  dgdp1 <- sample(dgdp, size=1)
  oilshock <- sample(fit.oil$resids, size=1)
  gdpshock <- sample(fit.gdp$resids, size=1)
  return(gdp.fcsts(dgdp1, doil1, 1.0, gdpshock, oilshock) - gdp.fcsts(dgdp1, doil1, 0.0, gdpshock, oilshock))
}, simplify = "array")

# Average IRF

rowMeans(irfs)



# Local Projections

rhs <- ts.combine(lags(doil, 0:1), lags(doil^2, 0:1), lags(doil^3, 0:1), lags(dgdp, 1), lags(dgdp^2, 1), lags(dgdp^3, 1))

fit0 <- tsreg(gdp, rhs)