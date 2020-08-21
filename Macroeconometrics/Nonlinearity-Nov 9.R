library(tstools)
data.raw <- read.csv("u-inf.csv", header=TRUE)
u <- ts(data.raw[,1], start=c(1948,1), frequency=12)

fit.linear <- tsreg(u, lags(u, 1))

d <- lags(u, 1) < 6

rhs.nonlinear <- ts.combine(lags(u,1), d, d*lags(u,1))
fit.nonlinear <- tsreg(u, rhs.nonlinear)

udiff <- u - lags(u,1)
rhs.linear <- lags(udiff,1)
fit.linear <- tsreg(udiff, rhs.linear)
fit.linear

rhs.nonlinear <- ts.combine(lags(udiff,1), d, d*lags(udiff,1))
fit.nonlinear <- tsreg(udiff, rhs.nonlinear)


# Threshold Forecasting

thrpred <- function(u, du) {
  if (u < 6) {
    return(0.01 + 0.15*du)
  } else {
    return(-0.015+0.10*du)
  }
}

# Previous U = 5, Change in U = 0.4

thrpred(5, 0.4)

# Previous U = 7, Change in U = 0.4

thrpred(7, 0.4)


# Make Forecasts for 3 months into the future

# Create a vector with 3 elements, each one is NA

predictions <- rep(NA,3)

unem <- 5
dunem <- 0.4

# Make the first prediction and save it

dunem.hat <- thrpred(unem, dunem)
predictions[1] <- dunem.hat

# Second prediction

res <- fit.nonlinear$resids

# Update the unemployment regime

unem2 <- unem + dunem.hat + sample(res, size=1)
dunem.hat2 <- thrpred(unem2, dunem.hat)
predictions[2] <- dunem.hat2

# Third prediction

unem3 <- unem2 + dunem.hat2 + sample(res, size=1)
dunem.hat3 <- thrpred(unem3, dunem.hat2)
predictions[3] <- dunem.hat3
predictions

# Put this in a function

makepred <- function(unem, dunem) {
  predictions <- rep(NA,3)
  dunem.hat <- thrpred(unem, dunem)
  predictions[1] <- dunem.hat
  
  unem2 <- unem + dunem.hat + sample(res, size=1)
  dunem.hat2 <- thrpred(unem2, dunem.hat)
  predictions[2] <- dunem.hat2
  
  unem3 <- unem2 + dunem.hat2 + sample(res, size=1)
  dunem.hat3 <- thrpred(unem3, dunem.hat2)
  predictions[3] <- dunem.hat3
  
  return(predictions)
}

makepred(5, 0.4)

# Earlier we were randomly picking from residuals and treating
# that as a shock.

# Pass Shocks 

makepred <- function(unem, dunem, shocks) {
  dunem.hat <- thrpred(unem, dunem)
  predictions <- dunem.hat
  
  for (shock in shocks) {
    unem <- unem + dunem.hat + shock
    dunem.hat <- thrpred(unem, dunem.hat)
    predictions <- c(predictions, dunem.hat)
  }
  
  return(predictions)
}

# Shocks are 1.2 and 0.8

makepred(5, 0.4, c(1.2, 0.8))

# Make 12 Step Ahead Forecasts

makepred(5, 0.4, sample(res, size=11))

# Averaging over 1000 replications 

set.seed(500)
forecasts <- replicate(1000, {
  makepred(5, 0.4, sample(res, size=11))
}, simplify = "array")

forecasts
class(forecasts)
dim(forecasts)
rowMeans(forecasts)

sum(rowMeans(forecasts))

# After one year the increase in unemployment is 0.19 percent. 

# GDP Model

rgdp.raw <- read.csv("realgdp.csv", header = TRUE)
rgdp <- ts(rgdp.raw[,2], start=c(1947, 1), frequency = 4)

library(tstools)
drgdp <- pctChange(rgdp)
fit.linear <- tsreg(drgdp, lags(drgdp, 1))
fit.linear

# Recession Dummy

rec <- (rgdp < lags(rgdp, 1))
rhs.nonlinear <- ts.combine(lags(drgdp, 1), rec, rec*lags(drgdp, 1))
fit.nonlinear <- tsreg(drgdp, rhs.nonlinear)
fit.nonlinear

thrpred <- function(growth) {
  if (growth < 0) {
    -0.008 + 0.08*growth 
  } else {
    0.008 + 0.21*growth
  }
}
