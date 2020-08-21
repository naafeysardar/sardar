library(tstools)

data.raw <- read.csv("u-inf.csv", header = TRUE)
dataset <- ts(data.raw, start=c(1948, 1), frequency = 12)

colnames(dataset) <- c("u", "inf")
u <- dataset[,"u"]
inf <- dataset[,"inf"]

# Estimate RFVAR 

fit.u <- tsreg(u, lags(dataset, 1))

fit.inf <- tsreg(inf, lags(dataset, 1))

resu <- fit.u$resids
resinf <- fit.inf$resids

# Take residuals and split them into regimes

# Regime 1 (Pre-Greenspan)

resu.regime1 <- window(resu, end=c(1987, 8))
resinf.regime1 <- window(resinf, end=c(1987, 8))

# Variance

var.u1 <- var(resu.regime1)
var.inf1 <- var(resinf.regime1)

# Covariance

cov.1 <- cov(resu.regime1, resinf.regime1)

# Regime 2 (Greenspan Era)

resu.regime2 <- window(resu, start=c(1987, 9))
resinf.regime2 <- window(resinf, start=c(1987, 9))

# Variance

var.u2 <- var(resu.regime2)
var.inf2 <- var(resinf.regime2)

# Covariance

cov.2 <- cov(resu.regime2, resinf.regime2)


# Nonlinear equation solver (Just Identified Model)

objfnc.nleq <- function(par){
  c <- par[1]
  su <- par[2]
  sinf <- par[3]
  dev1 <- 0.04350375 - su
  dev2 <- 0.07218054 - c^2*su - sinf
  dev3 <-  -0.003451332 - c*su
  return(c(dev1, dev2, dev3))
}

library(nleqslv)
nleqslv(c(0.0, 0.1, 0.1), objfnc.nleq) 
