library(tstools)

data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

colnames(dataset) <- c("u","inf")

u <- dataset[,"u"]

inf <- dataset[,"inf"]


# Estimate RFVAR

fit.u <- tsreg(u, lags(dataset, 1))
fit.u

fit.inf <- tsreg(inf, lags(dataset, 1))
fit.inf

# Find RFVAR residuals covariances

res1 <- fit.u$resids
res2 <- fit.inf$resids

s2.1 <- var(res1)
s2.2 <- var(res2)

s12 <- cov(res1, res2)

s2.1
s2.2
s12

# Numerical Optimization (Just Identified Model)

objfnc <- function(par){
  c <- par[1]
  su <- par[2]
  sinf <- par[3]
  dev1 <- 0.04350375 - su
  dev2 <- 0.07218054 - c^2*su - sinf
  dev3 <-  -0.003451332 - c*su
  return(dev1^2 + dev2^2 + dev3^2)
}

objfnc(c(0, 0.1, 0.1))

# Specifying starting values 

optim(c(0, 0.1, 0.1), objfnc)

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

# GMM 

# Just Identified Model

res <- ts.combine(res1, res2)

objfnc.gmm <- function(par, data){
  c <- par[1]
  su <- par[2]
  sinf <- par[3]
  dev1 <- data[,"res1"]^2 - su
  dev2 <- data[,"res2"]^2 - c^2*su - sinf
  dev3 <- data[,"res1"]*data[,"res2"] - c*su
  return(cbind(dev1, dev2, dev3))
}

library(gmm)

gmm(objfnc.gmm, res, t0=c(-0.1, 0.1, 0.1))

# Over Identified Model (Use GMM when over identified)

objfnc.oid <- function(par, data){
  c <- par[1]
  su <- par[2]
  dev1 <- data[,"res1"]^2 - su
  dev2 <- data[,"res2"]^2 - (c^2+1)*su
  dev3 <- data[,"res1"]*data[,"res2"] - c*su
  return(cbind(dev1, dev2, dev3))
}

gmm(objfnc.oid, res, t0=c(-0.1, 0.1))


