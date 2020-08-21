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


objfnc.efficient <- function(par, data) {
  a0 <- par[1]
  a1 <- par[2]
  a2 <- par[3]
  b0 <- par[4]
  b1 <- par[5]
  b2 <- par[6]
  c <- par[7]
  su <- par[8]
  
  u <- data[,"u"]
  inf <- data[,"inf"]
  
  res.u <- u - a0 - a1*lags(u, 1) - a2*lags(inf, 1)
  res.inf <- inf - b0 - b1*lags(u, 1) - b2*lags(inf, 1)
  
  # Moment Conditions for OLS
  
  dev1 <- res.u
  dev2 <- res.u*lags(u, 1)
  dev3 <- res.u*lags(inf, 1)
  dev4 <- res.inf
  dev5 <- res.inf*lags(u, 1)
  dev6 <- res.inf*lags(inf, 1)
  
  # Moment Conditions for SVAR
  
  dev7 <- res.u^2 - su
  dev8 <- res.inf^2 - (c^2+1)*su
  dev9 <- res.u*res.inf - c*su
  
  return(cbind(dev1, dev2, dev3, dev4, dev5, dev6, dev7, dev8, dev9))
  
}

objfnc.efficient(c(coefficients(fit.u), coefficients(fit.inf), 
                 0, 0.1), dataset)


# Default Iterations are 500

gmm(objfnc.efficient, dataset, t0=c(coefficients(fit.u), coefficients(fit.inf), 
                                    0, 0.1))

# 10,000 Iterations

gmm(objfnc.efficient, dataset, t0=c(coefficients(fit.u), coefficients(fit.inf), 
                                    0, 0.1), control=list(maxit=10000))


# Recover the shocks

fit.u <- tsreg(u, lags(dataset, 1))
fit.inf <- tsreg(inf, lags(dataset, 1))

res.u <- fit.u$resids
res.inf <- fit.inf$resids

res.u[1]
res.inf[1]

A <- matrix(c(1, -0.08, 0, 1), ncol=2)
res <- cbind(res.u, res.inf)

structural.shocks <- solve(A, t(res))
structural.shocks

eu <- structural.shocks[1,]
einf <- structural.shocks[2,]

eu <- ts(eu, start = start(res), frequency = 12)

einf <- ts(einf, start = start(res), frequency = 12)

plot(eu)
sum(window(eu, start=c(2008, 1), end=c(2009,12)))

plot(einf)