library(tstools)
data.raw <- read.csv("u-inf.csv", header=TRUE)
dataset <- ts(data.raw, start=c(1948,1), frequency=12)
colnames(dataset) <- c("u", "inf")
u <- dataset[,"u"]
inf <- dataset[,"inf"]
fit.u <- tsreg(u, lags(dataset, 1))
fit.inf <- tsreg(inf, lags(dataset, 1))

res1 <- fit.u$resids
res2 <- fit.inf$resids
s2.1 <- var(res1)
s2.2 <- var(res2)
s12 <- cov(res1, res2)

# Moment Conditions
# Assumption b=c

objfnc.nleq <- function(par) {
  c <- par[1]
  su <- par[2]
  sinf <- par[3]
  dev1 <- 0.04350375 - su - c^2*sinf
  dev2 <- 0.07218054 - c^2*su - sinf
  dev3 <- -0.003451332 - c*su - c*sinf 
  return(c(dev1, dev2, dev3))
}

library(nleqslv)
sol <- nleqslv(c(0, 0.1, 0.1), objfnc.nleq)
b <- sol$x[1]
b
irf0 <- matrix(c(1.0, b))
irf0
library(vars)
varfit <- VAR(dataset, p=1)
phimat <- Phi(varfit, nstep=12)
irf1 <- phimat[,,2] %*% matrix(irf0)
irf1
irf2 <- phimat[,,3] %*% matrix(irf0)
irf2
irf.calc <- function(h) {
  return(phimat[,,h+1] %*% matrix(irf0))
}
irfs <- lapply(1:12, irf.calc)
irfs

get.inf <- function(irf) {
  return(irf[2])
}
response.inf <- sapply(irfs, get.inf)
response.inf <- c(b, response.inf)
plot(ts(response.inf), main="Response of Inflation to Demand Shock")

set.seed(100)
irfboot <- replicate(100, {
  z <- sample(c(-1,1), size=length(res1), replace=TRUE)
  u.sim <- fit.u$fitted + z*res1
  inf.sim <- fit.inf$fitted + z*res2
  fit.usim <- tsreg(u.sim, lags(dataset, 1))
  fit.infsim <- tsreg(inf.sim, lags(dataset, 1))
  resusim <- fit.usim$resids
  resinfsim <- fit.infsim$resids
  
  s2.1 <- var(resusim)
  s2.2 <- var(resinfsim)
  s12 <- cov(resusim, resinfsim)
  objfnc.nleq <- function(par) {
    c <- par[1]
    su <- par[2]
    sinf <- par[3]
    dev1 <- s2.1 - su - c^2*sinf
    dev2 <- s2.2 - c^2*su - sinf
    dev3 <- s12 - c*su - c*sinf 
    return(c(dev1, dev2, dev3))
  }
  
  library(nleqslv)
  solsim <- nleqslv(c(0, 0.1, 0.1), objfnc.nleq)
  bsim <- solsim$x[1]
  
  irf0 <- matrix(c(1.0, bsim))
  
  dataset.sim <- ts.combine(u.sim, inf.sim)
  varfit <- VAR(dataset.sim, p=1)
  phimat <- Phi(varfit, nstep=12)
  
  irf.calc <- function(h) {
    return(phimat[,,h+1] %*% matrix(irf0))
  }
  irfs <- lapply(1:12, irf.calc)
  response.inf <- sapply(irfs, get.inf)
  c(bsim, response.inf)
}, simplify="array")

class(irfboot)
dim(irfboot)

irflower <- function(h) {
  return(quantile(irfboot[h+1,], probs=0.025))
}
lower <- sapply(0:12, irflower)
lower

irfupper <- function(h) {
  return(quantile(irfboot[h+1,], probs=0.975))
}
upper <- sapply(0:12, irfupper)
upper

irfdata <- cbind(ts(lower, start=0), ts(response.inf, start=0), 
                 ts(upper, start=0))
colnames(irfdata) <- c("lower", "response", "upper")
irfdata
plot(irfdata, plot.type="single", 
     main="Response of Inflation to a Demand Shock", 
     lty=c(2,1,2))
