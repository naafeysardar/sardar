library(tstools)
data.raw <- read.csv("u-inf.csv", header=TRUE)
dataset <- ts(data.raw, start=c(1948,1), frequency=12)
colnames(dataset) <- c("u", "inf")
u <- dataset[, "u"]
inf <- dataset[, "inf"]

# Estimate RF VAR
fit.u <- tsreg(u, lags(dataset,1))
fit.inf <- tsreg(inf, lags(dataset,1))
resu <- fit.u$resids
resinf <- fit.inf$resids

# Regime 1
resu.regime1 <- window(resu, end=c(1987,8))
resinf.regime1 <- window(resinf, end=c(1987,8))
var.u1 <- var(resu.regime1)
var.inf1 <- var(resinf.regime1)
cov.1 <- cov(resu.regime1, resinf.regime1)

# Regime 2
resu.regime2 <- window(resu, start=c(1987,9))
resinf.regime2 <- window(resinf, start=c(1987,9))
var.u2 <- var(resu.regime2)
var.inf2 <- var(resinf.regime2)
cov.2 <- cov(resu.regime2, resinf.regime2)
var.u1
var.inf1
cov.1
var.u2
var.inf1
cov.2

# Set up the obj function to solve the moment equations
dev <- function(par) {
  b <- par[1]
  c <- par[2]
  su1 <- par[3]
  su2 <- par[4]
  sinf1 <- par[5]
  sinf2 <- par[6]
  dev1 <- 0.05906124 - su1 - b^2*sinf1
  dev2 <- 0.08342547 - c^2*su1 - sinf1
  dev3 <- -0.006096401 - c*su1 - b*sinf1
  dev4 <- 0.02367125 - su2 - b^2*sinf2
  dev5 <- 0.05656814 - c^2*su2 - sinf2
  dev6 <- -0.0002381365 - c*su2 - b*sinf2
  return(c(dev1, dev2, dev3, dev4, dev5, dev6))
}
dev(c(0, 0, 0.2, 0.2, 0.1, 0.2))
library(nleqslv)
nleqslv(c(0, 0, 0.2, 0.2, 0.1, 0.2), dev)

ufit <- fit.u$fitted
inffit <- fit.inf$fitted

# Draw from pick distribution
z <- sample(c(-1,1), size=length(resu), replace=TRUE)

# Create simulated data
u.sim <- ufit + z*resu
inf.sim <- inffit + z*resinf

# Calculate IRFs for the simulated data
dataset.sim <- ts.combine(u.sim, inf.sim)
library(vars)
varfit.sim <- VAR(dataset.sim, lag.max=13, ic="SC")
irf.sim <- irf(varfit.sim, boot=FALSE, n.ahead=12,
               impulse="u.sim", response="inf.sim")
plot(irf.ts(irf.sim))

set.seed(100)
irfboot <- replicate(100, {
  z <- sample(c(-1,1), size=length(resu), replace=TRUE)
  u.sim <- ufit + z*resu
  inf.sim <- inffit + z*resinf
  dataset.sim <- ts.combine(u.sim, inf.sim)
  varfit.sim <- VAR(dataset.sim, lag.max=13, ic="SC")
  irf.sim <- irf(varfit.sim, boot=FALSE, n.ahead=12,
                 impulse="u.sim", response="inf.sim")
  irf.ts(irf.sim)
}, simplify="array")
class(irfboot)
dim(irfboot)
# 95% band for IRF0
irf0 <- sort(irfboot[1,])
quantile(irf0, probs=c(0.025, 0.975))

irf1 <- sort(irfboot[2,])
quantile(irf1, probs=c(0.025, 0.975))
