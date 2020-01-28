library(tis, exclude = c("lags", "month", "quarter", "year"))
library(systemfit)
library(msm)
library(tstools)
library(sandwich)

par(mfrow=c(3,2))

dataset.gas <- ts.intersect(auerbach, gdp, gas, inter)

# Define Initial Response Vectors
x0.rec <- c(10.00, 10.00)
x0.exp <- c(10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.gas, h:(h+1), type="bylag")
  coefficients(tsreg(gas, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:5] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:5] * x0.rec)
})

exp <- unlist(list(10.00, response.exp))
rec <- unlist(list(10.00, response.rec))

gas.rec <- (cumsum(unlist(rec)))
gas.exp <- (cumsum(unlist(exp)))

#####################################################################################################################################################

dataset.pce <- ts.intersect(auerbach, gdp, pce, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(pce, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(pce, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.21*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.pce, h:(h+1), type="bylag")
  coefficients(tsreg(pce, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

pce.rec <- (cumsum(unlist(rec)))
pce.exp <- (cumsum(unlist(exp)))

pce.rec <- 10*(pce.rec/gas.rec)
pce.exp <- 10*(pce.exp/gas.exp)

d.pce <- pce.rec - pce.exp

x <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, pce, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.pce, h:(h+1), type="bylag")
  tsreg(pce, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                           -1.32*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                           -1.32*x13+18.96*x14+26.63*x15-0.67*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                                 -1.32*x13+18.96*x14+26.63*x15-0.67*x16+8.46*x17+15.99*x18-0.75*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.12*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                           -1.32*x13+18.96*x14+26.63*x15-0.67*x16+8.46*x17+15.99*x18-0.75*x19+9.83*x20+17.09*x21-1.08*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.pce.upper <- d.pce + 1.96*s
d.pce.lower <- d.pce - 1.96*s

#####################################################################################################################################################

dataset.durable <- ts.intersect(auerbach, gdp, durable, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(durable, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(durable, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.15*x1+10*x2)), theta.diff, v.diff)

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.durable, h:(h+1), type="bylag")
  coefficients(tsreg(durable, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

durable.rec <- (cumsum(unlist(rec)))
durable.exp <- (cumsum(unlist(exp)))

durable.rec <- 10*(durable.rec/gas.rec)
durable.exp <- 10*(durable.exp/gas.exp)

d.durable <- durable.rec - durable.exp

x <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, durable, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.durable, h:(h+1), type="bylag")
  tsreg(durable, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.84*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.84*x7+7.87*x8+15.30*x9-1.33*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.84*x7+7.87*x8+15.30*x9-1.33*x10+13.11*x11+20.51*x12
                           -1.86*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.84*x7+7.87*x8+15.30*x9-1.33*x10+13.11*x11+20.51*x12
                           -1.86*x13+18.96*x14+26.63*x15-0.89*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.84*x7+7.87*x8+15.30*x9-1.33*x10+13.11*x11+20.51*x12
                           -1.86*x13+18.96*x14+26.63*x15-0.89*x16+8.46*x17+15.99*x18-1.02*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.04*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.84*x7+7.87*x8+15.30*x9-1.33*x10+13.11*x11+20.51*x12
                           -1.86*x13+18.96*x14+26.63*x15-0.89*x16+8.46*x17+15.99*x18-1.02*x19+9.83*x20+17.09*x21-1.50*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.durable.upper <- d.durable + 1.96*s
d.durable.lower <- d.durable - 1.96*s

#####################################################################################################################################################

dataset.nondurable <- ts.intersect(auerbach, gdp, nondurable, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(nondurable, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(nondurable, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.24*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.nondurable, h:(h+1), type="bylag")
  coefficients(tsreg(nondurable, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

nondurable.rec <- (cumsum(unlist(rec)))
nondurable.exp <- (cumsum(unlist(exp)))

nondurable.rec <- 10*(nondurable.rec/gas.rec)
nondurable.exp <- 10*(nondurable.exp/gas.exp)

d.nondurable <- nondurable.rec - nondurable.exp

x <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, nondurable, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.nondurable, h:(h+1), type="bylag")
  tsreg(nondurable, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.58*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.58*x7+7.87*x8+15.30*x9-0.84*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.58*x7+7.87*x8+15.30*x9-0.84*x10+13.11*x11+20.51*x12
                           -1.14*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.58*x7+7.87*x8+15.30*x9-0.84*x10+13.11*x11+20.51*x12
                           -1.14*x13+18.96*x14+26.63*x15-0.61*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.58*x7+7.87*x8+15.30*x9-0.84*x10+13.11*x11+20.51*x12
                           -1.14*x13+18.96*x14+26.63*x15-0.61*x16+8.46*x17+15.99*x18-0.67*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.16*x1-1.03*x2+7.70*x3-0.33*x4+2.34*x5+11.10*x6-0.58*x7+7.87*x8+15.30*x9-0.84*x10+13.11*x11+20.51*x12
                           -1.14*x13+18.96*x14+26.63*x15-0.61*x16+8.46*x17+15.99*x18-0.67*x19+9.83*x20+17.09*x21-0.94*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.nondurable.upper <- d.nondurable + 1.96*s
d.nondurable.lower <- d.nondurable - 1.96*s

#####################################################################################################################################################

dataset.services <- ts.intersect(auerbach, gdp, services, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(services, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(services, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.24*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.services, h:(h+1), type="bylag")
  coefficients(tsreg(services, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

services.rec <- (cumsum(unlist(rec)))
services.exp <- (cumsum(unlist(exp)))

services.rec <- 10*(services.rec/gas.rec)
services.exp <- 10*(services.exp/gas.exp)

d.services <- services.rec - services.exp

x <- embed(cbind(auerbach, gdp, services, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, services, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.services, h:(h+1), type="bylag")
  tsreg(services, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.34*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.34*x4+2.34*x5+11.10*x6-0.62*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.34*x4+2.34*x5+11.10*x6-0.62*x7+7.87*x8+15.30*x9-0.92*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.34*x4+2.34*x5+11.10*x6-0.62*x7+7.87*x8+15.30*x9-0.92*x10+13.11*x11+20.51*x12
                           -1.25*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.34*x4+2.34*x5+11.10*x6-0.62*x7+7.87*x8+15.30*x9-0.92*x10+13.11*x11+20.51*x12
                           -1.25*x13+18.96*x14+26.63*x15-0.66*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.34*x4+2.34*x5+11.10*x6-0.62*x7+7.87*x8+15.30*x9-0.92*x10+13.11*x11+20.51*x12
                           -1.25*x13+18.96*x14+26.63*x15-0.66*x16+8.46*x17+15.99*x18-0.73*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.34*x4+2.34*x5+11.10*x6-0.62*x7+7.87*x8+15.30*x9-0.92*x10+13.11*x11+20.51*x12
                           -1.25*x13+18.96*x14+26.63*x15-0.66*x16+8.46*x17+15.99*x18-0.73*x19+9.83*x20+17.09*x21-1.03*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.services.upper <- d.services + 1.96*s
d.services.lower <- d.services - 1.96*s

#####################################################################################################################################################

dataset.motor <- ts.intersect(auerbach, gdp, motor, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(motor, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(motor, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.62*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.motor, h:(h+1), type="bylag")
  coefficients(tsreg(motor, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

motor.rec <- (cumsum(unlist(rec)))
motor.exp <- (cumsum(unlist(exp)))

motor.rec <- 10*(motor.rec/gas.rec)
motor.exp <- 10*(motor.exp/gas.exp)

d.motor <- motor.rec - motor.exp

x <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, motor, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.motor, h:(h+1), type="bylag")
  tsreg(motor, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3-0.91*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3-0.91*x4+2.34*x5+11.10*x6-1.71*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3-0.91*x4+2.34*x5+11.10*x6-1.71*x7+7.87*x8+15.30*x9-2.54*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3-0.91*x4+2.34*x5+11.10*x6-1.71*x7+7.87*x8+15.30*x9-2.54*x10+13.11*x11+20.51*x12
                           -3.49*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3-0.91*x4+2.34*x5+11.10*x6-1.71*x7+7.87*x8+15.30*x9-2.54*x10+13.11*x11+20.51*x12
                           -3.49*x13+18.96*x14+26.63*x15-1.81*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3-0.91*x4+2.34*x5+11.10*x6-1.71*x7+7.87*x8+15.30*x9-2.54*x10+13.11*x11+20.51*x12
                           -3.49*x13+18.96*x14+26.63*x15-1.81*x16+8.46*x17+15.99*x18-2.01*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.38*x1-1.03*x2+7.70*x3-0.91*x4+2.34*x5+11.10*x6-1.71*x7+7.87*x8+15.30*x9-2.54*x10+13.11*x11+20.51*x12
                           -3.49*x13+18.96*x14+26.63*x15-1.81*x16+8.46*x17+15.99*x18-2.01*x19+9.83*x20+17.09*x21-2.86*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.motor.upper <- d.motor + 1.96*s
d.motor.lower <- d.motor - 1.96*s

#####################################################################################################################################################

dataset.furnishing <- ts.intersect(auerbach, gdp, furnishing, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(furnishing, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(furnishing, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((+0.04*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.furnishing, h:(h+1), type="bylag")
  coefficients(tsreg(furnishing, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

furnishing.rec <- (cumsum(unlist(rec)))
furnishing.exp <- (cumsum(unlist(exp)))

furnishing.rec <- 10*(furnishing.rec/gas.rec)
furnishing.exp <- 10*(furnishing.exp/gas.exp)

d.furnishing <- furnishing.rec - furnishing.exp

x <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, furnishing, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.furnishing, h:(h+1), type="bylag")
  tsreg(furnishing, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3-0.09*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3-0.09*x4+2.34*x5+11.10*x6-0.39*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3-0.09*x4+2.34*x5+11.10*x6-0.39*x7+7.87*x8+15.30*x9-0.67*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3-0.09*x4+2.34*x5+11.10*x6-0.39*x7+7.87*x8+15.30*x9-0.67*x10+13.11*x11+20.51*x12
                           -0.98*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3-0.09*x4+2.34*x5+11.10*x6-0.39*x7+7.87*x8+15.30*x9-0.67*x10+13.11*x11+20.51*x12
                           -0.98*x13+18.96*x14+26.63*x15-0.42*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3-0.09*x4+2.34*x5+11.10*x6-0.39*x7+7.87*x8+15.30*x9-0.67*x10+13.11*x11+20.51*x12
                           -0.98*x13+18.96*x14+26.63*x15-0.42*x16+8.46*x17+15.99*x18-0.49*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(+0.09*x1-1.03*x2+7.70*x3-0.09*x4+2.34*x5+11.10*x6-0.39*x7+7.87*x8+15.30*x9-0.67*x10+13.11*x11+20.51*x12
                           -0.98*x13+18.96*x14+26.63*x15-0.42*x16+8.46*x17+15.99*x18-0.49*x19+9.83*x20+17.09*x21-0.77*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.furnishing.upper <- d.furnishing + 1.96*s
d.furnishing.lower <- d.furnishing - 1.96*s

#####################################################################################################################################################

dataset.otherdurable <- ts.intersect(auerbach, gdp, otherdurable, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(otherdurable, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(otherdurable, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((+0.24*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.otherdurable, h:(h+1), type="bylag")
  coefficients(tsreg(otherdurable, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

otherdurable.rec <- (cumsum(unlist(rec)))
otherdurable.exp <- (cumsum(unlist(exp)))

otherdurable.rec <- 10*(otherdurable.rec/gas.rec)
otherdurable.exp <- 10*(otherdurable.exp/gas.exp)

d.otherdurable <- otherdurable.rec - otherdurable.exp

x <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, otherdurable, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.otherdurable, h:(h+1), type="bylag")
  tsreg(otherdurable, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3+0.17*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3+0.17*x4+2.34*x5+11.10*x6+0.04*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3+0.17*x4+2.34*x5+11.10*x6+0.04*x7+7.87*x8+15.30*x9-0.06*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3+0.17*x4+2.34*x5+11.10*x6+0.04*x7+7.87*x8+15.30*x9-0.06*x10+13.11*x11+20.51*x12
                           -0.15*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3+0.17*x4+2.34*x5+11.10*x6+0.04*x7+7.87*x8+15.30*x9-0.06*x10+13.11*x11+20.51*x12
                           -0.15*x13+18.96*x14+26.63*x15+0.03*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3+0.17*x4+2.34*x5+11.10*x6+0.04*x7+7.87*x8+15.30*x9-0.06*x10+13.11*x11+20.51*x12
                           -0.15*x13+18.96*x14+26.63*x15+0.03*x16+8.46*x17+15.99*x18+0.00*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(+0.24*x1-1.03*x2+7.70*x3+0.17*x4+2.34*x5+11.10*x6+0.04*x7+7.87*x8+15.30*x9-0.06*x10+13.11*x11+20.51*x12
                           -0.15*x13+18.96*x14+26.63*x15+0.03*x16+8.46*x17+15.99*x18+0.00*x19+9.83*x20+17.09*x21-0.08*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.otherdurable.upper <- d.otherdurable + 1.96*s
d.otherdurable.lower <- d.otherdurable - 1.96*s

#####################################################################################################################################################

dataset.food <- ts.intersect(auerbach, gdp, food, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(food, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(food, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.37*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.food, h:(h+1), type="bylag")
  coefficients(tsreg(food, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

food.rec <- (cumsum(unlist(rec)))
food.exp <- (cumsum(unlist(exp)))

food.rec <- 10*(food.rec/gas.rec)
food.exp <- 10*(food.exp/gas.exp)

d.food <- food.rec - food.exp

x <- embed(cbind(auerbach, gdp, food, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, food, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.food, h:(h+1), type="bylag")
  tsreg(food, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3-0.43*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3-0.43*x4+2.34*x5+11.10*x6-0.63*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3-0.43*x4+2.34*x5+11.10*x6-0.63*x7+7.87*x8+15.30*x9-0.87*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3-0.43*x4+2.34*x5+11.10*x6-0.63*x7+7.87*x8+15.30*x9-0.87*x10+13.11*x11+20.51*x12
                           -1.15*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3-0.43*x4+2.34*x5+11.10*x6-0.63*x7+7.87*x8+15.30*x9-0.87*x10+13.11*x11+20.51*x12
                           -1.15*x13+18.96*x14+26.63*x15-0.66*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3-0.43*x4+2.34*x5+11.10*x6-0.63*x7+7.87*x8+15.30*x9-0.87*x10+13.11*x11+20.51*x12
                           -1.15*x13+18.96*x14+26.63*x15-0.66*x16+8.46*x17+15.99*x18-0.72*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.27*x1-1.03*x2+7.70*x3-0.43*x4+2.34*x5+11.10*x6-0.63*x7+7.87*x8+15.30*x9-0.87*x10+13.11*x11+20.51*x12
                           -1.15*x13+18.96*x14+26.63*x15-0.66*x16+8.46*x17+15.99*x18-0.72*x19+9.83*x20+17.09*x21-0.97*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.food.upper <- d.food + 1.96*s
d.food.lower <- d.food - 1.96*s


#####################################################################################################################################################

dataset.clothing <- ts.intersect(auerbach, gdp, clothing, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(clothing, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(clothing, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((+0.04*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.clothing, h:(h+1), type="bylag")
  coefficients(tsreg(clothing, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

clothing.rec <- (cumsum(unlist(rec)))
clothing.exp <- (cumsum(unlist(exp)))

clothing.rec <- 10*(clothing.rec/gas.rec)
clothing.exp <- 10*(clothing.exp/gas.exp)

d.clothing <- clothing.rec - clothing.exp

x <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, clothing, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.clothing, h:(h+1), type="bylag")
  tsreg(clothing, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3+0.04*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3+0.04*x4+2.34*x5+11.10*x6-0.24*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3-0.04*x4+2.34*x5+11.10*x6-0.24*x7+7.87*x8+15.30*x9-0.42*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3-0.04*x4+2.34*x5+11.10*x6-0.24*x7+7.87*x8+15.30*x9-0.42*x10+13.11*x11+20.51*x12
                           -0.61*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3-0.04*x4+2.34*x5+11.10*x6-0.24*x7+7.87*x8+15.30*x9-0.42*x10+13.11*x11+20.51*x12
                           -0.61*x13+18.96*x14+26.63*x15-0.26*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3-0.04*x4+2.34*x5+11.10*x6-0.24*x7+7.87*x8+15.30*x9-0.42*x10+13.11*x11+20.51*x12
                           -0.61*x13+18.96*x14+26.63*x15-0.26*x16+8.46*x17+15.99*x18-0.30*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(+0.07*x1-1.03*x2+7.70*x3-0.04*x4+2.34*x5+11.10*x6-0.24*x7+7.87*x8+15.30*x9-0.42*x10+13.11*x11+20.51*x12
                           -0.61*x13+18.96*x14+26.63*x15-0.26*x16+8.46*x17+15.99*x18-0.30*x19+9.83*x20+17.09*x21-0.48*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.clothing.upper <- d.clothing + 1.96*s
d.clothing.lower <- d.clothing - 1.96*s

#####################################################################################################################################################

dataset.othernondurable <- ts.intersect(auerbach, gdp, othernondurable, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(othernondurable, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(othernondurable, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.20*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.othernondurable, h:(h+1), type="bylag")
  coefficients(tsreg(othernondurable, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

othernondurable.rec <- (cumsum(unlist(rec)))
othernondurable.exp <- (cumsum(unlist(exp)))

othernondurable.rec <- 10*(othernondurable.rec/gas.rec)
othernondurable.exp <- 10*(othernondurable.exp/gas.exp)

d.othernondurable <- othernondurable.rec - othernondurable.exp

x <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, othernondurable, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.othernondurable, h:(h+1), type="bylag")
  tsreg(othernondurable, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3-0.32*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3-0.32*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3-0.32*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.97*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3-0.32*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.97*x10+13.11*x11+20.51*x12
                           -1.33*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3-0.32*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.97*x10+13.11*x11+20.51*x12
                           -1.33*x13+18.96*x14+26.63*x15-0.68*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3-0.32*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.97*x10+13.11*x11+20.51*x12
                           -1.33*x13+18.96*x14+26.63*x15-0.68*x16+8.46*x17+15.99*x18-0.76*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.11*x1-1.03*x2+7.70*x3-0.32*x4+2.34*x5+11.10*x6-0.64*x7+7.87*x8+15.30*x9-0.97*x10+13.11*x11+20.51*x12
                           -1.33*x13+18.96*x14+26.63*x15-0.68*x16+8.46*x17+15.99*x18-0.76*x19+9.83*x20+17.09*x21-1.09*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.othernondurable.upper <- d.othernondurable + 1.96*s
d.othernondurable.lower <- d.othernondurable - 1.96*s

##############################################################################################################################################################

dataset.housing <- ts.intersect(auerbach, gdp, housing, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(housing, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(housing, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.24*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.housing, h:(h+1), type="bylag")
  coefficients(tsreg(housing, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

housing.rec <- (cumsum(unlist(rec)))
housing.exp <- (cumsum(unlist(exp)))

housing.rec <- 10*(housing.rec/gas.rec)
housing.exp <- 10*(housing.exp/gas.exp)

d.housing <- housing.rec - housing.exp

x <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, housing, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.housing, h:(h+1), type="bylag")
  tsreg(housing, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.65*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.65*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.65*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                           -1.32*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.65*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                           -1.32*x13+18.96*x14+26.63*x15-0.69*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.65*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                           -1.32*x13+18.96*x14+26.63*x15-0.69*x16+8.46*x17+15.99*x18-0.76*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.15*x1-1.03*x2+7.70*x3-0.35*x4+2.34*x5+11.10*x6-0.65*x7+7.87*x8+15.30*x9-0.96*x10+13.11*x11+20.51*x12
                           -1.32*x13+18.96*x14+26.63*x15-0.69*x16+8.46*x17+15.99*x18-0.76*x19+9.83*x20+17.09*x21-1.08*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.housing.upper <- d.housing + 1.96*s
d.housing.lower <- d.housing - 1.96*s

##############################################################################################################################################################

dataset.transport <- ts.intersect(auerbach, gdp, transport, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(transport, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(transport, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.11*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.transport, h:(h+1), type="bylag")
  coefficients(tsreg(transport, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

transport.rec <- (cumsum(unlist(rec)))
transport.exp <- (cumsum(unlist(exp)))

transport.rec <- 10*(transport.rec/gas.rec)
transport.exp <- 10*(transport.exp/gas.exp)

d.transport <- transport.rec - transport.exp

x <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, transport, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.transport, h:(h+1), type="bylag")
  tsreg(transport, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3-0.12*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3-0.12*x4+2.34*x5+11.10*x6-0.17*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3-0.12*x4+2.34*x5+11.10*x6-0.17*x7+7.87*x8+15.30*x9-0.23*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3-0.12*x4+2.34*x5+11.10*x6-0.17*x7+7.87*x8+15.30*x9-0.23*x10+13.11*x11+20.51*x12
                           -0.30*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3-0.12*x4+2.34*x5+11.10*x6-0.17*x7+7.87*x8+15.30*x9-0.23*x10+13.11*x11+20.51*x12
                           -0.30*x13+18.96*x14+26.63*x15-0.18*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3-0.12*x4+2.34*x5+11.10*x6-0.17*x7+7.87*x8+15.30*x9-0.23*x10+13.11*x11+20.51*x12
                           -0.30*x13+18.96*x14+26.63*x15-0.18*x16+8.46*x17+15.99*x18-0.19*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.08*x1-1.03*x2+7.70*x3-0.12*x4+2.34*x5+11.10*x6-0.17*x7+7.87*x8+15.30*x9-0.23*x10+13.11*x11+20.51*x12
                           -0.30*x13+18.96*x14+26.63*x15-0.18*x16+8.46*x17+15.99*x18-0.76*x19+9.83*x20+17.09*x21-1.08*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.transport.upper <- d.transport + 1.96*s
d.transport.lower <- d.transport - 1.96*s

##############################################################################################################################################################

dataset.other.services <- ts.intersect(auerbach, gdp, other.services, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(other.services, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(other.services, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.45*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.other.services, h:(h+1), type="bylag")
  coefficients(tsreg(other.services, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

other.services.rec <- (cumsum(unlist(rec)))
other.services.exp <- (cumsum(unlist(exp)))

other.services.rec <- 10*(other.services.rec/gas.rec)
other.services.exp <- 10*(other.services.exp/gas.exp)

d.other.services <- other.services.rec - other.services.exp

x <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)

x <- as.data.frame(x)

one.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 6:15]
two.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 11:20]
three.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 16:25] 
four.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 21:30] 
five.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 26:35] 
six.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 31:40] 
seven.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 36:45]
eight.step <- embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 3] ~ embed(cbind(auerbach, gdp, other.services, gas, inter), 10)[, 41:50]

system <- list(one.step, two.step, three.step, four.step, five.step, six.step, seven.step, eight.step)

fitsur <- systemfit(system, data = x)

v <- matrix(vcov(fitsur), nrow = 88, ncol = 88)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.other.services, h:(h+1), type="bylag")
  tsreg(other.services, rhs)
})

cov.mats <- lapply(rfvar, function(z) {
  NeweyWest(z)
})


replaceDiag <- function(a, mats, start=1) {
  if (length(mats) == 0) {
    return(a)
  } else {
    a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
      mats[[1]]
    Recall(a, mats[-1], start+nrow(mats[[1]]))
  }
}

x.newey <- replaceDiag(v, cov.mats)

# Standard Errors

# 1-Step #
theta <- coefficients(fitsur)[4:6]
v1 <- as.matrix(x.newey[c(4, 5, 6), c(4, 5, 6)])
std.err.1 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3), theta, v1)
# 2-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17)]
v2 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17), c(4, 5, 6, 15, 16, 17)])
std.err.2 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3-0.57*x4+2.34*x5+11.10*x6), theta, v2)
# 3-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28)]
v3 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28), c(4, 5, 6, 15, 16, 17, 26, 27, 28)])
std.err.3 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3-0.57*x4+2.34*x5+11.10*x6-0.94*x7+7.87*x8+15.30*x9), theta, v3)
# 4-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)]
v4 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39)])
std.err.4 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3-0.57*x4+2.34*x5+11.10*x6-0.94*x7+7.87*x8+15.30*x9-1.35*x10+13.11*x11+20.51*x12), theta, v4)
# 5-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)]
v5 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50), c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50)])
std.err.5 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3-0.57*x4+2.34*x5+11.10*x6-0.94*x7+7.87*x8+15.30*x9-1.35*x10+13.11*x11+20.51*x12
                           -1.81*x13+18.96*x14+26.63*x15), theta, v5)
# 6-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)]
v6 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61)])
std.err.6 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3-0.57*x4+2.34*x5+11.10*x6-0.94*x7+7.87*x8+15.30*x9-1.35*x10+13.11*x11+20.51*x12
                           -1.81*x13+18.96*x14+26.63*x15-0.99*x16+8.46*x17+15.99*x18), theta, v6)
# 7-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)]
v7 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72)])
std.err.7 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3-0.57*x4+2.34*x5+11.10*x6-0.94*x7+7.87*x8+15.30*x9-1.35*x10+13.11*x11+20.51*x12
                           -1.81*x13+18.96*x14+26.63*x15-0.99*x16+8.46*x17+15.99*x18-1.09*x19+9.83*x20+17.09*x21), theta, v7)
# 8-Step #
theta <- coefficients(fitsur)[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)]
v8 <- as.matrix(x.newey[c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83), 
                        c(4, 5, 6, 15, 16, 17, 26, 27, 28, 37, 38, 39, 48, 49, 50, 59, 60, 61, 70, 71, 72, 81, 82, 83)])
std.err.8 <- deltamethod(~(-0.31*x1-1.03*x2+7.70*x3-0.57*x4+2.34*x5+11.10*x6-0.94*x7+7.87*x8+15.30*x9-1.35*x10+13.11*x11+20.51*x12
                           -1.81*x13+18.96*x14+26.63*x15-0.99*x16+8.46*x17+15.99*x18-1.09*x19+9.83*x20+17.09*x21-1.50*x22+15.05*x23+22.59*x24), theta, v8)

s <- unlist(list(std.err.0, std.err.1, std.err.2, std.err.3, std.err.4, std.err.5))

d.other.services.upper <- d.other.services + 1.96*s
d.other.services.lower <- d.other.services - 1.96*s

##############################################################################################################################################################

par(mfrow=c(2,2))

plot(pce.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "PCE"); abline(h=0)
points(pce.exp , col = "black", type="l", lty=2)

plot(durable.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "Durables"); abline(h=0)
points(durable.exp, type = "l", col = "black", lty=2)

plot(nondurable.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "Nondurables"); abline(h=0)
points(nondurable.exp, type = "l", col = "black", lty=2)

plot(services.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "Services"); abline(h=0)
points(services.exp, type = "l", col = "black", lty=2)

par(mfrow=c(2,2))

plot(d.pce, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-6, 2), main = "PCE"); abline(h=0)
points(d.pce.upper , col = "black", type="l", lty=2)
points(d.pce.lower, col = "black", type="l", lty=2)

plot(d.durable, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-15, 5), main = "Durables"); abline(h=0)
points(d.durable.upper , col = "black", type="l", lty=2)
points(d.durable.lower, col = "black", type="l", lty=2)

plot(d.nondurable, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 2), main = "Nondurables"); abline(h=0)
points(d.nondurable.upper , col = "black", type="l", lty=2)
points(d.nondurable.lower, col = "black", type="l", lty=2)

plot(d.services, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 2), main = "Services"); abline(h=0)
points(d.services.upper , col = "black", type="l", lty=2)
points(d.services.lower, col = "black", type="l", lty=2)

par(mfrow=c(2,3))

plot(motor.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-20, 10), main = "Motor Vehicles"); abline(h=0)
points(motor.exp, type = "l", col = "black", lty=2)

plot(d.motor, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-20, 10), main = "Motor Vehicles"); abline(h=0)
points(d.motor.upper , col = "black", type="l", lty=2)
points(d.motor.lower, col = "black", type="l", lty=2)

plot(furnishing.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "Furnishing"); abline(h=0)
points(furnishing.exp, type = "l", col = "black", lty=2)

plot(d.furnishing, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "Furnishing"); abline(h=0)
points(d.furnishing.upper , col = "black", type="l", lty=2)
points(d.furnishing.lower, col = "black", type="l", lty=2)

plot(otherdurable.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "Other Durables"); abline(h=0)
points(otherdurable.exp, type = "l", col = "black", lty=2)

plot(d.otherdurable, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 5), main = "Other Durables"); abline(h=0)
points(d.otherdurable.upper , col = "black", type="l", lty=2)
points(d.otherdurable.lower, col = "black", type="l", lty=2)

plot(food.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 2), main = "Food"); abline(h=0)
points(food.exp, type = "l", col = "black", lty=2)

plot(d.food, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 2), main = "Food"); abline(h=0)
points(d.food.upper , col = "black", type="l", lty=2)
points(d.food.lower, col = "black", type="l", lty=2)

plot(clothing.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-6, 2), main = "Clothing"); abline(h=0)
points(clothing.exp, type = "l", col = "black", lty=2)

plot(d.clothing, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-6, 2), main = "Clothing"); abline(h=0)
points(d.clothing.upper , col = "black", type="l", lty=2)
points(d.clothing.lower, col = "black", type="l", lty=2)

plot(othernondurable.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 2), main = "Other Nondurables"); abline(h=0)
points(othernondurable.exp, type = "l", col = "black", lty=2)

plot(d.othernondurable, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-10, 2), main = "Other Nondurables"); abline(h=0)
points(d.othernondurable.upper , col = "black", type="l", lty=2)
points(d.othernondurable.lower, col = "black", type="l", lty=2)

plot(housing.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 5), main = "Housing & Utilities"); abline(h=0)
points(housing.exp, type = "l", col = "black", lty=2)

plot(d.housing, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 5), main = "Housing & Utilities"); abline(h=0)
points(d.housing.upper , col = "black", type="l", lty=2)
points(d.housing.lower, col = "black", type="l", lty=2)

plot(transport.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 5), main = "Transportation"); abline(h=0)
points(transport.exp, type = "l", col = "black", lty=2)

plot(d.transport, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 5), main = "Transportation"); abline(h=0)
points(d.transport.upper , col = "black", type="l", lty=2)
points(d.transport.lower, col = "black", type="l", lty=2)

plot(other.services.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 5), main = "Other Services"); abline(h=0)
points(other.services.exp, type = "l", col = "black", lty=2)

plot(d.other.services, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-5, 5), main = "Other Services"); abline(h=0)
points(d.other.services.upper , col = "black", type="l", lty=2)
points(d.other.services.lower, col = "black", type="l", lty=2)

#####################################################################################################################################################

dataset.saving <- ts.intersect(auerbach, gdp, saving, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(saving, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(saving, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

(10/gas.rec)*rec.contemp - (10/gas.exp)*exp.contemp

theta.diff <- coefficients(fit)[c(10, 13)]
v.diff <- vcov(fit)[c(10, 13), c(10, 13)]
std.err.0 <- deltamethod(~((-0.21*x1+10*x2)), theta.diff, v.diff)

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.saving, h:(h+1), type="bylag")
  coefficients(tsreg(saving, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

saving.exp <- (unlist(list(exp.contemp, response.exp)))
saving.rec <- (unlist(list(rec.contemp, response.rec)))

saving.rec <- 10*(saving.rec/gas.rec)
saving.exp <- 10*(saving.exp/gas.exp)










dataset.income <- ts.intersect(auerbach, gdp, income, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(income, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(income, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.income, h:(h+1), type="bylag")
  coefficients(tsreg(income, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

income.rec <- (cumsum(unlist(rec)))
income.exp <- (cumsum(unlist(exp)))

income.rec <- 10*(income.rec/gas.rec)
income.exp <- 10*(income.exp/gas.exp)










dataset.savings <- ts.intersect(auerbach, gdp, savings, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(savings, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(savings, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp


# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.savings, h:(h+1), type="bylag")
  coefficients(tsreg(savings, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

exp <- list(exp.contemp, response.exp)
rec <- list(rec.contemp, response.rec)

savings.rec <- (cumsum(unlist(rec)))
savings.exp <- (cumsum(unlist(exp)))

savings.rec <- 10*(savings.rec/gas.rec)
savings.exp <- 10*(savings.exp/gas.exp)
