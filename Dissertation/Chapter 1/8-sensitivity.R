library(systemfit)
library(msm)
library(tstools)
library(sandwich)

dataset.rpg <- ts.intersect(auerbach, gdp, rpg, inter)

# Define Initial Response Vectors
x0.rec <- c(10.00, 10.00)
x0.exp <- c(10.00, 0)

rfvar <- lapply(1:5, function(h) {
  rhs <- lags(dataset.rpg, h:(h+1), type="bylag")
  coefficients(tsreg(rpg, rhs))
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

rpg.rec <- (cumsum(unlist(rec)))
rpg.exp <- (cumsum(unlist(exp)))

t <- 0:5

#####################################################################################################################################################

dataset.pce <- ts.intersect(auerbach, gdp, pce, rpg, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(pce, 1:2), lags(rpg, 0:2), lags(inter, 0:2))
fit <- tsreg(pce, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

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

pce.rec <- 10*(pce.rec/rpg.rec)
pce.exp <- 10*(pce.exp/rpg.exp)

#####################################################################################################################################################

dataset.durable <- ts.intersect(auerbach, gdp, durable, rpg, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(durable, 1:2), lags(rpg, 0:2), lags(inter, 0:2))
fit <- tsreg(durable, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

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

durable.rec <- 10*(durable.rec/rpg.rec)
durable.exp <- 10*(durable.exp/rpg.exp)


#####################################################################################################################################################

dataset.nondurable <- ts.intersect(auerbach, gdp, nondurable, rpg, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(nondurable, 1:2), lags(rpg, 0:2), lags(inter, 0:2))
fit <- tsreg(nondurable, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

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

nondurable.rec <- 10*(nondurable.rec/rpg.rec)
nondurable.exp <- 10*(nondurable.exp/rpg.exp)


#####################################################################################################################################################

dataset.services <- ts.intersect(auerbach, gdp, services, rpg, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(services, 1:2), lags(rpg, 0:2), lags(inter, 0:2))
fit <- tsreg(services, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

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

services.rec <- 10*(services.rec/rpg.rec)
services.exp <- 10*(services.exp/rpg.exp)

par(mfrow=c(2,2))

plot(t, pce.rec, type = "l", col = "black", xlab = "Quarters", ylab = "Percent", ylim= c(-6, 6), main = "PCE"); abline(h=0)
points(t, pce.exp, type = "l", col = "black", lty = 2)

plot(t, durable.rec, type = "l", col = "black", xlab = "Quarters", ylab = "Percent", ylim= c(-6, 6), main = "Durables"); abline(h=0)
points(t, durable.exp, type = "l", col = "black", lty = 2)

plot(t, nondurable.rec, type = "l", col = "black", xlab = "Quarters", ylab = "Percent", ylim= c(-6, 6), main = "Nondurables"); abline(h=0)
points(t, nondurable.exp, type = "l", col = "black", lty = 2)

plot(t, services.rec, type = "l", col = "black", xlab = "Quarters", ylab = "Percent", ylim= c(-6, 6), main = "Services"); abline(h=0)
points(t, services.exp, type = "l", col = "black", lty = 2)

