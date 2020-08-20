library(systemfit)
library(msm)
library(tstools)
library(sandwich)
library(seasonal)

t <- 0:4

par(mfrow=c(2,2))

dataset.gas <- ts.intersect(auerbach, gdp, gas, inter)

# Define Initial Response Vectors
x0.rec <- c(10.00, 10.00)
x0.exp <- c(10.00, 0)

rfvar <- lapply(1:4, function(h) {
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

gas.exp <- cumsum(unlist(list(10.00, response.exp)))
gas.rec <- cumsum(unlist(list(10.00, response.rec)))

plot(t, gas.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(0, 20), main = "Gasoline Expenditures"); abline(h=0)
points(t, gas.exp, type = "l", col = "black", lty = 2)

#####################################################################################################################################################

dataset.rpg <- ts.intersect(auerbach, gdp, rpg, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(rpg, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(rpg, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:4, function(h) {
  rhs <- lags(dataset.rpg, h:(h+1), type="bylag")
  coefficients(tsreg(rpg, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

rpg.exp <- cumsum(unlist(list(exp.contemp, response.exp)))
rpg.rec <- cumsum(unlist(list(rec.contemp, response.rec)))

plot(t, rpg.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(0, 20), main = "Real Price of Gas"); abline(h=0)
points(t, rpg.exp, type = "l", col = "black", lty = 2)

#####################################################################################################################################################

mydata <- read.csv("gasoline-primary.csv", header=TRUE)
monthly <- ts(mydata, start=c(1973,1),frequency=12)

residential <- monthly[,"residential"]
residential <- aggregate(residential, nfrequency=4)

adjresidential <- seas(residential)

res <- final(adjresidential)

res <- 100*(res/stats::lag(res,-1) - 1)

dataset.res <- ts.intersect(auerbach, gdp, res, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(res, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(res, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:4, function(h) {
  rhs <- lags(dataset.res, h:(h+1), type="bylag")
  coefficients(tsreg(res, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

residential.exp <- cumsum(unlist(list(exp.contemp, response.exp)))
residential.rec <- cumsum(unlist(list(rec.contemp, response.rec)))

plot(t, residential.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-2, 2), main = "Residential Consumption"); abline(h=0)
points(t, residential.exp, type = "l", col = "black", lty = 2)

#####################################################################################################################################################

transport <- monthly[,"transport"]
transport <- aggregate(transport, nfrequency=4)

adjtransport <- seas(transport)

trans <- final(adjtransport)

trans <- 100*(trans/stats::lag(trans,-1) - 1)

dataset.trans <- ts.intersect(auerbach, gdp, trans, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(trans, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(trans, rhs)

# Contemporaneous transponse

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

# Define Initial transponse Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:4, function(h) {
  rhs <- lags(dataset.trans, h:(h+1), type="bylag")
  coefficients(tsreg(trans, rhs))
})

# Impulse transponse in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse transponse in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

trans.exp <- cumsum(unlist(list(exp.contemp, response.exp)))
trans.rec <- cumsum(unlist(list(rec.contemp, response.rec)))

plot(t, trans.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-2, 2), main = "Transport Consumption"); abline(h=0)
points(t, trans.exp, type = "l", col = "black", lty = 2)


#####################################################################################################################################################

par(mfrow=c(1,1))

dataset.miles <- ts.intersect(auerbach, gdp, miles, gas, inter)

rhs <- ts.intersect(lags(auerbach, 0:2), lags(gdp, 0:2), lags(miles, 1:2), lags(gas, 0:2), lags(inter, 0:2))
fit <- tsreg(miles, rhs)

# Contemporaneous Response

rec.contemp <- (coefficients(fit)[10]+ coefficients(fit)[13])*10.00
exp.contemp <- (coefficients(fit)[10])*10.00
rec.contemp
exp.contemp

# Define Initial Response Vectors
x0.rec <- c(rec.contemp, 10.00, 10.00)
x0.exp <- c(exp.contemp, 10.00, 0)

rfvar <- lapply(1:4, function(h) {
  rhs <- lags(dataset.miles, h:(h+1), type="bylag")
  coefficients(tsreg(miles, rhs))
})

# Impulse Response in Expansion
response.exp <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.exp)
})
# Impulse Response in Recession
response.rec <- lapply(rfvar, function(b) {
  sum(b[4:6] * x0.rec)
})

miles.exp <- cumsum(unlist(list(exp.contemp, response.exp)))
miles.rec <- cumsum(unlist(list(rec.contemp, response.rec)))

plot(t, miles.rec, type = "l", col = "black", xlab = "Quarters", ylab = "% Change", ylim= c(-2, 2), main = "Miles Travelled"); abline(h=0)
points(t, miles.exp, type = "l", col = "black", lty = 2)
