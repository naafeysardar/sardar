project.raw <- read.csv("Book .csv", header=TRUE)
dataset <- ts(project.raw, start=c(1970,1), frequency=12)
y <- dataset[,"y"]
nopi <- dataset[,"nopi"]
inter <- dataset[,"inter"]

dataset

library(vars)
library(tstools)
rfvar <- lapply(1:12, function(h) {
  k <- 15
  rhs <- lags(dataset, h:(h+k-1), type="bylag")
  coefficients(tsreg(y, rhs))
})

x0.exp <- c(0, 1, 1)
x0.rec <- c(0, 1, 0)

response.exp <- apply(rfvar, function(b) {
  sum(b[2:4] * x0.exp)
})

response.rec <- apply(rfvar, function(b) {
  sum(b[2:4] * x0.rec)
})

Map(function(irf1, irf2) { irf1-irf2 }, response.exp, response.rec)
