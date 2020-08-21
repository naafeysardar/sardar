dcpi.raw <- read.csv("dcpi.csv", header=TRUE)
dcpi <- ts(dcpi.raw[,2], start=c(1947,2), frequency=12)
u.raw <- read.csv("unrate.csv", header=TRUE)
u <- ts(u.raw[,2], start=c(1948,1), frequency=12)
library(tstools)
dataset <- ts.intersect(u, dcpi)

library(vars)
fit <- VAR(dataset, lag.max=13, ic="SC")
fit

E <- residuals(fit)
head(E, 5)

irfs <- chol(cov(E))
irfs

irfs %*% c(1,0)
irfs %*% c(0,1)

irf(fit, boot=FALSE)
irf(fit, boot=FALSE, impulse="u")
ds <- dataset[,c("dcpi","u")]

fit2 <- VAR(ds, lag.max = 13, ic="SC")
fit2

irf(fit2, boot=FALSE, impulse="u")
plot(irf(fit2, boot=FALSE, impulse="u", n.ahead=24))

plot(irf(fit2, boot=FALSE, impulse="u", response="dcpi", n.ahead=24))