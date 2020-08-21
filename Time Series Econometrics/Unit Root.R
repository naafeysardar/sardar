library(urca)
data("Raotbl3")
attach(Raotbl3)

## Work with log of consumption & log of income ##

lc <- ts(lc, start=c(1966,4), frequency=4)
li <- ts(li, start=c(1966,4), frequency=4)

plot(lc)
plot(li)

## There is going to be linear time trend and constant ##

## type none means no deterministic or time trend ##

lc.df <- ur.df(y=lc, type="none", lags=0)
summary(lc.df)

lc.df <- ur.df(y=lc, type="drift", lags=0)
summary(lc.df)

lc.df <- ur.df(y=lc, type="trend", lags=0)
summary(lc.df)

## Haven't rejected H0 in any of the three cases so take difference ##

lc.df <- ur.df(y=lc, type="trend", selectlags = "AIC")
summary(lc.df)

lc.df <- ur.df(y=lc, type="trend", selectlags = "BIC")
summary(lc.df)

pp.lc <- ur.pp(lc)
summary(pp.lc)

## Test for Conitegration ##

library(tstools)

fit <- tsreg(lc, li)
fit

z <- fit$resids
plot(z)


summary(ur.df(z, type="drift", lags=4))

fit2 <- tsreg(li, lc)

