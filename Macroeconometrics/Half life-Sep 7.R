# AR(1)

ysim <- arima.sim(model=list(ar=0.99), n=60)
 
plot(ysim)


# ARMA(3,2)

ysim1 <- arima.sim(model=list(ar=c(0.3, 0.2, 0.2), ma=c(0.1,0.4)), n=60)

# Estimate AR(1)

library(tstools)
fit <- tsreg(ysim, lags(ysim, 1))
fit

# Do the simulation many times
# Capture the second estimated coefficient each time

simoutput <- replicate(100, {
  ysim <- arima.sim(model = list(ar=0.99), n=60)
  fit <- tsreg(ysim, lags(ysim, 1))
  coef(fit)[2]
})
simoutput
median(simoutput)

plot(density(simoutput))

# Most of the weight centered around the true value.

mean(simoutput>0.99)
mean(simoutput<0.99)

range(simoutput)

# If we repeat this exercise, we get different answers. 
# So for that purpose we set seed. 

set.seed(200)
simoutput <- replicate(1000, {
  ysim <- arima.sim(model = list(ar=0.99), n=60)
  fit <- tsreg(ysim, lags(ysim, 1))
  coef(fit)[2]
})
simoutput
median(simoutput)

# Traditional Confirdence Interval 
se <- sd(simoutput)
se

CI_UB <- median(simoutput)+1.96*se
CI_UB

CI_LB <- median(simoutput)-1.96*se
CI_LB

mean(simoutput<0.78)

mean(simoutput>1.06)

# Empirical Distribution

quantile(simoutput, probs = c(0.025, 0.975))

# 97.5% value is below 1. 2.5% value is below 0.74. [0.74, 1.00] instead of [0.78, 1.06].

# Find CI of half life.
  
set.seed(200)
halflife <- replicate(1000, {
  ysim <- arima.sim(model = list(ar=0.9), n=60)
  fit <- tsreg(ysim, lags(ysim, 1))
  log(0.5) / log(coef(fit)[2])
})

quantile(halflife, probs=c(0.025,0.975))

median(halflife)
mean(halflife)





