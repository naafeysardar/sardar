data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]

library(tstools)

possible.breaks <- dates(c(1958,7), c(2007,11), 12)

f.calc <- function(break.date){
  d <- time.dummy(u, break.date, "end")
  rhs <- ts.intersect(lags(u, 1:2), d, d*lags(u, 1:2))
  rhs1 <- rhs[, 1:2]
  fit.large <- tsreg(u, rhs)
  fit.small <- tsreg(u, rhs1)
  test <- waldtest(fit.large, fit.small)
  return((test[2,3]))                
}

# Calculate F statistics for all possible break dates 

ftest.stats <- sapply(possible.breaks, f.calc)   
ftest.stats

# Sup F tells us the largest of the F statistics

supF <- max(ftest.stats)
supF

which.max(ftest.stats)

# Largest F Statistic occurs at 85th observation

possible.breaks[85]

# First break date occurs at July 1965
# Second break is from 212 to 593.

dates2 <- dates(c(1976,2), c(2007,11), 12)

f.calc2 <- function(break2){
  d1 <- time.dummy(u, c(1965,7), "end")
  d2 <- time.dummy(u, break2, "end")
  rhs <- ts.intersect(lags(u, 1:2), d1, d1*lags(u, 1:2), d2, d2*lags(u, 1:2))
  
  rhs1 <- rhs[, 1:2]
  fit.large <- tsreg(u, rhs)
  fit.small <- tsreg(u, rhs1)
  test <- waldtest(fit.large, fit.small)
  return((test[2,3]))                
}

ftest.stats2 <- sapply(dates2, f.calc2)
ftest.stats2

plot(ftest.stats2, type="l")

which.max(ftest.stats2)

dates2[84]

# Break 1: July 1965
# Break 2: January 1983

supF2 <- max(ftest.stats2)
supF2

# Can't compare supF2 to F statistic. Need special critical values for this.





# Bootstrapping 

inf <- dataset[,"inf"]

# Estimate model using our sample

fit <- tsreg(inf, u)
fit

# Generate new sample of inflation

res <- residuals(fit)

res.sim <- sample(res, replace = TRUE)

inf.sim <- 0.24 + 0.007*u + res.sim
fit.sim <- tsreg(inf.sim, u)
fit.sim

set.seed(100)
simcoef <- replicate(1000, {
  res.sim <- sample(res, replace = TRUE)
  inf.sim <- 0.24 + 0.007*u + res.sim
  fit.sim <- tsreg(inf.sim, u)
  coefficients(fit.sim)[1]^0.3/coefficients(fit.sim)[2]
})

# fit.sim[1] is alpha, fit.sim[2] is beta.

sd(simcoef)

plot(density(simcoef))
quantile(simcoef, probs = c(0.025, 0.975))

median(simcoef)

set.seed(100)
simcoef <- replicate(1000, {
  res.sim <- sample(res, replace = TRUE)
  inf.sim <- 0.24 + 0.007*u + res.sim
  fit.sim <- tsreg(inf.sim, u)
  coefficients(fit.sim)[2]
})


sd(simcoef)

plot(density(simcoef))
quantile(simcoef, probs = c(0.025, 0.975))



