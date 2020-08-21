# Bootstrapping 

data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]

inf <- dataset[,"inf"]

# Estimate model using our sample

fit <- tsreg(inf, u)
summary(fit)

coefficients(summary(fit))

# Generate new sample of inflation

res <- residuals(fit)

# Bootstrapping Standard Errors

set.seed(100)
simcoef <- replicate(1000, {
  res.sim <- sample(res, replace = TRUE)
  inf.sim <- 0.24 + 0.007*u + res.sim
  fit.sim <- tsreg(inf.sim, u)
  coefficients(fit.sim)[2]
})


sd(simcoef)

# Compare the outcome of Standard Errors 
# through OLS and bootstrapping method

# IID Bootstrapping Test Statistic

set.seed(100)
tstats <- replicate(1000, {
  res.sim <- sample(res, replace = TRUE)
  inf.sim <- 0.24 + 0.007*u + res.sim
  fit.sim <- tsreg(inf.sim, u)
  coefficients(summary(fit.sim))[2,3]
})

quantile(tstats, probs=c(0.025, 0.975))     # We want CI to be [-1.96, 1.96]

# Recentered t-stat

set.seed(100)
tstats <- replicate(1000, {
  res.sim <- sample(res, replace = TRUE)
  inf.sim <- 0.24 + 0.007*u + res.sim
  fit.sim <- tsreg(inf.sim, u)
  coefficients(summary(fit.sim))[2,3] - 1.05  #1.05 is the t-stat from OLS
})

# Critical Values for recentered t-stat

quantile(tstats, probs=c(0.025, 0.975))

# We compare the original value of t-stat to the recentered
# t-statistic CI. We don't reject the null in this case because
# our original value lies within the CI. 

# Pairs Bootstrap

res <- residuals(fit)

set.seed(100)
simcoef <- replicate(1000, {
  obs <- sample(1:length(u), replace = TRUE)     # Resampling data
  inf.new <- ts(inf[obs])
  u.new <- ts(u[obs])
  fit.new <- tsreg(inf.new, u.new)
  coefficients(fit.new)[2]
})

# SE calculated through the Pairs Bootstrap

sd(simcoef)

# Confidence Interval of Beta calculated through Pairs

quantile(simcoef, probs=c(0.025, 0.975))

