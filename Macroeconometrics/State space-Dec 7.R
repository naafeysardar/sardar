library(tstools)
library(mvtnorm)

gdp <- readRDS("gdp.RDS")
dgdp <- pctChange(gdp)

fit <- tsreg(dgdp, lags(dgdp, 1))
s2 <- summary(fit)$sigma^2

nu <- length(dgdp)

# Omitting NA Values

xx <- na.omit(cbind(1, lags(dgdp, 1)))
xx.inv <- solve(crossprod(xx))

theta <- coefficients(fit)

draws <- replicate(1000, {
  tau.draw <- rgamma(1, nu/2, nu*s2/2)      # Drawing one value
  rmvnorm(1, theta, (1/tau.draw)*xx.inv)
}, simplify = "matrix")

dim(draws)

# 2 rows & 1000 columns
# Rows are values of alpha and beta taken over 1000 replications 

rowMeans(draws)

theta

sd(draws[1,])
sd(draws[2,])

summary(fit)

0.7*sd(draws[2,])
theta[2]+0.039
theta[2]-0.039

plot(density(draws[2,]))
