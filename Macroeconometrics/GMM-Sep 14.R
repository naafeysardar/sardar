# GMM Estimation


data.raw <  - read.csv("u-inf.csv", header=TRUE)
dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]
inf <- dataset[,"inf"]

plot(u)
plot(inf)

# Compute Moment Condition
# First Moment Condtion

dev1 <- function(a, b, u, inf){
  ds <- ts.intersect(u, inf, lag(inf, 1))
  return((ds[,1] - a - b*ds[,3])*ds[,2])
}

# (Dev1|alpha = 0,beta=0.1)

dev1(0.0, 0.1, u, inf)

# Mean(Dev1|alpha = 0,beta=0.1)

mean(dev1(0.0, 0.1, u, inf))

# Second Moment Condition

library(tstools)

dev2 <- function(a, b, u, inf){
  ds <- ts.intersect(u, lags(inf,1), lag(inf, 1))
  return((ds[,1] - a - b*ds[,3])*ds[,2])
}

dev2(0.0, 0.1, u, inf)
mean(dev2(0.0, 0.1, u, inf))

# Moment condition fails because g is not equal to 0.
# Define the objective function g=min[{mean(dev1)}^2+{mean(dev2)}^2]

objfun <- function(parameters) {
  a <- parameters[1]
  b <- parameters[2]
  return(mean(dev1(a, b, u, inf))^2 
         + mean(dev2(a, b, u, inf))^2)
  }

# Compute Moment Condition such that starting values are a=0.0, b=0.0

optim(c(0.0, 0.0), objfun)

# 10,000 Iterations

optim(c(0.0, 0.0), objfun, control=list(maxit=10000))
