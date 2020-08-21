data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]

library(tstools)

fit.nochange <- tsreg(u, lags(u, 1:2))
fit.nochange

# If the sum of the coefficient on the AR terms less than 1,
# then non-stationery. 

# Structural break date is Jan 1980.

d <- time.dummy(u, c(1980, 1), "end")
d

rhs <- ts.intersect(lags(u, 1:2), d, d*lags(u, 1:2))

colnames <- c("u1","u2","d", "d_u1", "d_u2")
rhs

# Estimating the Structural Break Model
# All three parameters are changing

fit.change <- tsreg(u, rhs)
fit.change

# Wald, F Test for no structural change

library(lmtest)
lrtest(fit.change, fit.nochange)

# p-value = 0.4, don't reject H0.
# That means no structural break at Jan 1980.


# Next we estimate the restricted model.
# Then the unrestricted model. This is 
# followed by a Wald Test of both these
# models. 


# Restricted Model

rhs1 <- rhs[, 1:2]
fit.large <- tsreg(u, rhs)

# Unrestricted Model 

fit.small <- tsreg(u, rhs1)

waldtest(fit.large, fit.small)

waldtest(fit.large, fit.small, test="Chisq")

