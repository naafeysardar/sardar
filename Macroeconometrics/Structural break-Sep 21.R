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

# First, need to find set of break dates.

# 15% of the sample in each regime. Pre-break sample:

floor(0.15*length(u))

# Withnin our sample, the earliest break date would be after 127 observations.

# 15% of the sample in each regime. Post-break sample:

floor(0.85*length(u))

# Withnin our sample, the latest break date would be after 719 observations.

time(u)[127]
time(u)[719]

# All possible sturctural break dates are between observations 127 & 719.

possible.breaks <- dates(c(1958,7), c(2007,11), 12)
possible.breaks

f.calc <- function(break.date){
  d <- time.dummy(u, break.date, "end")
  rhs <- ts.intersect(lags(u, 1:2), d, d*lags(u, 1:2))
  rhs1 <- rhs[, 1:2]
  fit.large <- tsreg(u, rhs)
  fit.small <- tsreg(u, rhs1)
  test <- waldtest(fit.large, fit.small)
  return(test[2,3])                 # Return 3rd column of 2nd row
}

f.calc(2000.0)
f.calc(1980)

ftest.stats <- sapply(possible.breaks, f.calc)    # Running all break dates through the defined function
ftest.stats

plot(ftest.stats, type="l")

# We want higher F-statistic.

which.max(ftest.stats)
possible.breaks[85]

# July 1965 is estimated as the best structural break test. 

ftest.stats[85]


# Alternative Method of finding break date

library(strucchange)

u1 <- rhs[,1]
u2 <- rhs[,2]
 
model <- u ~ u1 + u2

dataset <- ts.intersect(u, u1, u2)

fs <- Fstats(model, from=c(1958,7), to=c(2007,11), data = dataset)

plot(fs, aveF=TRUE)

sctest(fs, type = "supF")
sctest(fs, type = "expF")
sctest(fs, type = "aveF")

sctest(model, type="expF", from = c(1958,7), to=c(2007,11), data = dataset)
