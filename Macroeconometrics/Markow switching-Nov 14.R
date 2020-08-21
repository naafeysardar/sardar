gdp <- readRDS("gdp.RDS")
plot(gdp)

oil <- readRDS("oilppi.RDS")
plot(oil)

plot(window(oil, start=c(1970,1), end=c(1982,4)))

library(tstools)
dgdp <- pctChange(gdp)
doil.raw <- pctChange(oil)
doil <- toQuarterly(doil.raw, sum)

# Default option is the mean

plot(doil)

# Define threshold

d <- (doil < 0)

# Estimating whether increases and decreases in the price of oil have differnet effects

inter <- lags(d, 1:4) * lags(doil, 1:4)

colnames(inter) <- c("d1", "d2", "d3", "d4")

rhs.nonlinear <- cbind(lags(dgdp, 1:4), lags(d, 1), lags(doil, 1:4), inter)
fit.nonlinear <- tsreg(dgdp, rhs.nonlinear)
fit.nonlinear

# Markov Switching Model

library(MSwM)
data(example)
mod <- lm(y ~ x, example)
summary(mod)

mod.mswm <- msmFit(object = mod, k=2, sw=c(TRUE, TRUE, TRUE, TRUE), p=1)
summary(mod.mswm)