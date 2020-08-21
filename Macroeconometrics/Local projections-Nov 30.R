gdp <- readRDS("gdp.RDS")
oil <- readRDS("oilppi.RDS")

library(tstools)
dgdp <- pctChange(gdp)
doil.raw <- pctChange(oil)
doil <- toQuarterly(doil.raw, sum)

# Define threshold

d <- (doil < 0)
inter <- lags(d, 1) * lags(doil, 1)

# Local Projections

# Impact Effect (initial response equation)

rhs <- ts.combine(lags(doil, 0:1), doil^2, doil^3, lags(dgdp, 1))

fit0 <- tsreg(dgdp, rhs)

mean(dgdp)

# GDP forecast when oil shock

sum(coefficients(fit0) * c (1, 0.1, 0.1, 0.01, 0.001, 0.008))

rhs1 <- ts.combine(lags(doil, 1), lags(doil^2, 1), lags(doil^3, 1), lags(dgdp, 1), lags(dgdp^2, 1), lags(dgdp^3, 1))
fit1 <- tsreg(dgdp, rhs1)

# GDP forecast 1 period after the oil shock takes place

sum(coefficients(fit1)*c(1, 0.1, 0.01, 0.001, 0.0076, 0.0076^2, 0.0076^3))

# GDP forecast 1 period after no oil shock 

sum(coefficients(fit1)*c(1, 0.0, 0.0, 0.0, 0.0079, 0.0079^2, 0.0079^3))


# Make foreacasts for 2 quarters in to the future

rhs2 <- ts.combine(lags(doil, 2), lags(doil^2, 2), lags(doil^3, 2), lags(dgdp, 2), lags(dgdp^2, 2), lags(dgdp^3, 2))
fit2 <- tsreg(dgdp, rhs2)

sum(coefficients(fit2)*c(1, 0.1, 0.01, 0.001, 0.0076, 0.0076^2, 0.0076^3))
sum(coefficients(fit2)*c(1, 0.0, 0.0, 0.0, 0.0079, 0.0079^2, 0.0079^3))


