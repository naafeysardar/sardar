df <- read.csv("Emerging.csv", header = TRUE, check.names=F, fileEncoding = 'UTF-8-BOM')

pca1 <- princomp(df, scores=TRUE, cor=TRUE)
summary(pca1)
loadings(pca1)
pca <- pca1$scores[,1]/10

##### Upload Packages #####

project.raw <- read.csv("Quarterly Employment Data.csv", header=TRUE)
data <- ts(project.raw, start=c(1991,1), frequency=4)

library(tstools)
library(vars)
library(lmtest)
library(sandwich)
library(car)

oildd <- data[,"oildd"]
oilad <- data[,"oilad"]
emerging.factor <- data[,"factor"]

apparel <- data[,"apparel"]
chemicals <- data[,"chemicals"]
computer <- data[,"computer"]
electrical <- data[,"electrical"]
fabricated <- data[,"fabricated"]
food <- data[,"food"]
furniture <- data[,"furniture"]
machinery <- data[,"machinery"]
misc.dur <- data[,"misc.durable"]
misc.nondur <- data[,"misc.nondurable"]
nonmetallic <- data[,"nonmetallic"]
paper <- data[,"paper"]
petroleum <- data[,"petroleum"]
plastic <- data[,"plastics"]
primary <- data[,"primary"]
printing <- data[,"printing"]
textile.mills <- data[,"textile.mills"]
textile.product <- data[,"textile.product"]
transport <- data[,"transport"]
wood <- data[,"wood"]
dnm <- data[,"dnm"]

varinfo <- list(
  list(apparel, "Apparel"),
  list(chemicals, "Chemicals"),
  list(computer, "Computer"),
  list(electrical, "Electrical"),
  list(fabricated, "Fabricated Metals"),
  list(food, "Food"),
  list(furniture, "Furniture"),
  list(machinery, "Machinery"),
  list(misc.dur, "Misc. Durable"),
  list(misc.nondur, "Misc. Nondurable"),
  list(nonmetallic, "Nonmetallic"),
  list(paper, "Paper"),
  list(petroleum, "Petroleum"),
  list(plastic, "Plastic"),
  list(primary, "Primary Metals"),
  list(printing, "Printing"),
  list(textile.mills, "Textile Mills"),
  list(textile.product, "Textile Products"),
  list(transport, "Transportation"),
  list(wood, "Wood"),
  list(dnm, "Non-Man"))

fit <- tsreg(oilad, lags(emerging.factor, 0:4))
oilad.hat <- ts(fitted(fit), start=c(1992,1), frequency = 4)

##### Lag Selection #####

lag.test <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  output <- lapply(1:8, function(k) {
    rhs <- ts.intersect(lags(oildd, 0:k), lags(oilad.hat, 0:k))
    fit <- tsreg(macro, rhs)
    return(list(fit=fit, lags=k, aic=AIC(fit), sic=BIC(fit)))
  })
  
  # Returns the lag length with the lowest criteria value
  lowest.crit <- function(fits, min.crit=Inf, min.lags=0, crit="aic") {
    if (length(fits) < 1) {
      return(min.lags)
    } else {
      current.crit <- if (fits[[1]][[crit]] < min.crit) { fits[[1]][[crit]] } else { min.crit }
      current.lags <- if (fits[[1]][[crit]] < min.crit) { fits[[1]]$lags } else { min.lags }
      return(lowest.crit(fits[-1], current.crit, current.lags))
    }
  }
  k <- lowest.crit(output)
  return(list(k, varname=varname))
}

lagtest.results <- lapply(varinfo, lag.test)
names(lagtest.results) <- sapply(lagtest.results, function(z) { z$varname })

##### Wald Test #####

cons <- ts(1, start=c(1994,1), end=c(2007,4), frequency = 4)

ols.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  rhs.u <- ts.combine(lags(oildd, 0:8), lags(oilad.hat, 0:8))
  fit.u <- tsreg(macro, rhs.u)
  
  fit.r <- tsreg(macro, cons)
  
  return(list(durbinWatsonTest(fit.u), waldtest(fit.u, fit.r, vcov = NeweyWest), varname=varname))
}

linear <- lapply(varinfo, ols.linear)
names(linear) <- sapply(linear, function(z) { z$varname })

##### Estimate the first-stage equation #####

ols.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  rhs <- ts.combine(lags(oildd, 0:8), lags(oilad.hat, 0:8))
  fit <- tsreg(macro, rhs)
  return(list(dlhat_2000=sum(fitted(fit)[1:24]), dlhat_2007=sum(fitted(fit)[25:56]), 
              dlhat = sum(fitted(fit)), varname=varname))
}

linear <- lapply(varinfo, ols.linear)
names(linear) <- sapply(linear, function(z) { z$varname })

ols.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  rhs <- ts.combine(lags(oildd, 0:8), lags(oilad.hat, 0:8))
  fit <- tsreg(macro, rhs)
  f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
  return(list(f, macro, varname=varname))
}

linear <- lapply(varinfo, ols.linear)
names(linear) <- sapply(linear, function(z) { z$varname })

###############################################################

rhs <- ts.combine(lags(oildd, 0:8), lags(oilad.hat, 0:8))
fit <- tsreg(apparel, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(chemicals, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(computer, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(electrical, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(fabricated, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(food, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(furniture, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(machinery, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(misc.dur, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(misc.nondur, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(nonmetallic, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(paper, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(petroleum, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(plastic, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(primary, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(printing, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(textile.mills, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(textile.product, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(transport, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(wood, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)

fit <- tsreg(nonman, rhs)
f <- (ts(fitted(fit), start=c(1994,1), frequency=4))
perc <- (100+f)/100
View(perc)