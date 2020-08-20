rm(list=ls(all=TRUE))
library(vars)
library(nleqslv)
library(tstools)
library(msm)
library(functional)
library(systemfit)

source('functions.R')


Data <- read.csv('OTHER-COMPANIES-WEEKLY.csv')
Data <- ts(Data, start=c(1990,1), frequency=52) # Setting data as time series by WEEK
SSTA <- Data[,"SSTA"]
Data <- diff(log(Data)) # Growth rate by taking first difference of log of variables
DATA <- cbind(SSTA,Data)
DATA <- DATA[,c("SSTA","Data.ADM","Data.CONAGRA","Data.CPB", "Data.GMILLS", "Data.HORMEL","Data.SMUCKER",
                "Data.MKC","Data.HSY","Data.TYSONS","Data.FMC","Data.MOSAIC", "Data.SYSCO",
                "Data.ALICO", "Data.INGREDION", "Data.ANDERSONS", "Data.BUNGE", "Data.ALTRIA")]

colnames(DATA) <- c("SSTA","ADM","CONAGRA","CPB", "GMILLS","HORMEL","SMUCKER","MKC","HSY","TYSONS",
                    "FMC","MOSAIC", "SYSCO", "ALICO", "INGREDION", "ANDERSONS", "BUNGE", "ALTRIA")

ssta <- DATA[,"SSTA"]

adm <- DATA[,"ADM"]
conagra <- DATA[,"CONAGRA"]
cpb <- DATA[,"CPB"]
gmills <- DATA[,"GMILLS"]
hormel <- DATA[,"HORMEL"]
smucker <- DATA[,"SMUCKER"]
mkc <- DATA[,"MKC"]
hsy <- DATA[,"HSY"]
tysons <- DATA[,"TYSONS"]
fmc <- DATA[,"FMC"]
mosaic <- DATA[,"MOSAIC"]
sysco <- DATA[,"SYSCO"]
alico <- DATA[,"ALICO"]
ingredion <- DATA[,"INGREDION"]
andersons <- DATA[,"ANDERSONS"]
bunge <- DATA[,"BUNGE"]
altria <- DATA[,"ALTRIA"]

##### Generate Dummy Variable #####

d <- (SSTA < -0.5)
inter <- d*SSTA

par(mfrow=c(2,3))

irf.asymmetry <- function(obj) {
  stock <- obj[[1]]
  varname <- obj[[2]]
  dataset.stock <- ts.intersect(stock, ssta, inter)
  
  # Define Initial Response Vectors
  x0.nino <- c(0, 0.92, 0)
  x0.nina <- c(0, 0.92, 0.92)
  
  rfvar <- lapply(1:12, function(h) {
    rhs <- lags(dataset.stock, h, type="bylag")
    coefficients(tsreg(stock, rhs))
  })
  
  # Impulse Response in nina
  response.nina <- lapply(rfvar, function(b) {
    sum(b[2:4] * x0.nina)
  })
  # Impulse Response in Nino
  response.nino <- lapply(rfvar, function(b) {
    sum(b[2:4] * x0.nino)
  })
  
  stock.nino <- unlist(cumsum(response.nino))
  stock.nina <- unlist(cumsum(response.nina))
  
  # Difference in response 
  Map(function(irf1, irf2) { irf1-irf2 }, stock.nino, stock.nina)
  
  x <- embed(cbind(stock, ssta, inter), 13)
  
  x <- as.data.frame(x)
  
  one.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 2:4]
  two.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 7:9]
  three.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 10:12] 
  four.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 13:15] 
  five.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 16:18] 
  six.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 19:21] 
  seven.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 22:24]
  eight.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 25:27]
  nine.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 28:30] 
  ten.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 31:33] 
  eleven.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 34:36] 
  twelve.step <- embed(cbind(stock, ssta, inter), 13)[, 1] ~ embed(cbind(stock, ssta, inter), 13)[, 37:39] 
  
  system <- list(one.step, two.step, three.step, four.step, five.step, six.step,
                 seven.step, eight.step, nine.step, ten.step, eleven.step, twelve.step)
  
  fitsur <- systemfit(system, data = x)
  
  v <- matrix(vcov(fitsur), nrow = 48, ncol = 48)
  
  rfvar <- lapply(1:12, function(h) {
    rhs <- lags(dataset.stock, h, type="bylag")
    tsreg(stock, rhs)
  })
  
  cov.mats <- lapply(rfvar, function(z) {
    NeweyWest(z)
  })
  
  
  replaceDiag <- function(a, mats, start=1) {
    if (length(mats) == 0) {
      return(a)
    } else {
      a[start:(start+nrow(mats[[1]])-1), start:(start+nrow(mats[[1]])-1)] <- 
        mats[[1]]
      Recall(a, mats[-1], start+nrow(mats[[1]]))
    }
  }
  
  x.newey <- replaceDiag(v, cov.mats)

  
  # Asymmetry
  
  d.stock <- stock.nino - stock.nina
  
  # 1-Step #
  theta <- coefficients(fitsur)[c(4)]
  v1 <- as.matrix(x.newey[c(4), c(4)])
  std.err.1 <- deltamethod(~(-x1), theta, v1)
  # 2-Step #
  theta <- coefficients(fitsur)[c(4, 8)]
  v2 <- as.matrix(x.newey[c(4, 8), c(4, 8)])
  std.err.2 <- deltamethod(~(-(x1+x2)), theta, v2)
  # 3-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12)]
  v3 <- as.matrix(x.newey[c(4, 8, 12), c(4, 8, 12)])
  std.err.3 <- deltamethod(~(-(x1+x2+x3)), theta, v3)
  # 4-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16)]
  v4 <- as.matrix(x.newey[c(4, 8, 12, 16), c(4, 8, 12, 16)])
  std.err.4 <- deltamethod(~(-(x1+x2+x3+x4)), theta, v4)
  # 5-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20)]
  v5 <- as.matrix(x.newey[c(4, 8, 12, 16, 20), c(4, 8, 12, 16, 20)])
  std.err.5 <- deltamethod(~(-(x1+x2+x3+x4+x5)), theta, v5)
  # 6-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20, 24)]
  v6 <- as.matrix(x.newey[c(4, 8, 12, 16, 20, 24), c(4, 8, 12, 16, 20, 24)])
  std.err.6 <- deltamethod(~(-(x1+x2+x3+x4+x5+x6)), theta, v6)
  # 7-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20, 24, 28)]
  v7 <- as.matrix(x.newey[c(4, 8, 12, 16, 20, 24, 28), 
                          c(4, 8, 12, 16, 20, 24, 28)])
  std.err.7 <- deltamethod(~(-(x1+x2+x3+x4+x5+x6+x7)), theta, v7)
  # 8-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20, 24, 28, 32)]
  v8 <- as.matrix(x.newey[c(4, 8, 12, 16, 20, 24, 28, 32), 
                          c(4, 8, 12, 16, 20, 24, 28, 32)])
  std.err.8 <- deltamethod(~(-(x1+x2+x3+x4+x5+x6+x7+x8)), theta, v8)
  # 9-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20, 24, 28, 32, 36)]
  v9 <- as.matrix(x.newey[c(4, 8, 12, 16, 20, 24, 28, 32, 36), 
                          c(4, 8, 12, 16, 20, 24, 28, 32, 36)])
  std.err.9 <- deltamethod(~(-(x1+x2+x3+x4+x5+x6+x7+x8+x9)), theta, v9)
  # 10-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40)]
  v10 <- as.matrix(x.newey[c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40), 
                           c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40)])
  std.err.10 <- deltamethod(~(-(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)), theta, v10)
  # 11-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44)]
  v11 <- as.matrix(x.newey[c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44), 
                           c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44)])
  std.err.11 <- deltamethod(~(-(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11)), theta, v11)
  # 12-Step #
  theta <- coefficients(fitsur)[c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48)]
  v12 <- as.matrix(x.newey[c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48), 
                           c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48)])
  std.err.12 <- deltamethod(~(-(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12)), theta, v12)
  
  s <- unlist(list(std.err.1, std.err.2, std.err.3, std.err.4, std.err.5, std.err.6,
                   std.err.7, std.err.8, std.err.9, std.err.10, std.err.11, std.err.12))
  d.stock.upper <- d.stock + 0.92*1.64*s
  d.stock.lower <- d.stock - 0.92*1.64*s
  
  plot(d.stock, type = "l", col = "black", xlab = "Weeks Afer Shock", ylab = "Percent", ylim= c(-0.1, 0.1), main = varname); abline(h=0)
  points(d.stock.upper , col = "black", type="l", lty=2)
  points(d.stock.lower, col = "black", type="l", lty=2)
  
}

varinfo <- list(
  list(adm, "ADM"),
  list(conagra, "Conagra"),
  list(cpb, "CPB"),
  list(gmills, "G Mills"),
  list(hormel, "Hormel"),
  list(smucker, "Smucker"),
  list(mkc, "MKC"),
  list(hsy, "HSY"),
  list(tysons, "Tysons"),
  list(fmc, "FMC"),
  list(mosaic, "Mosaic"),
  list(sysco, "Sysco")
)
irfs <- lapply(varinfo, irf.asymmetry)
