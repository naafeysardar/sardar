rm(list=ls(all=TRUE))
library(vars)
library(nleqslv)
library(tstools)
library(msm)
library(functional)
library(systemfit)

source('functions.R')


Data <- read.csv('OTHER-COMPANIES-1990.csv')
Data <- ts(Data, start=c(1990,1), frequency=12) # Setting data as time series by month
SSTA <- Data[,"SSTA"]
Data <- diff(log(Data)) # Growth rate by taking first difference of log of variables
DATA <- cbind(SSTA,Data)
DATA <- DATA[,c("SSTA","Data.SOI","Data.SP","Data.CPI","Data.INDPRO","Data.PCOM","Data.ADM",
                "Data.CONAGRA","Data.CPB", "Data.GMILLS", "Data.HORMEL","Data.SMUCKER",
                "Data.MKC","Data.HSY","Data.TYSONS","Data.FMC","Data.MOSAIC", "Data.SYSCO")]

colnames(DATA) <- c("SSTA","SOI", "SP","CPI","INDPRO","PCOM","ADM","CONAGRA","CPB",
                    "GMILLS","HORMEL","SMUCKER","MKC","HSY","TYSONS","FMC","MOSAIC", "SYSCO")

ssta <- DATA[,"SSTA"]
ind <- DATA[,"INDPRO"]
cpi <- DATA[,"CPI"]
sp <- DATA[,"SP"]

pcom <- DATA[,"PCOM"]
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

##### Generate Dummy Variable #####

d <- (SSTA > 0.5)
inter <- d*SSTA

par(mfrow=c(2,3))

irf.asymmetry <- function(obj) {
  stock <- obj[[1]]
  varname <- obj[[2]]
  dataset.stock <- ts.intersect(cpi, ind, stock, ssta, inter)
  
  # Define Initial Response Vectors
  x0.nino <- c(0, 0.88, 0.88)
  x0.other <- c(0, 0.88, 0)
  
  rfvar <- lapply(1:12, function(h) {
    rhs <- lags(dataset.stock, h, type="bylag")
    coefficients(tsreg(stock, rhs))
  })
  
  # Impulse Response in Other
  response.other <- lapply(rfvar, function(b) {
    sum(b[4:6] * x0.other)
  })
  # Impulse Response in Nino
  response.nino <- lapply(rfvar, function(b) {
    sum(b[4:6] * x0.nino)
  })
  
  stock.nino <- unlist(cumsum(response.nino))
  stock.other <- unlist(cumsum(response.other))
  
  # Difference in response 
  Map(function(irf1, irf2) { irf1-irf2 }, stock.nino, stock.other)
  
  x <- embed(cbind(cpi, ind, stock, ssta, inter), 13)
  
  x <- as.data.frame(x)
  
  one.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 6:10]
  two.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 11:15]
  three.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 16:20] 
  four.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 21:25] 
  five.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 26:30] 
  six.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 31:35] 
  seven.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 36:40]
  eight.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 41:45]
  nine.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 46:50] 
  ten.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 51:55] 
  eleven.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 56:60] 
  twelve.step <- embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 3] ~ embed(cbind(cpi, ind, stock, ssta, inter), 13)[, 61:65] 
  
  system <- list(one.step, two.step, three.step, four.step, five.step, six.step,
                 seven.step, eight.step, nine.step, ten.step, eleven.step, twelve.step)
  
  fitsur <- systemfit(system, data = x)
  
  v <- matrix(vcov(fitsur), nrow = 72, ncol = 72)
  
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
  
  d.stock <- stock.nino - stock.other
  
  # 1-Step #
  theta <- coefficients(fitsur)[c(6)]
  v1 <- as.matrix(x.newey[c(6), c(6)])
  std.err.1 <- deltamethod(~(x1), theta, v1)
  # 2-Step #
  theta <- coefficients(fitsur)[c(6, 12)]
  v2 <- as.matrix(x.newey[c(6, 12), c(6, 12)])
  std.err.2 <- deltamethod(~(x1+x2), theta, v2)
  # 3-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18)]
  v3 <- as.matrix(x.newey[c(6, 12, 18), c(6, 12, 18)])
  std.err.3 <- deltamethod(~(x1+x2+x3), theta, v3)
  # 4-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24)]
  v4 <- as.matrix(x.newey[c(6, 12, 18, 24), c(6, 12, 18, 24)])
  std.err.4 <- deltamethod(~(x1+x2+x3+x4), theta, v4)
  # 5-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30)]
  v5 <- as.matrix(x.newey[c(6, 12, 18, 24, 30), c(6, 12, 18, 24, 30)])
  std.err.5 <- deltamethod(~(x1+x2+x3+x4+x5), theta, v5)
  # 6-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30, 36)]
  v6 <- as.matrix(x.newey[c(6, 12, 18, 24, 30, 36), c(6, 12, 18, 24, 30, 36)])
  std.err.6 <- deltamethod(~(x1+x2+x3+x4+x5+x6), theta, v6)
  # 7-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30, 36, 42)]
  v7 <- as.matrix(x.newey[c(6, 12, 18, 24, 30, 36, 42), 
                          c(6, 12, 18, 24, 30, 36, 42)])
  std.err.7 <- deltamethod(~(x1+x2+x3+x4+x5+x6+x7), theta, v7)
  # 8-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30, 36, 42, 48)]
  v8 <- as.matrix(x.newey[c(6, 12, 18, 24, 30, 36, 42, 48), 
                          c(6, 12, 18, 24, 30, 36, 42, 48)])
  std.err.8 <- deltamethod(~(x1+x2+x3+x4+x5+x6+x7+x8), theta, v8)
  # 9-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30, 36, 42, 48, 54)]
  v9 <- as.matrix(x.newey[c(6, 12, 18, 24, 30, 36, 42, 48, 54), 
                          c(6, 12, 18, 24, 30, 36, 42, 48, 54)])
  std.err.9 <- deltamethod(~(x1+x2+x3+x4+x5+x6+x7+x8+x9), theta, v9)
  # 10-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60)]
  v10 <- as.matrix(x.newey[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60), 
                           c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60)])
  std.err.10 <- deltamethod(~(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10), theta, v10)
  # 11-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66)]
  v11 <- as.matrix(x.newey[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66), 
                           c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66)])
  std.err.11 <- deltamethod(~(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11), theta, v11)
  # 12-Step #
  theta <- coefficients(fitsur)[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72)]
  v12 <- as.matrix(x.newey[c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72), 
                           c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72)])
  std.err.12 <- deltamethod(~(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12), theta, v12)
  
  s <- unlist(list(std.err.1, std.err.2, std.err.3, std.err.4, std.err.5, std.err.6,
                   std.err.7, std.err.8, std.err.9, std.err.10, std.err.11, std.err.12))
  d.stock.upper <- d.stock + 0.88*1.64*s
  d.stock.lower <- d.stock - 0.88*1.64*s
  
  plot(d.stock, type = "l", col = "black", xlab = "Months Afer Shock", ylab = "Percent", ylim= c(-0.45, 0.45), main = varname); abline(h=0)
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
