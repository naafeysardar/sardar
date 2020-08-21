library(tstools)
library(FKF)

gdp <- readRDS("gdp.RDS")
dgdp <- pctChange(gdp)

plot(dgdp)

# Put data in rows 

data <- matrix(dgdp, nrow=1)

dt <- matrix(0)
ct <- matrix(0)
Zt <- matrix(1)
Tt <- matrix(1)

a0 <- first(dgdp)
P0 <- matrix(100)

# Optim minimizes the function so we write negative
# sign so it is maximized

# Maximize the likelihood function

objfnc <- function(par) {
  -fkf(HHt=matrix(par[1]), GGt=matrix(par[2]), 
     yt=data, a0=a0, P0=P0,  
     dt=dt, ct=ct, Zt=Zt, Tt=Tt)$logLik
}

fit.fkf <- optim(c(0.5*var(dgdp), 0.5*var(dgdp)), objfnc)

# Find variances of both equations

fit.fkf$par


fkf.final <- fkf(a0, P0, dt, ct, Tt, Zt, HHt=matrix(fit.fkf$par[1]), GGt=matrix(fit.fkf$par[2]), yt=data)

names(attributes(fkf.final))
names(fkf.final)

# att is the first state variable at each point in time 

mean.tvp <- ts(fkf.final$att[1,], start=start(dgdp), frequency = 4)
plot(mean.tvp, main="Filtered Estiamte of the Mean 
     Real GDP Growth Rate")

