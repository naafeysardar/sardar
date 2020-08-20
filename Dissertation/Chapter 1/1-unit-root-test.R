rm(list=ls(all=TRUE))
library(urca)

Data <- read.csv('full-dataset.csv')
Data <- ts(Data, start=c(1972,4), frequency=4) 
price <- Data[,"price"]
consumption <- (na.omit(Data[,c("n.pce","n.durable","n.motor","n.furnishing","n.otherdurable","n.nondurable","n.food","n.clothing","n.gasoline","n.othernondurable","n.services")])/price)

############################ ADF Tests For Levels ##############################

stockunitroottests = matrix(NA,nrow=ncol(consumption),ncol=4)

for (ii in 1:ncol(consumption))
{  
  testadf <-ur.df(consumption[,ii], type=c("trend"), selectlags="AIC")
  testadfgls <-ur.ers(consumption[,ii], type="DF-GLS", model="trend")
  testpp <-ur.pp(consumption[,ii], type=("Z-tau"), model=c("trend"))
  testkpss <-ur.kpss(consumption[,ii], type=c("tau"))
  stockunitroottests[ii,1] <- testadf@teststat[1] 
  stockunitroottests[ii,2] <- testadfgls@teststat[1] 
  stockunitroottests[ii,3] <- testpp@teststat[1] 
  stockunitroottests[ii,4] <- testkpss@teststat[1]
}


#####################################################################################################
###################################### ADF Tests For Growth #########################################
rm(list=ls(all=TRUE))
library(urca)

Data <- read.csv('full-dataset.csv')
Data <- ts(Data, start=c(1972,4), frequency=4) 

price <- Data[,"price"]
consumption <- 100*diff(log(na.omit(Data[,c("n.pce","n.durable","n.motor","n.furnishing","n.otherdurable","n.nondurable","n.food","n.clothing","n.gasoline","n.othernondurable","n.services")])/price))

############################ ADF Tests For Levels ##############################

stockunitroottests = matrix(NA,nrow=ncol(consumption),ncol=4)

for (ii in 1:ncol(consumption))
{  
  testadf <-ur.df(consumption[,ii], type=c("trend"), selectlags="AIC")
  testadfgls <-ur.ers(consumption[,ii], type="DF-GLS", model="trend")
  testpp <-ur.pp(consumption[,ii], type=("Z-tau"), model=c("trend"))
  testkpss <-ur.kpss(consumption[,ii], type=c("tau"))
  stockunitroottests[ii,1] <- testadf@teststat[1] 
  stockunitroottests[ii,2] <- testadfgls@teststat[1] 
  stockunitroottests[ii,3] <- testpp@teststat[1] 
  stockunitroottests[ii,4] <- testkpss@teststat[1]
}

