rm(list=ls(all=TRUE))
library(urca)

Data <- read.csv('OTHER-COMPANIES-WEEKLY.csv')
Data <- ts(Data, start=c(1990,1), frequency=52) # Setting data as time series by month
SSTA <- Data[,"SSTA"]; Smucker <- na.omit(Data[,"SMUCKER"])
Stocks <- na.omit(Data[,c("ADM","CPB","CONAGRA","FMC","GMILLS","HSY","HORMEL","MKC","MOSAIC","SYSCO","TYSONS")])

################################### For Variables of Different start Dates #################

#### SSTA
SSTA_adf <-ur.df(SSTA, type=c("trend"), selectlags="AIC")
SSTA_adfgls <-ur.ers(SSTA, type="DF-GLS", model="trend")
SSTA_pp <-ur.pp(SSTA, type=("Z-tau"), model=c("trend"))
SSTA_kpss <-ur.kpss(SSTA, type=c("tau"))
SSTAstats <- cbind(SSTA_adf@teststat[1],SSTA_adfgls@teststat[1],SSTA_pp@teststat[1], SSTA_kpss@teststat[1])
SSTAcv <- cbind(SSTA_adf@cval[1,2],SSTA_adfgls@cval[2],SSTA_pp@cval[2], SSTA_kpss@cval[2])
SSTAtests <- rbind(SSTAstats,SSTAcv)

#### SSTA without trend
SSTA_adf <-ur.df(SSTA, type=c("none"), selectlags="AIC")
SSTA_adfgls <-ur.ers(SSTA, type="DF-GLS", model="constant")
SSTA_pp <-ur.pp(SSTA, type=("Z-tau"), model=c("constant"))
SSTA_kpss <-ur.kpss(SSTA, type = c("mu"))
SSTAstats <- cbind(SSTA_adf@teststat[1],SSTA_adfgls@teststat[1],SSTA_pp@teststat[1], SSTA_kpss@teststat[1])
SSTAcv <- cbind(SSTA_adf@cval[1,2],SSTA_adfgls@cval[2],SSTA_pp@cval[2], SSTA_kpss@cval[2])
SSTAtests <- rbind(SSTAstats,SSTAcv)


#### Smucker
smucker_adf <-ur.df(Smucker, type=c("trend"), selectlags="AIC")
smucker_adfgls <-ur.ers(Smucker, type="DF-GLS", model="trend")
smucker_pp <-ur.pp(Smucker, type=("Z-tau"), model=c("trend"))
smucker_kpss <-ur.kpss(Smucker, type=c("tau"))
smuckerstats <- cbind(smucker_adf@teststat[1],smucker_adfgls@teststat[1],smucker_pp@teststat[1], smucker_kpss@teststat[1])
smuckercv <- cbind(smucker_adf@cval[1,2],smucker_adfgls@cval[2],smucker_pp@cval[2], smucker_kpss@cval[2])
smuckertests <- rbind(smuckerstats,smuckercv)


#### Smucker without intercept and/or trend

smucker_adf <-ur.df(Smucker, type=c("none"), selectlags="AIC")
smucker_adfgls <-ur.ers(Smucker, type="DF-GLS", model="constant")
smucker_pp <-ur.pp(Smucker, type=("Z-tau"), model=c("constant"))
smucker_kpss <-ur.kpss(Smucker, type=c("mu"))
smuckerstats <- cbind(smucker_adf@teststat[1],smucker_adfgls@teststat[1],smucker_pp@teststat[1], smucker_kpss@teststat[1])
smuckercv <- cbind(smucker_adf@cval[1,2],smucker_adfgls@cval[2],smucker_pp@cval[2], smucker_kpss@cval[2])
smuckertests <- rbind(smuckerstats,smuckercv)

############################ ADF Tests For Levels ##############################

stockunitroottests = matrix(NA,nrow=ncol(Stocks),ncol=4)

for (ii in 1:ncol(Stocks))
{  
  testadf <-ur.df(Stocks[,ii], type=c("trend"), selectlags="AIC")
  testadfgls <-ur.ers(Stocks[,ii], type="DF-GLS", model="trend")
  testpp <-ur.pp(Stocks[,ii], type=("Z-tau"), model=c("trend"))
  testkpss <-ur.kpss(Stocks[,ii], type=c("tau"))
  stockunitroottests[ii,1] <- testadf@teststat[1] 
  stockunitroottests[ii,2] <- testadfgls@teststat[1] 
  stockunitroottests[ii,3] <- testpp@teststat[1] 
  stockunitroottests[ii,4] <- testkpss@teststat[1]
}

#### ADF tests without trend

stockunitroottests = matrix(NA,nrow=ncol(Stocks),ncol=4)

for (ii in 1:ncol(Stocks))
{  
  testadf <-ur.df(Stocks[,ii], type=c("none"), selectlags="AIC")
  testadfgls <-ur.ers(Stocks[,ii], type="DF-GLS", model="constant")
  testpp <-ur.pp(Stocks[,ii], type=("Z-tau"), model=c("constant"))
  testkpss <-ur.kpss(Stocks[,ii], type=c("mu"))
  stockunitroottests[ii,1] <- testadf@teststat[1] 
  stockunitroottests[ii,2] <- testadfgls@teststat[1] 
  stockunitroottests[ii,3] <- testpp@teststat[1] 
  stockunitroottests[ii,4] <- testkpss@teststat[1]
}


#####################################################################################################
###################################### ADF Tests For Growth #########################################
rm(list=ls(all=TRUE))
library(urca)

Data <- read.csv('OTHER-COMPANIES-WEEKLY.csv')
Data <- ts(Data, start=c(1990,1), frequency=52) # Setting data as time series by month
SSTA <- diff(Data[,"SSTA"]); Smucker <- diff(log(na.omit(Data[,"SMUCKER"])))
Stocks <- diff(log(na.omit(Data[,c("ADM","CPB","CONAGRA","FMC","GMILLS","HSY","HORMEL","MKC","MOSAIC","SYSCO","TYSONS")])))

################################### For Variables of Different start Dates ###########################
#### SSTA
SSTA_adf <-ur.df(SSTA, type=c("trend"), selectlags="AIC")
SSTA_adfgls <-ur.ers(SSTA, type="DF-GLS", model="trend")
SSTA_pp <-ur.pp(SSTA, type=("Z-tau"), model=c("trend"))
SSTA_kpss <-ur.kpss(SSTA, type=c("tau"))
SSTAstats <- cbind(SSTA_adf@teststat[1],SSTA_adfgls@teststat[1],SSTA_pp@teststat[1], SSTA_kpss@teststat[1])
SSTAcv <- cbind(SSTA_adf@cval[1,2],SSTA_adfgls@cval[2],SSTA_pp@cval[2], SSTA_kpss@cval[2])
SSTAtests <- rbind(SSTAstats,SSTAcv)


#### SSTA without trend
SSTA_adf <-ur.df(SSTA, type=c("none"), selectlags="AIC")
SSTA_adfgls <-ur.ers(SSTA, type="DF-GLS", model="constant")
SSTA_pp <-ur.pp(SSTA, type=("Z-tau"), model=c("constant"))
SSTA_kpss <-ur.kpss(SSTA, type = c("mu"))
SSTAstats <- cbind(SSTA_adf@teststat[1],SSTA_adfgls@teststat[1],SSTA_pp@teststat[1], SSTA_kpss@teststat[1])
SSTAcv <- cbind(SSTA_adf@cval[1,2],SSTA_adfgls@cval[2],SSTA_pp@cval[2], SSTA_kpss@cval[2])
SSTAtests <- rbind(SSTAstats,SSTAcv)


#### Smucker
smucker_adf <-ur.df(Smucker, type=c("trend"), selectlags="AIC")
smucker_adfgls <-ur.ers(Smucker, type="DF-GLS", model="trend")
smucker_pp <-ur.pp(Smucker, type=("Z-tau"), model=c("trend"))
smucker_kpss <-ur.kpss(Smucker, type=c("tau"))
smuckerstats <- cbind(smucker_adf@teststat[1],smucker_adfgls@teststat[1],smucker_pp@teststat[1], smucker_kpss@teststat[1])
smuckercv <- cbind(smucker_adf@cval[1,2],smucker_adfgls@cval[2],smucker_pp@cval[2], smucker_kpss@cval[2])
smuckertests <- rbind(smuckerstats,smuckercv)

#### Smucker without intercept and/or trend

smucker_adf <-ur.df(Smucker, type=c("none"), selectlags="AIC")
smucker_adfgls <-ur.ers(Smucker, type="DF-GLS", model="constant")
smucker_pp <-ur.pp(Smucker, type=("Z-tau"), model=c("constant"))
smucker_kpss <-ur.kpss(Smucker, type=c("mu"))
smuckerstats <- cbind(smucker_adf@teststat[1],smucker_adfgls@teststat[1],smucker_pp@teststat[1], smucker_kpss@teststat[1])
smuckercv <- cbind(smucker_adf@cval[1,2],smucker_adfgls@cval[2],smucker_pp@cval[2], smucker_kpss@cval[2])
smuckertests <- rbind(smuckerstats,smuckercv)


############################ ADF Tests For Levels ##############################

stockunitroottests = matrix(NA,nrow=ncol(Stocks),ncol=4)

for (ii in 1:ncol(Stocks))
{  
  testadf <-ur.df(Stocks[,ii], type=c("trend"), selectlags="AIC")
  testadfgls <-ur.ers(Stocks[,ii], type="DF-GLS", model="trend")
  testpp <-ur.pp(Stocks[,ii], type=("Z-tau"), model=c("trend"))
  testkpss <-ur.kpss(Stocks[,ii], type=c("tau"))
  stockunitroottests[ii,1] <- testadf@teststat[1] 
  stockunitroottests[ii,2] <- testadfgls@teststat[1] 
  stockunitroottests[ii,3] <- testpp@teststat[1] 
  stockunitroottests[ii,4] <- testkpss@teststat[1]
}

#### ADF tests without trend

stockunitroottests = matrix(NA,nrow=ncol(Stocks),ncol=4)

for (ii in 1:ncol(Stocks))
{  
  testadf <-ur.df(Stocks[,ii], type=c("none"), selectlags="AIC")
  testadfgls <-ur.ers(Stocks[,ii], type="DF-GLS", model="constant")
  testpp <-ur.pp(Stocks[,ii], type=("Z-tau"), model=c("constant"))
  testkpss <-ur.kpss(Stocks[,ii], type=c("mu"))
  stockunitroottests[ii,1] <- testadf@teststat[1] 
  stockunitroottests[ii,2] <- testadfgls@teststat[1] 
  stockunitroottests[ii,3] <- testpp@teststat[1] 
  stockunitroottests[ii,4] <- testkpss@teststat[1]
}

