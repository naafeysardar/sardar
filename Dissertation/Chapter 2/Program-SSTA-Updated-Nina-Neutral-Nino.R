rm(list=ls(all=TRUE))
library(vars)
library(nleqslv)
library(functional)
source('functions.R')
library(tsDyn)

################################# Responses to El Nino Shocks ####################################

Data <- read.csv('Data-Monthly-Updated.csv')
Data <- Data[Data[,"SSTA"]>= 0.5,] ## Observations for El Nino
Data <- ts(Data, start=c(1990,1), frequency=12) # Setting data as time series by month
SSTA <- Data[,"SSTA"]
Data <- 100*diff(log(Data)) # Growth rate by taking first difference of log of variables
DATA <- cbind(SSTA,Data)
DATA <- DATA[,c("SSTA","Data.SP","Data.Energy","Data.CPI","Data.INDPRO","Data.Alico","Data.Conagra",
                "Data.Tyson","Data.Monsanto", "Data.Ingredion", "Data.Hormel","Data.Andersons",
                "Data.Agrium","Data.Potash","Data.Mosaic","Data.Altria","Data.ArcherDaniels")]

colnames(DATA) <- c("SSTA","SP","Energy","CPI","INDPRO","Alico","Conagra","Tyson","Monsanto",
                    "Ingredion","Hormel","Andersons","Agrium","Potash","Mosaic","Altria","ArcherDaniels")

#######################################################################################################

par(mfrow=c(2,2))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","SP")])
fit <- VAR(vardata, 6)
fevd1 <- fevd(fit,n.ahead=36)
HD1 <- VARhd(Estimation=fit)
res_ssta <- irf(fit,impulse=c("SSTA"), response="SSTA", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
res_sp <- irf(fit,impulse=c("SSTA"), response="SP", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
res_cpi <- irf(fit,impulse=c("SSTA"), response="CPI", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
res_gdp <- irf(fit,impulse=c("SSTA"), response="INDPRO", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)

ts.plot(cbind(res_ssta$Lower$SSTA,res_ssta$irf$SSTA,res_ssta$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of SSTA",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

ts.plot(cbind(res_sp$Lower$SSTA,res_sp$irf$SSTA,res_sp$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of S&P Index",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

ts.plot(cbind(res_cpi$Lower$SSTA,res_cpi$irf$SSTA,res_cpi$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of CPI Inflation",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

ts.plot(cbind(res_gdp$Lower$SSTA,res_gdp$irf$SSTA,res_gdp$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Output",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

#######################################################################################################

par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Alico")])
fit <- VAR(vardata, 6)
fevd1 <- fevd(fit,n.ahead=36)
HD1 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Alico", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Alico",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Conagra")])
fit <- VAR(vardata, 6)
fevd2 <- fevd(fit,n.ahead=36)
HD2 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Conagra", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Conagra",col=c("red","black", "red"),lwd=c(2,2,2),lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Tyson")])
fit <- VAR(vardata, 6)
fevd3 <- fevd(fit,n.ahead=36)
HD3 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Tyson", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Tyson",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Monsanto")])
fit <- VAR(vardata, 6)
fevd4 <- fevd(fit,n.ahead=36)
HD4 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Monsanto", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Monsanto",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Ingredion")])
fit <- VAR(vardata, 6)
fevd5 <- fevd(fit,n.ahead=36)
HD5 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Ingredion", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Ingredion",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Hormel")])
fit <- VAR(vardata, 6)
fevd6 <- fevd(fit,n.ahead=36)
HD6 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Hormel", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Hormel",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)


par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Andersons")])
fit <- VAR(vardata, 6)
fevd7 <- fevd(fit,n.ahead=36)
HD7 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Andersons", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Andersons",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Agrium")])
fit <- VAR(vardata, 6)
fevd8 <- fevd(fit,n.ahead=36)
HD8 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Agrium", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Agrium",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Potash")])
fit <- VAR(vardata, 6)
fevd9 <- fevd(fit,n.ahead=36)
HD9 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Potash", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Potash",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Mosaic")])
fit <- VAR(vardata, 6)
fevd10 <- fevd(fit,n.ahead=36)
HD10 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Mosaic", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Mosaic",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Altria")])
fit <- VAR(vardata, 6)
fevd11 <- fevd(fit,n.ahead=36)
HD11 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Altria", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Altria",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","ArcherDaniels")])
fit <- VAR(vardata, 6)
fevd12 <- fevd(fit,n.ahead=36)
HD12 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="ArcherDaniels", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Archer Daniels",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

allfevd <- cbind(fevd1$SSTA[,"Alico"],fevd2$SSTA[,"Conagra"],fevd3$SSTA[,"Tyson"],fevd4$SSTA[,"Monsanto"],
                 fevd5$SSTA[,"Ingredion"],fevd6$SSTA[,"Hormel"],fevd7$SSTA[,"Andersons"],fevd8$SSTA[,"Agrium"],
                 fevd9$SSTA[,"Potash"],fevd10$SSTA[,"Mosaic"],fevd11$SSTA[,"Altria"],fevd12$SSTA[,"ArcherDaniels"])

write.csv(allfevd, file="fevd_Nino.csv")



#################################################################################################
############################ Now we will examine La Nina Responses ##############################
#################################################################################################

rm(list=ls(all=TRUE))
library(vars)
library(nleqslv)
library(functional)
source('functions.R')


################################# Responses to La Nino Shocks ####################################

Data <- read.csv('Data-Monthly-Updated.csv')
Data <- Data[Data[,"SSTA"]<= -0.5,] ## Observations for El Nino
Data <- ts(Data, start=c(1990,1), frequency=12) # Setting data as time series by month
SSTA <- Data[,"SSTA"]
Data <- 100*diff(log(Data)) # Growth rate by taking first difference of log of variables
DATA <- cbind(SSTA,Data)
DATA <- DATA[,c("SSTA","Data.SP","Data.Energy","Data.CPI","Data.INDPRO","Data.Alico","Data.Conagra",
                "Data.Tyson","Data.Monsanto", "Data.Ingredion", "Data.Hormel","Data.Andersons",
                "Data.Agrium","Data.Potash","Data.Mosaic","Data.Altria","Data.ArcherDaniels")]

colnames(DATA) <- c("SSTA","SP","Energy","CPI","INDPRO","Alico","Conagra","Tyson","Monsanto",
                    "Ingredion","Hormel","Andersons","Agrium","Potash","Mosaic","Altria","ArcherDaniels")


#######################################################################################################

#######################################################################################################

par(mfrow=c(2,2))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","SP")])
fit <- VAR(vardata, 6)
fevd1 <- fevd(fit,n.ahead=36)
HD1 <- VARhd(Estimation=fit)
res_ssta <- irf(fit,impulse=c("SSTA"), response="SSTA", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)
res_sp <- irf(fit,impulse=c("SSTA"), response="SP", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)
res_cpi <- irf(fit,impulse=c("SSTA"), response="CPI", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)
res_gdp <- irf(fit,impulse=c("SSTA"), response="INDPRO", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)

ts.plot(cbind(res_ssta$Lower$SSTA,res_ssta$irf$SSTA,res_ssta$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of SSTA",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

ts.plot(cbind(res_sp$Lower$SSTA,res_sp$irf$SSTA,res_sp$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of S&P Index",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-12,12));abline(h=0)

ts.plot(cbind(res_cpi$Lower$SSTA,res_cpi$irf$SSTA,res_cpi$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of CPI Inflation",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

ts.plot(cbind(res_gdp$Lower$SSTA,res_gdp$irf$SSTA,res_gdp$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Output",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)


par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Alico")])
fit <- VAR(vardata, 6)
fevd1 <- fevd(fit,n.ahead=36)
HD1 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Alico", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Alico",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Conagra")])
fit <- VAR(vardata, 6)
fevd2 <- fevd(fit,n.ahead=36)
HD2 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Conagra", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Conagra",col=c("red","black", "red"),lwd=c(2,2,2),lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Tyson")])
fit <- VAR(vardata, 6)
fevd3 <- fevd(fit,n.ahead=36)
HD3 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Tyson", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Tyson",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Monsanto")])
fit <- VAR(vardata, 6)
fevd4 <- fevd(fit,n.ahead=36)
HD4 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Monsanto", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA/2.3,res_alico$Upper$SSTA/2),xlab="Months", ylab="percent",
        main="Response of Monsanto",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Ingredion")])
fit <- VAR(vardata, 6)
fevd5 <- fevd(fit,n.ahead=36)
HD5 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Ingredion", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Ingredion",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Hormel")])
fit <- VAR(vardata, 6)
fevd6 <- fevd(fit,n.ahead=36)
HD6 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Hormel", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Hormel",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)


par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Andersons")])
fit <- VAR(vardata, 6)
fevd7 <- fevd(fit,n.ahead=36)
HD7 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Andersons", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Andersons",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Agrium")])
fit <- VAR(vardata, 6)
fevd8 <- fevd(fit,n.ahead=36)
HD8 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Agrium", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Agrium",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Potash")])
fit <- VAR(vardata, 6)
fevd9 <- fevd(fit,n.ahead=36)
HD9 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Potash", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Potash",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Mosaic")])
fit <- VAR(vardata, 6)
fevd10 <- fevd(fit,n.ahead=36)
HD10 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Mosaic", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Mosaic",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Altria")])
fit <- VAR(vardata, 6)
fevd11 <- fevd(fit,n.ahead=36)
HD11 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Altria", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Altria",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","ArcherDaniels")])
fit <- VAR(vardata, 6)
fevd12 <- fevd(fit,n.ahead=36)
HD12 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="ArcherDaniels", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Archer Daniels",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

allfevd <- cbind(fevd1$SSTA[,"Alico"],fevd2$SSTA[,"Conagra"],fevd3$SSTA[,"Tyson"],fevd4$SSTA[,"Monsanto"],
                 fevd5$SSTA[,"Ingredion"],fevd6$SSTA[,"Hormel"],fevd7$SSTA[,"Andersons"],fevd8$SSTA[,"Agrium"],
                 fevd9$SSTA[,"Potash"],fevd10$SSTA[,"Mosaic"],fevd11$SSTA[,"Altria"],fevd12$SSTA[,"ArcherDaniels"])

write.csv(allfevd, file="fevd_Nina.csv")

#################################################################################################
############################ Now we will examine Neutral Responses ##############################
#################################################################################################

rm(list=ls(all=TRUE))
library(vars)
library(nleqslv)
library(functional)
source('functions.R')


################################# Responses to La Nino Shocks ####################################

Data <- read.csv('Data-Monthly-Updated.csv')
Data <- Data[Data[,"SSTA"] > -0.5 & Data[,"SSTA"]<=0.5,]  ## Observations for Neutral
Data <- ts(Data, start=c(1990,1), frequency=12) # Setting data as time series by month
SSTA <- Data[,"SSTA"]
Data <- 100*diff(log(Data)) # Growth rate by taking first difference of log of variables
DATA <- cbind(SSTA,Data)
DATA <- DATA[,c("SSTA","Data.SP","Data.Energy","Data.CPI","Data.INDPRO","Data.Alico","Data.Conagra",
                "Data.Tyson","Data.Monsanto", "Data.Ingredion", "Data.Hormel","Data.Andersons",
                "Data.Agrium","Data.Potash","Data.Mosaic","Data.Altria","Data.ArcherDaniels")]

colnames(DATA) <- c("SSTA","SP","Energy","CPI","INDPRO","Alico","Conagra","Tyson","Monsanto",
                    "Ingredion","Hormel","Andersons","Agrium","Potash","Mosaic","Altria","ArcherDaniels")


#######################################################################################################

par(mfrow=c(2,2))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","SP")])
fit <- VAR(vardata, 6)
fevd1 <- fevd(fit,n.ahead=36)
HD1 <- VARhd(Estimation=fit)
res_ssta <- irf(fit,impulse=c("SSTA"), response="SSTA", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)
res_sp <- irf(fit,impulse=c("SSTA"), response="SP", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)
res_cpi <- irf(fit,impulse=c("SSTA"), response="CPI", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)
res_gdp <- irf(fit,impulse=c("SSTA"), response="INDPRO", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.95)

ts.plot(cbind(res_ssta$Lower$SSTA,res_ssta$irf$SSTA,res_ssta$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of SSTA",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

ts.plot(cbind(res_sp$Lower$SSTA,res_sp$irf$SSTA,res_sp$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of S&P Index",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

ts.plot(cbind(res_cpi$Lower$SSTA,res_cpi$irf$SSTA,res_cpi$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of CPI Inflation",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

ts.plot(cbind(res_gdp$Lower$SSTA,res_gdp$irf$SSTA,res_gdp$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Output",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

#######################################################################################################

par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Alico")])
fit <- VAR(vardata, 6)
fevd1 <- fevd(fit,n.ahead=36)
HD1 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Alico", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Alico",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Conagra")])
fit <- VAR(vardata, 6)
fevd2 <- fevd(fit,n.ahead=36)
HD2 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Conagra", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Conagra",col=c("red","black", "red"),lwd=c(2,2,2),lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Tyson")])
fit <- VAR(vardata, 6)
fevd3 <- fevd(fit,n.ahead=36)
HD3 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Tyson", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Tyson",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Monsanto")])
fit <- VAR(vardata, 6)
fevd4 <- fevd(fit,n.ahead=36)
HD4 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Monsanto", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Monsanto",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Ingredion")])
fit <- VAR(vardata, 6)
fevd5 <- fevd(fit,n.ahead=36)
HD5 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Ingredion", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Ingredion",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Hormel")])
fit <- VAR(vardata, 6)
fevd6 <- fevd(fit,n.ahead=36)
HD6 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Hormel", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Hormel",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)


par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Andersons")])
fit <- VAR(vardata, 6)
fevd7 <- fevd(fit,n.ahead=36)
HD7 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Andersons", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Andersons",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Agrium")])
fit <- VAR(vardata, 6)
fevd8 <- fevd(fit,n.ahead=36)
HD8 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Agrium", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Agrium",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Potash")])
fit <- VAR(vardata, 6)
fevd9 <- fevd(fit,n.ahead=36)
HD9 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Potash", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Potash",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Mosaic")])
fit <- VAR(vardata, 6)
fevd10 <- fevd(fit,n.ahead=36)
HD10 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Mosaic", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Mosaic",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Altria")])
fit <- VAR(vardata, 6)
fevd11 <- fevd(fit,n.ahead=36)
HD11 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Altria", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Altria",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","ArcherDaniels")])
fit <- VAR(vardata, 6)
fevd12 <- fevd(fit,n.ahead=36)
HD12 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="ArcherDaniels", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Archer Daniels",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-40,40));abline(h=0)

allfevd <- cbind(fevd1$SSTA[,"Alico"],fevd2$SSTA[,"Conagra"],fevd3$SSTA[,"Tyson"],fevd4$SSTA[,"Monsanto"],
                 fevd5$SSTA[,"Ingredion"],fevd6$SSTA[,"Hormel"],fevd7$SSTA[,"Andersons"],fevd8$SSTA[,"Agrium"],
                 fevd9$SSTA[,"Potash"],fevd10$SSTA[,"Mosaic"],fevd11$SSTA[,"Altria"],fevd12$SSTA[,"ArcherDaniels"])

write.csv(allfevd, file="fevd_Neutral.csv")

##################################################################################################################################################




