rm(list=ls(all=TRUE))
library(vars)
library(nleqslv)
library(functional)
source('functions.R')

Data <- read.csv('Data-Monthly-Updated.csv')
Data <- ts(Data, start=c(1990,1), frequency=12) # Setting data as time series by month
SSTA <- Data[,"SSTA"]
Data <- 100*diff(log(Data)) # Growth rate by taking first difference of log of variables
DATA <- cbind(SSTA,Data)
DATA <- DATA[,c("SSTA","Data.SP","Data.Energy","Data.CPI","Data.INDPRO","Data.Alico","Data.Conagra",
                "Data.Tyson","Data.Monsanto", "Data.Ingredion", "Data.Hormel","Data.Andersons",
                "Data.Agrium","Data.Potash","Data.Mosaic","Data.Altria","Data.ArcherDaniels")]

colnames(DATA) <- c("SSTA","SP","Energy","CPI","INDPRO","Alico","Conagra","Tyson","Monsanto",
                    "Ingredion","Hormel","Andersons","Agrium","Potash","Mosaic","Altria","ArcherDaniels")

############################# Plot El Nino and La Nina thresholds ################################

plot(SSTA, type="l", xlab="Year", ylab="3-Month Niño Region 3.4 Average",lwd=2,yaxt="n")
axis(side=2, at=seq(-2.5,2.5,0.5))

abline(h=0, lwd=2); abline(h=-0.5,col="blue", lty=2); abline(h=-1,col="blue", lty=2)
abline(h=-1.5,col="blue", lty=2); abline(h=-2,col="blue", lty=2)

abline(h=0.5,col="red", lty=2); abline(h=1,col="red", lty=2)
abline(h=1.5,col="red", lty=2); abline(h=2,col="red", lty=2)

text(1990,2.2, "Very Strong El Niño", col = "red", adj = c(-.1, -.1))
text(2000,1.65, "Strong El Niño", col = "red", adj = c(-.1, -.1))
text(2003,1.2, "Moderate El Niño", col = "red", adj = c(-.1, -.1))
text(2009.8,0.65, "Weak El Niño", col = "red", adj = c(-.1, -.1))
text(1998.7,0.2, "Neutral", col = "green", adj = c(-.1, -.1))
text(1998.7,-0.25, "Neutral", col = "green", adj = c(-.1, -.1))
text(1990,-1.85, "Strong La Niña", col = "blue", adj = c(-.1, -.1))
text(1990,-1.3, "Moderate La Niña", col = "blue", adj = c(-.1, -.1))
text(1990,-0.8, "Weak La Niña", col = "blue", adj = c(-.1, -.1))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","SP")])

############# Impulse Response Functions
 
fit <- VAR(vardata, 6)

res_ssta <- irf(fit,impulse=c("SSTA"), response=c("SSTA"), n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
res_sp <- irf(fit,impulse=c("SSTA"), response="SP", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
res_cpi <- irf(fit,impulse=c("SSTA"), response="CPI", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
res_indpro <- irf(fit,impulse=c("SSTA"), response="INDPRO", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
res_energy <- irf(fit,impulse=c("SSTA"), response="Energy", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)

par(mfrow=c(2,2))

ts.plot(cbind(res_ssta$Lower$SSTA[,"SSTA"],res_ssta$irf$SSTA[,"SSTA"],res_ssta$Upper$SSTA[,"SSTA"]),xlab="Months", ylab="percent",
        main="Response of SSTA",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(0,4));abline(h=0)

ts.plot(cbind(res_sp$Lower$SSTA,res_sp$irf$SSTA,res_sp$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of S&P Index",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-4,4));abline(h=0)

ts.plot(cbind(res_cpi$Lower$SSTA,res_cpi$irf$SSTA,res_cpi$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of CPI Inflation",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-0.2,0.2));abline(h=0)

ts.plot(cbind(res_indpro$Lower$SSTA,res_indpro$irf$SSTA,res_indpro$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Output",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-1,1));abline(h=0)

ts.plot(cbind(res_energy$Lower$SSTA,res_energy$irf$SSTA,res_energy$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Oil Prices",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-8,4));abline(h=0)

#######################################################################################################

par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Alico")])
fit <- VAR(vardata, 6)
fevd1 <- fevd(fit,n.ahead=36)
HD1 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Alico", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Alico",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Conagra")])
fit <- VAR(vardata, 6)
fevd2 <- fevd(fit,n.ahead=36)
HD2 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Conagra", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Conagra",col=c("red","black", "red"),lwd=c(2,2,2),lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Tyson")])
fit <- VAR(vardata, 6)
fevd3 <- fevd(fit,n.ahead=36)
HD3 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Tyson", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Tyson",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Monsanto")])
fit <- VAR(vardata, 6)
fevd4 <- fevd(fit,n.ahead=36)
HD4 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Monsanto", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Monsanto",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Ingredion")])
fit <- VAR(vardata, 6)
fevd5 <- fevd(fit,n.ahead=36)
HD5 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Ingredion", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Ingredion",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Hormel")])
fit <- VAR(vardata, 6)
fevd6 <- fevd(fit,n.ahead=36)
HD6 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Hormel", n.ahead=18, cumulative = TRUE, runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Hormel",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)


par(mfrow=c(2,3))

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Andersons")])
fit <- VAR(vardata, 6)
fevd7 <- fevd(fit,n.ahead=36)
HD7 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Andersons", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Andersons",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Agrium")])
fit <- VAR(vardata, 6)
fevd8 <- fevd(fit,n.ahead=36)
HD8 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Agrium", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Agrium",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Potash")])
fit <- VAR(vardata, 6)
fevd9 <- fevd(fit,n.ahead=36)
HD9 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Potash", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Potash",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Mosaic")])
fit <- VAR(vardata, 6)
fevd10 <- fevd(fit,n.ahead=36)
HD10 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Mosaic", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Mosaic",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","Altria")])
fit <- VAR(vardata, 6)
fevd11 <- fevd(fit,n.ahead=36)
HD11 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="Altria", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Altria",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

vardata <- na.omit(DATA[,c("SSTA","INDPRO","CPI","ArcherDaniels")])
fit <- VAR(vardata, 6)
fevd12 <- fevd(fit,n.ahead=36)
HD12 <- VARhd(Estimation=fit)
res_alico <- irf(fit,impulse=c("SSTA"), response="ArcherDaniels", n.ahead=18, cumulative = TRUE,runs=1000, ci=0.90)
ts.plot(cbind(res_alico$Lower$SSTA,res_alico$irf$SSTA,res_alico$Upper$SSTA),xlab="Months", ylab="percent",
        main="Response of Archer Daniels",col=c("red","black", "red"),lwd=c(2,2,2), lty=c(2,1,2),ylim=c(-10,10));abline(h=0)

allfevd <- cbind(fevd1$SSTA[,"Alico"],fevd2$SSTA[,"Conagra"],fevd3$SSTA[,"Tyson"],fevd4$SSTA[,"Monsanto"],
                 fevd5$SSTA[,"Ingredion"],fevd6$SSTA[,"Hormel"],fevd7$SSTA[,"Andersons"],fevd8$SSTA[,"Agrium"],
                 fevd9$SSTA[,"Potash"],fevd10$SSTA[,"Mosaic"],fevd11$SSTA[,"Altria"],fevd12$SSTA[,"ArcherDaniels"])

write.csv(allfevd, file="fevd.csv")

##### Make historical decompositions time-series

HD1 <- ts(na.omit(HD1[,,1][,4]), start=c(1990,8), frequency=12); HD1 <- na.omit(cbind(HD1, DATA[,"Alico"])) 
HD2 <- ts(na.omit(HD2[,,1][,4]), start=c(1990,8), frequency=12); HD2 <- na.omit(cbind(HD2, DATA[,"Conagra"])) 
HD3 <- ts(na.omit(HD3[,,1][,4]), start=c(1990,8), frequency=12); HD3 <- na.omit(cbind(HD3, DATA[,"Tyson"])) 
HD4 <- ts(na.omit(HD4[,,1][,4]), start=c(2001,5), frequency=12); HD4 <- na.omit(cbind(HD4, DATA[,"Monsanto"]))  
HD5 <- ts(na.omit(HD5[,,1][,4]), start=c(1998,7), frequency=12); HD5 <- na.omit(cbind(HD5, DATA[,"Ingredion"])) 
HD6 <- ts(na.omit(HD6[,,1][,4]), start=c(1990,8), frequency=12); HD6 <- na.omit(cbind(HD6, DATA[,"Hormel"]))  
HD7 <- ts(na.omit(HD7[,,1][,4]), start=c(1996,9), frequency=12); HD7 <- na.omit(cbind(HD7, DATA[,"Andersons"])) 
HD8 <- ts(na.omit(HD8[,,1][,4]), start=c(1995,12), frequency=12); HD8 <- na.omit(cbind(HD8, DATA[,"Agrium"]))  
HD9 <- ts(na.omit(HD9[,,1][,4]), start=c(1990,8), frequency=12); HD9 <- na.omit(cbind(HD9, DATA[,"Potash"]))  
HD10 <- ts(na.omit(HD10[,,1][,4]), start=c(1990,8), frequency=12); HD10 <- na.omit(cbind(HD10, DATA[,"Mosaic"]))  
HD11 <- ts(na.omit(HD11[,,1][,4]), start=c(1990,8), frequency=12); HD11 <- na.omit(cbind(HD11, DATA[,"Altria"])) 
HD12 <- ts(na.omit(HD12[,,1][,4]), start=c(1990,8), frequency=12); HD12 <- na.omit(cbind(HD12, DATA[,"ArcherDaniels"])) 


#######   Plot Quarterly Averages
par(mfrow=c(3,2))
ts.plot(HD1,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Alico",ylab="", xlab="")
ts.plot(HD2,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Conagra",ylab="", xlab="")
ts.plot(HD3,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Tyson",ylab="", xlab="")
ts.plot(HD4,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Monsanto",ylab="", xlab="")
ts.plot(HD5,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Ingredion",ylab="", xlab="")
ts.plot(HD6,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Hormel",ylab="", xlab="")

par(mfrow=c(3,2))

ts.plot(HD7,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Andersons",ylab="", xlab="")
ts.plot(HD8,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Agrium",ylab="", xlab="")
ts.plot(HD9,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Potash",ylab="", xlab="")
ts.plot(HD10,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Mosaic",ylab="", xlab="")
ts.plot(HD11,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Altria",ylab="", xlab="")
ts.plot(HD12,gpars=list(lty=c(1,2), lwd=c(2,1)),ylim=c(-30,30),
        main="Contribution to Archer Daniels",ylab="", xlab="")


