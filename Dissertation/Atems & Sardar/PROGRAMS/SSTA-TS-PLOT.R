rm(list=ls(all=TRUE))

Data <- read.csv('OTHER-COMPANIES-WEEKLY.csv')
Data <- ts(Data, start=c(1990,1), frequency=52) # Setting data as time series by month
SSTA <- Data[,"SSTA"]

############################# Plot El Nino and La Nina thresholds ################################

plot(SSTA, type="l", xlab="Year", ylab="3-Month Niño Region 3.4 Average",lwd=2,yaxt="n")
axis(side=2, at=seq(-3,3,0.5))
abline(h=0, lwd=2); abline(h=-0.5,col="blue", lty=2); abline(h=0.5,col="red", lty=2)

text(2003,2, "El Niño", col = "red", adj = c(-.1, -.1))
text(1998.7,0.2, "Neutral", col = "green", adj = c(-.1, -.1))
text(1998.7,-0.25, "Neutral", col = "green", adj = c(-.1, -.1))
text(2002,-1.5, "La Niña", col = "blue", adj = c(-.1, -.1))
