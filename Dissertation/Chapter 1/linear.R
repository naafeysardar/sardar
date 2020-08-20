library(tstools)
library(vars)

project.raw <- read.csv("full-dataset.csv", header=TRUE)
tsdata <- lapply(names(project.raw), function(n) {
  ts(project.raw[, n], start=c(1972,4), frequency=4)
})
names(tsdata) <- names(project.raw)

newtsdata <- within(tsdata, {
  gas <- 100*pctChange(n.gasoline/price)
  inter <- auerbach*gas
  pce <- 100*pctChange(n.pce/price)
  durable <- 100*pctChange(n.durable/price)
  nondurable <- 100*pctChange(n.nondurable/price)
  services <- 100*pctChange(n.services/price)
  motor <- 100*pctChange(n.motor/price)
  recreational <- 100*pctChange(n.recreational/price)
  furnishing <- 100*pctChange(n.furnishing/price)
  otherdurable <- 100*pctChange(n.otherdurable/price)
  food <- 100*pctChange(n.food/price)
  clothing <- 100*pctChange(n.clothing/price)
  othernondurable <- 100*pctChange(n.othernondurable/price)
  housing <- 100*pctChange(n.housing/price)
  utilities <- 100*pctChange(n.utilities/price)
  outpatient <- 100*pctChange(n.outpatient/price)
  hospital <- 100*pctChange(n.hospital/price)
  transport <- 100*pctChange(n.transport/price)
  public <- 100*pctChange(n.public.transport/price)
  motor.services <- 100*pctChange(n.motor.services/price)
  recreation <- 100*pctChange(n.recreation/price)
  food.services <- 100*pctChange(n.food.services/price)
  financial <- 100*pctChange(n.financial/price)
  other.services <- 100*pctChange(n.other.services/price)
})

par(mfrow=c(1,4))

irf.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  dataset <- ts.intersect(newtsdata$gas, macro)
  rf <- VAR(dataset, p=2)
  res <- irf(rf, impulse="newtsdata.gas", response="macro", n.ahead=5, cumulative = TRUE, runs=100, ci=0.95)
  n=length(res$irf$newtsdata.gas)
  for(i in 1:n){res$irf$newtsdata.gas[i]=res$irf$newtsdata.gas[i]*1.49}
  for(i in 1:n){res$Lower$newtsdata.gas[i]=res$Lower$newtsdata.gas[i]*1.49}
  for(i in 1:n){res$Upper$newtsdata.gas[i]=res$Upper$newtsdata.gas[i]*1.49}
  ts.plot(cbind(res$Lower$newtsdata.gas, res$irf$newtsdata.gas, res$Upper$newtsdata.gas), xlab="Quarters", ylab="% Change",
          main=varname, col=c("black","black", "black"), lwd=c(1,1,1), lty=c(2,1,2))
  abline(h=0)
  return(list(irf=res$irf$newtsdata.gas, corr=cor(residuals(rf))[1,2], varname=varname))
}

varinfo <- list(
  list(newtsdata$pce, "PCE"),
  list(newtsdata$durable, "Durables"),
  list(newtsdata$nondurable, "Nondurables"),
  list(newtsdata$services, "Services"))
irfs <- lapply(varinfo, irf.linear)
names(irfs) <- sapply(irfs, function(z) { z$varname })

local({
  services <- 0.64*irfs$Services$irf
  nondurables <- 0.23*irfs$Nondurables$irf
  durables <- 0.13*irfs$Durables$irf
  agg <- services + nondurables + durables
  print(agg)
  
  lapply(list(
    list(agg, "Implied PCE"),
    list(services, "Services"),
    list(nondurables, "Nondurables"),
    list(durables, "Durables")), function(z) {
      ts.plot(z[[1]], xlab="Quarters", ylab="% Change", main=z[[2]], lwd=1.3, ylim=c(-1,1))
      abline(h=0)
    })
  mtext("Contribution to Response of PCE", outer=TRUE, cex=1.5)
})

par(mfrow=c(3,3))

varinfo <- list(
  list(newtsdata$motor, "Motor Vehicles"),
  list(newtsdata$recreation, "Recreational"),
  list(newtsdata$furnishing, "Furnishing"),
  list(newtsdata$otherdurable, "Other Durables"),
  list(newtsdata$food, "Food"),
  list(newtsdata$clothing, "Clothing"),
  list(newtsdata$othernondurable, "Other Nondurables"),
  list(newtsdata$housing, "Housing"),
  list(newtsdata$utilities, "Utilities"),
  list(newtsdata$public, "Public Transport"),
  list(newtsdata$motor.services, "Motor Services"),
  list(newtsdata$recreation, "Recreation Services"),
  list(newtsdata$food.services, "Food Services"),
  list(newtsdata$outpatient, "Outpatient Services"),
  list(newtsdata$hospital, "Hospital & Nursing Homes"),
  list(newtsdata$financial, "Financial"),
  list(newtsdata$other.services, "Other Services"))
irfs <- lapply(varinfo, irf.linear)
names(irfs) <- sapply(irfs, function(z) { z$varname })
