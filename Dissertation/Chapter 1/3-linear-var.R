library(tstools)
library(vars)
library(tsDyn)

par(mfrow=c(2,2))

irf.linear <- function(obj) {
	macro <- obj[[1]]
	varname <- obj[[2]]
	dataset <- ts.intersect(gas, macro)
	rf <- VAR(dataset, p=2)
	res <- irf(rf, impulse="gas", response="macro", n.ahead=5, cumulative = TRUE, runs=100, ci=0.95)
	n=length(res$irf$gas)
	for(i in 1:n){res$irf$gas[i]=res$irf$gas[i]*1.63}
	for(i in 1:n){res$Lower$gas[i]=res$Lower$gas[i]*1.63}
	for(i in 1:n){res$Upper$gas[i]=res$Upper$gas[i]*1.63}
	ts.plot(cbind(res$Lower$gas, res$irf$gas, res$Upper$gas), xlab="Quarters", ylab="% Change",
		main=varname, col=c("black","black", "black"), lwd=c(1,1,1), lty=c(2,1,2), ylim=c(-5,5))
	abline(h=0)
	return(list(irf=res$irf$gas, rf, corr=cor(residuals(rf))[1,2], varname=varname))
}

varinfo <- list(
	list(pce, "PCE"),
	list(durable, "Durables"),
  list(nondurable, "Nondurables"),
	list(services, "Services"))
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

###################################################################################################################################################################

par(mfrow=c(2,3))

irf.linear <- function(obj) {
  macro <- obj[[1]]
  varname <- obj[[2]]
  dataset <- ts.intersect(gas, macro)
  rf <- VAR(dataset, p=2)
  res <- irf(rf, impulse="gas", response="macro", n.ahead=5, cumulative = TRUE, runs=100, ci=0.95)
  n=length(res$irf$gas)
  for(i in 1:n){res$irf$gas[i]=res$irf$gas[i]*1.63}
  for(i in 1:n){res$Lower$gas[i]=res$Lower$gas[i]*1.63}
  for(i in 1:n){res$Upper$gas[i]=res$Upper$gas[i]*1.63}
  ts.plot(cbind(res$Lower$gas, res$irf$gas, res$Upper$gas), xlab="Quarters", ylab="% Change",
          main=varname, col=c("black","black", "black"), lwd=c(1,1,1), lty=c(2,1,2))
  abline(h=0)
  return(list(irf=res$irf$gas, corr=cor(residuals(rf))[1,2], varname=varname))
}

varinfo <- list(
  list(motor, "Motor Vehicles"),
  list(furnishing, "Furnishing"),
  list(otherdurable, "Other Durables"),
  list(food, "Food"),
  list(clothing, "Clothing"),
  list(othernondurable, "Other Nondurables"),
  list(utilities, "Housing & Utilities"),
  list(public, "Public Transportation"),
  list(motor.services, "Motor Services"),
  list(food.services, "Food Services"),
  list(other.services, "Other Services"))
irfs <- lapply(varinfo, irf.linear)
names(irfs) <- sapply(irfs, function(z) { z$varname })
