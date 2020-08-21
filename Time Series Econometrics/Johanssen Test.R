library(urca)
data("Raotbl3")
vec <- ca.jo(Raotbl3[, c("li","lc")])

## ca.jo stands for conintegration analysis johannsen test ##

summary(vec)
dim(Raotbl3)

vec <- ca.jo(Raotbl3[, c("li","lc")], ecdet = "const")
summary(vec)

vec <- ca.jo(Raotbl3[, c("li","lc")], ecdet = "trend")
summary(vec)

## We estimated cointegration ##

lc <- ts(Raotbl3$lc, start=c(1966,4), frequency=4)
li <- ts(Raotbl3$li, start=c(1966,4), frequency=4)

## Take differenced VAR ##

data.diff <- ts.intersect(diff(lc),diff(li))
z <- li - 0.98*lc
plot(z)

## Z looks stationery ##

tsp(z)

library(tstools)
z1 <- lags(z,1)     ## Create lags for z1 ##
z1

## Put z term in VAR to estimate VECM ##

est.data <- ts.intersect(data.diff, z1)
colnames(est.data)

library(vars)

VAR(est.data[, c("data.diff.diff(lc)","data.diff.diff(li)")], exogen = est.data[,"z1"])

## Create 2x2 Matrix ##

a <- matrix(NA, ncol=2, nrow = 2)
a

a[2,1] <- 0   ##Imposing restrcition##

b <- matrix(NA, ncol=2, nrow=2)

svec.fit <- SVEC(vec, SR=a, LR=b)
summary(svec.fit)
irf(svec.fit)
plot(irf(svec.fit))
