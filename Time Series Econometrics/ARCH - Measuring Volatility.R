n <-  10500 ## lots of trials
z <- rnorm(n) ## sample standard normal distribution variates
e <-  z ## store variates
y <-  z ## store again in a different place
sig2 <-  z^2 ## create volatility series
omega <-  1 ## base variance
alpha <-  0.55 ## Markov dependence on previous variance
phi <-  0.8 ## mMarkov dependence on previous period
mu <-  0.1 ## average return
omega/(1-alpha) ; sqrt(omega/(1-alpha))

set.seed("1012")
for (t in 2:n) ## Because of lag start at second date
{
  e[t] <- sqrt(sig2[t])*z[t]          ## 1. e is conditional on sig
  y[t] <-  mu + phi*(y[t-1]-mu) + e[t] ## 2. generate returns
  sig2[t+1] <-  omega + alpha * e[t]^2 ## 3. generate new sigma^2 to feed 1.
}

par(mfrow = c(2, 4))
plot(z[10001:n], type = "l", xlab = "t", 
     ylab = expression(epsilon), main = "1. Simple noise")
plot(sqrt(sig2[10000:n]), type = "l", 
     xlab = "t", ylab = expression(sigma[t]), 
     main = "2. Conditional sigma")
plot(e[10001:n], type = "l", xlab = "t", 
     ylab = "a", main = "3. ARCH")
plot(y[10001:n], type = "l", xlab = "t", 
     ylab = "y", main = "4. AR+ARCH")
acf(e[10001:n], main = "5. ARCH")
acf(abs(e[10001:n]), main = "6. Absolute ARCH value")
acf(y[10001:n], main = "7. AR+ARCH")
acf(abs(y[10001:n]), main = "8. Absolute AR+ARCH value")

require(rugarch)
require(qrmdata)
require(xts)

data("EUR_USD")
data("GBP_USD")