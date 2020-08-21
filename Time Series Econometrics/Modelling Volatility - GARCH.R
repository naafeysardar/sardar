library(quantmod)
library(ggplot2)
library(broom)
getSymbols(c('SPY'))

SPY              <- SPY$SPY.Adjusted
SPYRet           <- log(SPY) - log(lag(SPY))
SPYRet_xts       <- SPYRet
colnames(SPYRet) <- c('SPY')
SPYRet           <- tidy(SPYRet)

ggplot(SPYRet, aes(x = index, y = value, color = series)) + 
  geom_line() + 
  theme_bw() +
  labs(title = "SPY Returns Returns from 2007 to 2017", x = "")

library(tibble)
SPYRet <- add_column(SPYRet, SquaredReturns = SPYRet$value^2, AbsoluteReturns = abs(SPYRet$value))

ggplot(SPYRet, aes(x = index, y = AbsoluteReturns, color = series)) + 
  geom_line() + 
  theme_bw() +
  labs(title = "SPY Absolute Value of Returns from 2007 to 2017", x = "")

ggplot(SPYRet, aes(x = index, y = SquaredReturns, color = series)) + 
  geom_line() + 
  theme_bw() +
  labs(title = "SPY Squared Returns from 2007 to 2017", x = "")

library(rugarch)

garch11        <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                             mean.model = list(armaOrder = c(1, 1), include.mean = TRUE), 
                             distribution.model = "norm")

garchfit       <- ugarchfit(spec = garch11, data = SPYRet_xts["2007-02-01/"], solver = "hybrid")
garchfit

spec           <- getspec(garchfit)
setfixed(spec) <- as.list(coef(garchfit))
garchforecast1 <- ugarchforecast(spec, n.ahead = 1, n.roll = 2499, data = SPYRet_xts["2007-02-01/"], out.sample = 2500)

plot(garchforecast1, which = 4)