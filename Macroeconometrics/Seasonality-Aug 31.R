library(seasonal)

unemp

plot(unemp)

adjunemp <- seas(unemp)
adjunemp$data

window(unemp, start=c(1990, 1), end=c(1990,12))
  
window(adjunemp$data, start=c(1990, 1), end=c(1990,12))

plot(adjunemp)
