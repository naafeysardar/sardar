# Running an AR(3) Model #

ARMAtoMA(ar=c(1.2, -0.2, -0.2), lag.max = 10)

# Computes IRFs 10 periods in to the future. Half-life #
#               is 5 time periods.                     #


# Running an AR(1) Mu_1odel #

ARMAtoMA(ar=0.9, lag.max = 30)

ARMAtoMA(ar=0.7, lag.max = 30)

# AR(3) #

ARMAtoMA(ar=c(1.2, -0.2, -0.6), lag.max = 10) 

# Half life is 4 time periods. #


