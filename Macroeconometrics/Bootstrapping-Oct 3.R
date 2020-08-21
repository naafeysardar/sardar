# Block Bootstrap 

block.obs <- function(block, k) {
  return( (block*k - (k-1)): (block*k) )      # k=block length, block=block number
}

block.obs(3, 4)   # block no. 3 and block length is 4. 
block.obs(10, 4)

# Step 1: Choose 50 Blocks

block.numbers <- sample(1:50, replace = TRUE)
block.numbers

# step 2: Get the Observation Numbers for all Blocks

obs.numbers <- lapply(block.numbers, block.obs, k=4)
obs.numbers

# Step 3: Convert list to vector

obs <- unlist(obs.numbers)
obs

# Write above 3 steps as one function

obs.numbers <- function(nblocks, len){
  block.numbers <- sample(1:nblocks, replace = TRUE)
  unlist(lapply(block.numbers, block.obs, k=len))
}

obs.numbers(50, 4)

# Upload Data
 
data.raw <- read.csv("u-inf.csv", header = TRUE)

dataset <- ts(data.raw, start = c(1948,1), frequency = 12)

u <- dataset[,"u"]

inf <- dataset[,"inf"]

# Estimate model using our sample

fit <- tsreg(inf, u)
summary(fit)

coefficients(summary(fit))

# Generate new sample of inflation

res <- residuals(fit)

length(inf) # Total observations are 847

# Generate observation number for our sampling

obs <- obs.numbers(211, 4)
res.first <- res[obs]

res.sim <- c(res.first, res[845:847])
res.sim

# Convert res.sim into time series

res.sim <- ts(res.sim, start = start(u), frequency=12)
plot(res.sim)

# Calculate the simulated inflation data

inf.sim <- 0.24 + 0.007*u + res.sim

# Reestimate the regression 

fit.sim <- tsreg(inf.sim, u)
fit.sim

library(tstools)

set.seed(200)
sim.betas <- replicate(100, {
  obs <- obs.numbers(211, 4)
  res.sim <- ts(c(res[obs], res[845:847]), start=start(u), 
                frequency=12)
  inf.sim <- 0.24 + 0.007*u + res.sim
  fit.sim <- tsreg(inf.sim, u)
  coefficients(fit.sim)[2]
})

sd(sim.betas)

# Block Length of 8

set.seed(200)
sim.betas <- replicate(100, {
  obs <- obs.numbers(105, 8)
  res.sim <- ts(c(res[obs], res[841:847]), start=start(u), 
                frequency=12)
  inf.sim <- 0.24 + 0.007*u + res.sim
  fit.sim <- tsreg(inf.sim, u)
  coefficients(fit.sim)[2]
})

sd(sim.betas)


# Block + Pairs (Deal with serial correlation & heteroskedasticity)

set.seed(200)
sim.betas <- replicate(100, {
  obs <- obs.numbers(211, 4)
  inf.sim <- ts(inf[obs])     # Resampling Data
  u.sim <- ts(u[obs])         # Resampling Data
  fit.sim <- tsreg(inf.sim, u.sim)
  coefficients(fit.sim)[2]
})

sd(sim.betas)


sample(c(-1, 1), size = length(u), replace = TRUE)






