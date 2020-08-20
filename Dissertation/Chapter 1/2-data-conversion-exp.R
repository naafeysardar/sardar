##### Upload Packages #####

rm(list = ls())

project.raw <- read.csv("full-dataset.csv", header=TRUE, fileEncoding="UTF-8-BOM")
data <- ts(project.raw, start=c(1972,4), frequency=4)

gdp <- data[,"rgdp"]
auerbach <- data[,"auerbach"]
hamilton <- data[,"hamilton"]

price <- data[,"price"]/100

miles <- data[,"miles"]
rpg <- data[,"rpg"]

n.pce <- data[,"n.pce"]
n.durable <- data[,"n.durable"]
n.nondurable <- data[,"n.nondurable"]
n.services <- data[,"n.services"]

n.motor <- data[,"n.motor"]

n.furnishing <- data[,"n.furnishing"]
n.recreational <- data[,"n.recreational"]
n.otherdurable <- data[,"n.otherdurable"]
n.food <- data[,"n.food"]
n.gasoline <- data[,"n.gasoline"]
n.clothing <- data[,"n.clothing"]
n.othernondurable <- data[,"n.othernondurable"]

n.housing.util <- data[,"n.housing.utilities"]
n.housing <- data[,"n.housing"]
n.utilities <- data[,"n.utilities"]
n.health <- data[,"n.healthcare"]
n.outpatient <- data[,"n.outpatient"]
n.hospital <- data[,"n.hospital"]
n.transport <- data[,"n.transport"]
n.public <- data[,"n.public.transport"]
n.motor.services <- data[,"n.motor.services"]
n.recreation <- data[,"n.recreation"]
n.food.services <- data[,"n.food.services"]
n.financial <- data[,"n.financial"]
n.other.services <- data[,"n.other.services"]

r.pce <- n.pce/price 
pce <- 100*(r.pce/stats::lag(r.pce,-1) - 1)

r.durable <- n.durable/price
durable <- 100*(r.durable/stats::lag(r.durable,-1) - 1)

r.nondurable <- n.nondurable/price
nondurable <- 100*(r.nondurable/stats::lag(r.nondurable,-1) - 1)

r.services <- n.services/price
services <- 100*(r.services/stats::lag(r.services,-1) - 1)

r.gas <- n.gasoline/price
gas <- 100*(r.gas/stats::lag(r.gas,-1) - 1)

inter <- auerbach*gas

r.motor <- n.motor/price
motor <- 100*(r.motor/stats::lag(r.motor,-1) - 1)

r.recreational <- n.recreational/price
recreational <- 100*(r.recreational/stats::lag(r.recreational,-1) - 1)

r.furnishing <- n.furnishing/price
furnishing <- 100*(r.furnishing/stats::lag(r.furnishing,-1) - 1)

r.otherdurable <- n.otherdurable/price
otherdurable <- 100*(r.otherdurable/stats::lag(r.otherdurable,-1) - 1)

r.food <- n.food/price
food <- 100*(r.food/stats::lag(r.food,-1) - 1)

r.clothing <- n.clothing/price
clothing <- 100*(r.clothing/stats::lag(r.clothing,-1) - 1)

r.othernondurable <- n.othernondurable/price
othernondurable <- 100*(r.othernondurable/stats::lag(r.othernondurable,-1) - 1)

r.housing.util <- n.housing.util/price
housing.util <- 100*(r.housing.util/stats::lag(r.housing.util,-1) - 1)

r.housing <- n.housing/price
housing <- 100*(r.housing/stats::lag(r.housing,-1) - 1)

r.utilities <- n.utilities/price
utilities <- 100*(r.utilities/stats::lag(r.utilities,-1) - 1)

r.health <- n.health/price
health <- 100*(r.health/stats::lag(r.health,-1) - 1)

r.outpatient <- n.outpatient/price
outpatient <- 100*(r.outpatient/stats::lag(r.outpatient,-1) - 1)

r.hospital <- n.hospital/price
hospital <- 100*(r.hospital/stats::lag(r.hospital,-1) - 1)

r.transport <- n.transport/price
transport <- 100*(r.transport/stats::lag(r.transport,-1) - 1)

r.public <- n.public/price
public <- 100*(r.public/stats::lag(r.public,-1) - 1)

r.motor.services <- n.motor.services/price
motor.services <- 100*(r.motor.services/stats::lag(r.motor.services, -1) - 1)

r.recreation <- n.recreation/price
recreation <- 100*(r.recreation/stats::lag(r.recreation,-1) - 1)

r.food.services <- n.food.services/price
food.services <- 100*(r.food.services/stats::lag(r.food.services,-1) - 1)

r.financial <- n.financial/price
financial <- 100*(r.financial/stats::lag(r.financial,-1) - 1)

r.other.services <- n.other.services/price
other.services <- 100*(r.other.services/stats::lag(r.other.services,-1) - 1)

n.savings <- data[,"n.savings"]

r.savings <- n.savings/price
savings <- 100*(r.savings/stats::lag(r.savings,-1) - 1)

length(n.pce)
length(n.gasoline)
t.pce <- n.pce + n.gasoline

s.pce <- n.pce/t.pce
s.gasoline <- n.gasoline/t.pce

par(mar = c(5, 5, 3, 5))
plot(s.pce, type ="l", ylab = "Non-Gasoline Share of Consumption", xlab = "Year", col = "blue", ylim=c(0.9, 1))
par(new = TRUE)
plot(s.gasoline, type = "l", xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "red", ylim=c(0.01, 0.1))
axis(side = 4)
mtext("Gasoline Share of Consumption", side = 4, line = 3)
legend("topleft",
       c("Share of Non-Gasoline Consumption","Share of Gasoline Consumption"),
       fill=c("blue","red"))
