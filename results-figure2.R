## Source nhppCDF function (CDF for non-stationary process)
(source("nhppCDF.R"))

## Source constant environment
(source("Environments/constant.R"))
## Source increasing1 environment
(source("Environments/increasing1.R"))
## Source increasing2 environment
(source("Environments/increasing2.R"))
## Source fluctuating1 environment
(source("Environments/fluctuating1.R"))
## Source fluctuating2 environment
(source("Environments/fluctuating2.R"))

# Exponential growth
## Source environment-dependent function for exponential growth
(source("Demographic functions/exponentialfun.R"))
## Source intensity function
(source("Intensity functions/exponentialINT.R"))
# Logistic growth
## Source environment-dependent function for logistic growth
(source("Demographic functions/logisticfun.R"))
## Source intensity function
(source("Intensity functions/logisticINT.R"))

# Global values
init <- 100 # initial population size
times <- seq(0, 10, .1) # time points

# Save plot
## pdf(file="Figures/fig2.pdf", width = 7, height = 7)
# Set layout
par(mfrow = c(2, 3), mar = c(5, 6, 2, 1), lend = 2, cex.axis = 1.3, cex.lab = 1.5)

##### EXPONENTIAL GROWTH
########## Increasing
# Parameters, exponential growth (increasing1)
param1 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = increasing1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Parameters, exponential growth (increasing2)
param2 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = increasing2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(exponentialINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhppCDF(times[i], init, param1, exponentialINT)
    ssaplus2.cdf[i] <- nhppCDF(times[i], init, param2, exponentialINT)
}
# Plots
plot(NA, xlab = "", ylab = "",
     xlim = c(0, 5), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
title(ylab = "Probability event has occurred", line = 4)
axis(1, at = seq(0, 5, 1), tcl = -.5)
axis(2, at = seq(0, 1, .2), tcl = -.5, las = 1)
points(ssa.cdf ~ times, type = "l", col = "black", lwd = 2)
points(ssaplus1.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dashed")
points(ssaplus2.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dotted")
mtext(side = 3, c("(a)"), adj = 0)
legend("bottomright",  c("SSAn", "SSA+ (inc. 1)", "SSA+ (inc. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)

########## Fluctuating
# Parameters, exponential growth (fluctuating1)
param1 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = fluctuating1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Parameters, exponential growth (fluctuating2)
param2 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = fluctuating2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(exponentialINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhppCDF(times[i], init, param1, exponentialINT)
    ssaplus2.cdf[i] <- nhppCDF(times[i], init, param2, exponentialINT)
}
## Plots
plot(NA, xlab = "", ylab = "",
     xlim = c(0, 5), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
axis(1, at = seq(0, 5, 1), tcl = -.5)
axis(2, at = seq(0, 1, .2), tcl = -.5, las = 1)
points(ssa.cdf ~ times, type = "l", col = "black", lwd = 2)
points(ssaplus1.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dashed")
points(ssaplus2.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dotted")
mtext(side = 3, c("(b)"), adj = 0)
legend("bottomright",  c("SSAn", "SSA+ (fluc. 1)", "SSA+ (fluc. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)

########## Random
## Create random1 environmental values
set.seed(20170731)
randenv1 <- rnorm(length(times), mean = 1.5, sd = .25)
## Create function that interpolates these values
randenvfun1 <- splinefun(times, randenv1)
## Create random2 environmental values
set.seed(20170731)
randenv2 <- rnorm(length(times), mean = 1.5, sd = .5)
## Create function that interpolates these values
randenvfun2 <- splinefun(times, randenv2)
# Parameters, exponential growth (random1)
param1 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = randenvfun1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Parameters, exponential growth (random2)
param2 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = randenvfun2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(exponentialINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhppCDF(times[i], init, param1, exponentialINT)
    ssaplus2.cdf[i] <- nhppCDF(times[i], init, param2, exponentialINT)
}
# Plots
plot(NA, xlab = "", ylab = "",
     xlim = c(0, 5), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
axis(1, at = seq(0, 5, 1), tcl = -.5)
axis(2, at = seq(0, 1, .2), tcl = -.5, las = 1)
points(ssa.cdf ~ times, type = "l", col = "black", lwd = 2)
points(ssaplus1.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dashed")
points(ssaplus2.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dotted")
mtext(side = 3, c("(c)"), adj = 0)
legend("bottomright",  c("SSAn", "SSA+ (rand. 1)", "SSA+ (rand. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)

##### LOGISTIC GROWTH
########## Increasing
# Parameters, logistic growth (increasing1)
param1 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = increasing1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Parameters, logistic growth (increasing2)
param2 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = increasing2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(logisticINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhppCDF(times[i], init, param1, logisticINT)
    ssaplus2.cdf[i] <- nhppCDF(times[i], init, param2, logisticINT)
}
# Plots
plot(NA, xlab = "Time", ylab = "",
     xlim = c(0, 5), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
title(ylab = "Probability event has occurred", line = 4)
axis(1, at = seq(0, 5, 1), tcl = -.5)
axis(2, at = seq(0, 1, .2), tcl = -.5, las = 1)
points(ssa.cdf ~ times, type = "l", col = "black", lwd = 2)
points(ssaplus1.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dashed")
points(ssaplus2.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dotted")
mtext(side = 3, c("(d)"), adj = 0)
legend("bottomright",  c("SSAn", "SSA+ (inc. 1)", "SSA+ (inc. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)

########## Fluctuating
# Parameters, logistic growth (fluctuating1)
param1 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = fluctuating1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Parameters, logistic growth (fluctuating2)
param2 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = fluctuating2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(logisticINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhppCDF(times[i], init, param1, logisticINT)
    ssaplus2.cdf[i] <- nhppCDF(times[i], init, param2, logisticINT)
}
# Plot
plot(NA, xlab = "Time", ylab = "",
     xlim = c(0, 5), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
axis(1, at = seq(0, 5, 1), tcl = -.5)
axis(2, at = seq(0, 1, .2), tcl = -.5, las = 1)
points(ssa.cdf ~ times, type = "l", col = "black", lwd = 2)
points(ssaplus1.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dashed")
points(ssaplus2.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dotted")
mtext(side = 3, c("(e)"), adj = 0)
legend("bottomright",  c("SSAn", "SSA+ (fluc. 1)", "SSA+ (fluc. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)

########## Random
# Parameters, logistic growth (random1)
param1 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = randenvfun1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Parameters, logistic growth (random2)
param2 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = randenvfun2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(logisticINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhppCDF(times[i], init, param1, logisticINT)
    ssaplus2.cdf[i] <- nhppCDF(times[i], init, param2, logisticINT)
}
# Plot
plot(NA, xlab = "Time", ylab = "",
     xlim = c(0, 5), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
axis(1, at = seq(0, 5, 1), tcl = -.5)
axis(2, at = seq(0, 1, .2), tcl = -.5, las = 1)
points(ssa.cdf ~ times, type = "l", col = "black", lwd = 2)
points(ssaplus1.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dashed")
points(ssaplus2.cdf ~ times, type = "l", col = "black", lwd = 2, lty = "dotted")
mtext(side = 3, c("(f)"), adj = 0)
legend("bottomright",  c("SSAn", "SSA+ (rand. 1)", "SSA+ (rand. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)

#dev.off()
