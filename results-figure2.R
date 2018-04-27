## Source nhpptimes function (CDF for non-stationary processes)
source("nhpptimes.R")
## View structure
nhpptimes

## Source constant environment
source("Environments/constant.R")
## View structure
constant
## Source increasing1 environment
source("Environments/increasing1.R")
## View structure
increasing1
## Source increasing2 environment
source("Environments/increasing2.R")
## View structure
increasing2
## Source fluctuating1 environment
source("Environments/fluctuating1.R")
## View structure
fluctuating1
## Source fluctuating2 environment
source("Environments/fluctuating2.R")
## View structure
fluctuating2

# Exponential growth
## Source environment-dependent function for exponential growth
source("birthdeath/Demographic functions/exponentialfun.R")
## View structure
exponentialfun
## Source intensity function
source("birthdeath/Intensity functions/exponentialINT.R")
## View structure
exponentialINT

# Logistic growth
## Source environment-dependent function for logistic growth
source("Demographic functions/logisticfun.R")
## View structure
logisticfun
## Source intensity function
source("Intensity functions/logisticINT.R")
## View structure
logisticINT


pdf(file="Figures/fig2.pdf", width = 7, height = 7)
par(mfrow = c(2, 3), mar = c(5, 6, 2, 1), lend = 2, cex.axis = 1.3, cex.lab = 1.5)
# Exponential initial values
init <- 100
# Exponential growth (increasing1)
param1 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = increasing1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Exponential growth (increasing2)
param2 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = increasing2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Create state change array
pproc <- matrix(c(1))
## Times
times <- seq(0, 10, .1)
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(exponentialINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhpptimes(times[i], init, param1, exponentialINT)
    ssaplus2.cdf[i] <- nhpptimes(times[i], init, param2, exponentialINT)
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
# Exponential growth (fluctuating1)
param1 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = fluctuating1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Exponential growth (fluctuating2)
param2 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = fluctuating2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Create state change array
pproc <- matrix(c(1))
## Times
times <- seq(0, 10, .1)
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(exponentialINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhpptimes(times[i], init, param1, exponentialINT)
    ssaplus2.cdf[i] <- nhpptimes(times[i], init, param2, exponentialINT)
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
# Exponential growth (irregular 1)
## Create irregular environmental values
set.seed(20170731)
randenv1 <- rnorm(length(times), mean = 1.5, sd = .25)
## View values
randenv1
# Exponential growth (irregular 2)
## Create function that interpolates these values
randenvfun1 <- splinefun(times, randenv1)
set.seed(20170731)
randenv2 <- rnorm(length(times), mean = 1.5, sd = .5)
## View values
randenv2
## Create function that interpolates these values
randenvfun2 <- splinefun(times, randenv2)
# Exponential growth (increasing1)
param1 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = randenvfun1)
param1 <- unlist(param1) # collapse to one-dimensional list
# Exponential growth (increasing2)
param2 <- list(b = .003,
               d = .0027,
               exponentialfun,
               env = randenvfun2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Create state change array
pproc <- matrix(c(1))
## Times
times <- seq(0, 10, .1)
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(exponentialINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhpptimes(times[i], init, param1, exponentialINT)
    ssaplus2.cdf[i] <- nhpptimes(times[i], init, param2, exponentialINT)
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
legend("bottomright",  c("SSAn", "SSA+ (irreg. 1)", "SSA+ (irreg. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)
## # Logistic growth (constant)
## param <- list(b1 = .003,
##               d1 = .0027,
##               d2 = .000003,
##               logisticfun,
##               env = constant)
## param <- unlist(param) # collapse to one-dimensional list
## ## Create state change array
## pproc <- matrix(c(1,
##                   -1),
##                 nrow = 2)
## ## Times
## times <- seq(0, 30, .1)
## ## Make CDFs
## ssa.cdf <- c()
## ssaplus.cdf <- c()
## for(i in 1:length(times)){
##     ssa.cdf[i] <- 1 - exp(-sum(logisticINT(0, 25, param)) * times[i])
##     ssaplus.cdf[i] <- nhpptimes(times[i], 25, param, logisticINT)
## }
## plot(NA, xlab = "Time", ylab = "Probability event has occurred",
##      xlim = c(0, 30), ylim = c(0, 1),
##      xaxt = "n", yaxt = "n", bty = "n",
##      main = "")
## axis(1, at = seq(0, 30, 5), tcl = -.25)
## axis(2, at = seq(0, 1, .2), tcl = -.25, las = 1)
## points(ssa.cdf ~ times, type = "l", col = "black", lwd = 2)
## points(ssaplus.cdf ~ times, type = "l", col = "black", lwd = 4, lty = "dashed")
## legend("bottomright",  c("SSA (constant)", "SSA+ (constant)"), bty = "n", inset = .05, adj = c(0, 0), lty = c("solid", "dashed"), lwd = 2)
# Initial values
init <- 100
# Logistic growth (increasing)
param1 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = increasing1)
param1 <- unlist(param1) # collapse to one-dimensional list
param2 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = increasing2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Create state change array
pproc <- matrix(c(1,
                  -1),
                nrow = 2)
## Times
times <- seq(0, 5, .1)
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(logisticINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhpptimes(times[i], init, param1, logisticINT)
    ssaplus2.cdf[i] <- nhpptimes(times[i], init, param2, logisticINT)
}
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
# Logistic growth (fluctuating)
param1 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = fluctuating1)
param1 <- unlist(param1) # collapse to one-dimensional list
param2 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = fluctuating2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Create state change array
pproc <- matrix(c(1,
                  -1),
                nrow = 2)
## Times
times <- seq(0, 5, .1)
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(logisticINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhpptimes(times[i], init, param1, logisticINT)
    ssaplus2.cdf[i] <- nhpptimes(times[i], init, param2, logisticINT)
}
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
# Exponential growth (irregular 1)
## Create irregular environmental values
set.seed(20170731)
randenv1 <- rnorm(length(times), mean = 1.5, sd = .25)
## View values
randenv1
# Exponential growth (irregular 2)
## Create function that interpolates these values
randenvfun1 <- splinefun(times, randenv1)
set.seed(20170731)
randenv2 <- rnorm(length(times), mean = 1.5, sd = .5)
## View values
randenv2
## Create function that interpolates these values
randenvfun2 <- splinefun(times, randenv2)
param1 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = randenvfun1)
param1 <- unlist(param1) # collapse to one-dimensional list
param2 <- list(b1 = .003,
              d1 = .0027,
              d2 = .000003,
              logisticfun,
              env = randenvfun2)
param2 <- unlist(param2) # collapse to one-dimensional list
## Create state change array
pproc <- matrix(c(1,
                  -1),
                nrow = 2)
## Times
times <- seq(0, 5, .1)
## Make CDFs
ssa.cdf <- c()
ssaplus1.cdf <- c()
ssaplus2.cdf <- c()
for(i in 1:length(times)){
    ssa.cdf[i] <- 1 - exp(-sum(logisticINT(0, init, param1)) * times[i])
    ssaplus1.cdf[i] <- nhpptimes(times[i], init, param1, logisticINT)
    ssaplus2.cdf[i] <- nhpptimes(times[i], init, param2, logisticINT)
}
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
legend("bottomright",  c("SSAn", "SSA+ (irreg. 1)", "SSA+ (irreg. 2)"), bty = "n", inset = .02, adj = c(0, 0), lty = c("solid", "dashed", "dotted"), lwd = 2)
dev.off()
