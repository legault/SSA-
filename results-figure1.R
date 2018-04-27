# Standard time sequence
times <- seq(0, 10, 1)
# Create time sequences for smooth function
times2 <- seq(0, 10, .01)

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
# Irregular1
## Create irregular environmental values
set.seed(20170731)
randenv1 <- rnorm(length(times), mean = 1.5, sd = .25)
randenv1 <- splinefun(times, randenv1)
# Irregular2
set.seed(20170731)
randenv2 <- rnorm(length(times), mean = 1.5, sd = .5)
randenv2 <- splinefun(times, randenv2)

pdf(file="Figures/fig1.pdf", height = 4, width = 8)
par(mfrow = c(1, 3), mar = c(5, 6, 2, 1), lend = 2, cex.axis = 1.3, cex.lab = 1.5)
# Plots
## Increasing
plot(NA, xlab = expression(paste("Time (", italic(t), ")")), ylab = "",
     xlim = c(0, 10), ylim = c(0.5, 2.5),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
title(ylab = expression(paste("Environment at ", italic("t"))),, line = 4)
axis(1, at = seq(0, 10, 2), tcl = -.5, cex.lab= 1.4)
axis(2, at = seq(0.5, 2.5, .5), tcl = -.5, las = 1, cex = 1.4)
points(increasing1(times) ~ times, type = "l", col = "black", lwd = 3, lty = "dashed")
points(increasing2(times) ~ times, type = "l", col = "black", lwd = 2, lty = "dotted")
legend("topright", c("increasing 1", "increasing 2"), bty = "n", inset = .02, adj = c(0, 0), lty = c("dashed", "dotted"), lwd = c(3, 2), cex = 1.3)
mtext(side = 3, c("(a)"), adj = 0)
## Fluctuating
plot(NA, xlab = expression(paste("Time (", italic(t), ")")), ylab = "",
     xlim = c(0, 10), ylim = c(0, 2),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
axis(1, at = seq(0, 10, 2), tcl = -.5, cex = 1.2)
axis(2, at = seq(0, 2, .5), tcl = -.5, las = 1)
points(fluctuating1(times2) ~ times2, type = "l", col = "black", lwd = 3, lty = "dashed")
points(fluctuating2(times2) ~ times2, type = "l", col = "black", lwd = 2, lty = "dotted")
legend("topright", c("fluctuating 1", "fluctuating 2"), bty = "n", inset = .02, adj = c(0, 0), lty = c("dashed", "dotted"), lwd = c(3, 2), cex = 1.3)
mtext(side = 3, c("(b)"), adj = 0)
## Irregular
plot(NA, xlab = expression(paste("Time (", italic(t), ")")), ylab = "",
     xlim = c(0, 10), ylim = c(0, 4),
     xaxt = "n", yaxt = "n", bty = "n",
     main = "")
axis(1, at = seq(0, 10, 2), tcl = -.5, cex = 1.2)
axis(2, at = seq(0, 4, 1), tcl = -.5, las = 1)
points(randenv1(times2) ~ times2, type = "l", col = "black", lwd = 3, lty = "dashed")
points(randenv2(times2) ~ times2, type = "l", col = "black", lwd = 2, lty = "dotted")
legend("topright", c("irregular 1", "irregular 2"), bty = "n", inset = .02, adj = c(0, 0), lty = c("dashed", "dotted"), lwd = c(3, 2), cex = 1.3)
mtext(side = 3, c("(c)"), adj = 0)
dev.off()
