# Case study 1
# Comparing response to 2 different drivers operating on different timescales
# Code to accompany Figure 2 in the manuscript
# The long and the short of it: 
# Decomposing synchrony and compensation across temporal scales

rm(list=ls())
library(tsvr) # library for the timescale specific variance ratio
library(here)

source(here("/Theory/Conceptual_fig_functions.R"))

time <- 150        # run model for 150 timesteps; use timesteps 50-150 for calculations
time_graph <- 150  # use timesteps 50-100 for visualizations

graph_min <- 450   # minimum and maximum abundances for graphing
graph_max <- 1550

# environmental drivers
x_fine <- seq(from=50, to=time_graph, by=0.1) # to make a smoother graph
x <- seq(from=1, to=time)                     # for our actual calculations
A1 <- .5
B1 <- 2*pi/3
C1 <- 2
D1 <- 0
env1 <- A1*sin(B1*x+C1)+D1
env1_fine <- A1*sin(B1*x_fine+C1)+D1

A2 <- 1
B2 <- 2*pi/20
C2 <- 0
D2 <- 0
env2 <- A2*sin(B2*x+C2)+D2
env2_fine <- A2*sin(B2*x_fine+C2)+D2

# graph it
#quartz(width=9, height=4)
pdf(here("Figures/Case_study_1.pdf"), width=9.5, height=4)
par(mfrow=c(2,4), mar=c(1,1,1,.5), oma=c(2,2,2,.5))
plot(x_fine, env1_fine, type="l", ylim=c(-1.5, 1.5), ylab="Environment", xlab="", xaxt="n", yaxt="n", col="blue")
text(x=50.25, y= 1.42, "a)", cex=1.25)
mtext("Driver", side=2, outer=FALSE, line=1)
mtext("Short Timescale", side=3, outer=FALSE, line=1)

plot(x_fine, env2_fine, type="l", ylim=c(-1.5, 1.5), ylab="Environment", xlab="", xaxt="n", yaxt="n", col="darkgreen")
text(x=50.25, y= 1.42, "b)", cex=1.25)
mtext("Long Timescale", side=3, outer=FALSE, line=1)

plot(x_fine, env1_fine, type="l", ylim=c(-1.5, 1.5), ylab="Environment", xlab="", xaxt="n", yaxt="n", col="blue")
lines(x_fine, env2_fine, col="darkgreen")
text(x=50.25, y= 1.42, "c)", cex=1.25)
mtext("Composite", side=3, outer=FALSE, line=1)

# legend
plot.new()
legend("center", c("Short-term Driver", "Long-term Driver", "Species 1", "Species 2", "Total Biomass"), 
       col=c("blue", "darkgreen", "black", "darkgrey", "red"), lwd=2, bty="n")

# Response to driver 1 only
# Run model -- synchronous response to driver 1
# Parameters
r1 <- .5
r2 <- .5
K1 <- 1000
K2 <- 1100
env1_sigma1 <- 0.2
env1_sigma2 <- 0.2
env2_sigma1 <- 0
env2_sigma2 <- 0
beta12 <- 0.5
beta21 <- 0.5

# starting conditions for species 1 (N1) and species 2 (N2)
N1 <- N2 <- rep(NA, time)
N1[1] <- K1
N2[1] <- K2

# run dynamical model
results1 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, env1, env2)

plot.model(time_graph, results1, "black", "darkgrey", "red", graph_min, graph_max)
text(x=50.25, y= graph_max-30, "d)", cex=1.25)
mtext("Abundance", side=2, outer=FALSE, line=1)

# Response to driver 2 ONLY
# Run model -- compensatory response to driver 2
# Parameters
r1 <- .5
r2 <- .5
K1 <- 1000
K2 <- 1100
env1_sigma1 <- 0
env1_sigma2 <- 0
env2_sigma1 <- 0.1
env2_sigma2 <- -0.1
beta12 <- 0.5
beta21 <- 0.5

# starting conditions for species 1 (N1) and species 2 (N2)
N1 <- N2 <- rep(NA, time)
N1[1] <- K1
N2[1] <- K2

results2 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, env1, env2)

plot.model(time_graph, results2, "black", "darkgrey", "red", graph_min, graph_max)
text(x=50.25, y= graph_max-30, "e)", cex=1.25)
mtext("Time", side=1, outer=FALSE, line=1)

# Run model -- response to both drivers now
# Parameters
r1 <- .5
r2 <- .5
K1 <- 1000
K2 <- 1100
env1_sigma1 <- 0.2
env1_sigma2 <- 0.2
env2_sigma1 <- 0.1
env2_sigma2 <- -0.1
beta12 <- 0.5
beta21 <- 0.5

# starting conditions
N1 <- N2 <- rep(NA, time)
N1[1] <- K1
N2[1] <- K2

results3 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, env1, env2)

plot.model(time_graph, results3, "black", "darkgrey", "red", graph_min, graph_max)
text(x=50.25, y= graph_max-30, "f)", cex=1.25)

# Calculate variance ratio for response to both drivers

vr.class <- vector("list", 1)
vr.trial.short <- vector("list", 1)
vr.trial.long <- vector("list", 1)

vr.trial <- tsvreq_classic(results3[,50:time])
aggresShort <- aggts(vr.trial, vr.trial$ts[vr.trial$ts<4])
aggresLong <- aggts(vr.trial, vr.trial$ts[vr.trial$ts>=4])
vr.trial.short[[1]] <- aggresShort[[3]]
vr.trial.long[[1]] <- aggresLong[[3]]
temp <-  vreq_classic(results3[,50:time])
vr.class[[1]] <- temp[[3]]

# Bar graph of variance ratios
plot(c(1,2,3), c(vr.trial.short[[1]], vr.trial.long[[1]], vr.class[[1]]), col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=FALSE, cex=1.25)
axis(side=2, at=c(0,1,2), labels=TRUE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
mtext("Variance Ratio", side=3, outer=FALSE, line=1)
abline(h=1, col="grey", lty=2, lwd=2)
text(x=0.72, y= 1.98, "g)", cex=1.25)

dev.off()
