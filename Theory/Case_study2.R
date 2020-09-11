# Case study 2
# Comparing growth forms/lag effects on synchrony dynamics and timescales of synchrony
# Code to accompany Figure 4 in the manuscript
# The long and the short of it: 
# Decomposing synchrony and compensation across temporal scales
# Conceptual figure 2 

rm(list=ls())
library(tsvr)
library(here)

source(here("/Theory/Conceptual_fig_functions.R"))

time <- 150        # run model for 150 timesteps; use timesteps 50-150 for calculations
time_graph <- 150  # use timesteps 50-100 for visualizations

graph_min <- 350
graph_max <- 2000

# environmental drivers
x_fine <- seq(from=50, to=time_graph, by=0.1)
x <- seq(from=1, to=time)
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

# Run model -- no environmental driver
# Parameters
r1 <- .15
r2 <- 1
r3 <- 1.8
K1 <- 1000
K2 <- 1000
K3 <- 1000

N1 <- N2 <- N3 <- rep(NA, time)
N1[1] <- 500
N2[1] <- 500
N3[1] <- 500

# run each species seperately to visualize growth functions
results_base <- run.model.growth(N1, N2, N3, r1, r2, r3, K1, K2, K3, time)

# set up layout plot
m2 <- matrix(c(1, 1, 1, 6, 6, 6, 11, 11,
               1, 1, 1, 6, 6, 6, 11, 11,
               1, 1, 1, 6, 6, 6, 11, 11,
               1, 1, 1, 6, 6, 6, 11, 11,
               2, 2, 2, 7, 7, 7, 11, 11,
               3, 3, 3, 8, 8, 8, 11, 11,
               3, 3, 3, 8, 8, 8, 11, 11,
               3, 3, 3, 8, 8, 8, 11, 11,
               3, 3, 3, 8, 8, 8, 11, 11,
               4, 4, 4, 9, 9, 9, 11, 11,
               4, 4, 4, 9, 9, 9, 11, 11,
               4, 4, 4, 9, 9, 9, 11, 11,
               4, 4, 4, 9, 9, 9, 11, 11,
               5, 5, 5, 10, 10, 10, 11, 11,
               5, 5, 5, 10, 10, 10, 11, 11,
               5, 5, 5, 10, 10, 10, 11, 11,
               5, 5, 5, 10, 10, 10, 11, 11), nrow=8)

pdf(here("Figures/Case_study_2.pdf"), width=10, height=4)
#quartz(width=10, height=4)
par(mar=c(.5,.5,.5,.5), oma=c(2,.5,2,.5))
layout(m2)
plot(x[1:30], results_base[1,1:30], type="l", xaxt="n", yaxt="n", 
     ylim=c(graph_min, 1350), ylab="", xlab="", col="orange3", lwd=1)
lines(x[1:30], results_base[2,1:30], col="black", lwd=1)
lines(x[1:30], results_base[3,1:30], col="darkorchid1", lwd=1)
lines(x[1:30], results_base[2,1:30], col="black", lwd=1)
text(x=1.25, y= 1329, "a)", cex=1.25)
mtext("Growth Functions", side=3, outer=FALSE, line=.5)

plot.new()

# Run model -- case 1 and 1: both species have the same growth rates
# Parameters
r1 <- 1
r2 <- 1
K1 <- 1100
K2 <- 1000
env1_sigma1 <- 0.1
env1_sigma2 <- 0.1
env2_sigma1 <- 0.2
env2_sigma2 <- 0.2
beta12 <- 0.5
beta21 <- 0.5

N1 <- N2 <- rep(NA, time)
N1[1] <- K1
N2[1] <- K2

results1 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, env1, env2)

plot.model(time_graph, results1, "black", "white", "red", graph_min, graph_max)
mtext("Abundance", side=2, outer=FALSE, line=2)
lines(x[50:time], results1[2,50:time], col="black", lty=2)
text(x=50.25, y= graph_max-30, "b)", cex=1.25)

# Run model -- case 1 and 2: one species has a lagged response to environmental change
# Parameters
r1 <- 1
r2 <- .15
K1 <- 1100
K2 <- 1000
env1_sigma1 <- 0.1
env1_sigma2 <- 0.1
env2_sigma1 <- 0.2
env2_sigma2 <- 0.2
beta12 <- 0.5
beta21 <- 0.5

N1 <- N2 <- rep(NA, time)
N1[1] <- K1
N2[1] <- K2

results2 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, env1, env2)

plot.model(time_graph, results2, "black", "orange3", "red", graph_min, graph_max)
text(x=50.25, y= graph_max-30, "c)", cex=1.25)

# Run model -- case 1 and 3: one species has dampening oscillatory dynamics 
# Parameters
r1 <- 1
r2 <- 1.8
K1 <- 1100
K2 <- 1000
env1_sigma1 <- 0.1
env1_sigma2 <- 0.1
env2_sigma1 <- 0.2
env2_sigma2 <- 0.2
beta12 <- 0.5
beta21 <- 0.5

N1 <- N2 <- rep(NA, time)
N1[1] <- K1
N2[1] <- K2

results3 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, env1, env2)

plot.model(time_graph, results3, "black", "darkorchid1", "red", graph_min, graph_max)
text(x=50.25, y= graph_max-30, "d)", cex=1.25)

# Plot Environmental Driver
plot(x_fine, env1_fine, type="l", ylim=c(-1.5, 1.5), ylab="Environment", xlab="", xaxt="n", yaxt="n", col="blue", lwd=1)
lines(x_fine, env2_fine, col="darkgreen", lwd=1)
text(x=50.25, y= 1.42, "e)", cex=1.25)
mtext("Driver", side=1, outer=FALSE, line=1)

plot.new()

# Calculate variance ratio for all of the scenarios

# Combine results
compiled <- list(results1[,50:time], results2[,50:time], results3[,50:time])

vr.class <- vector("list", 3)
vr.trial.short <- vector("list", 3)
vr.trial.long <- vector("list", 3)

for(i in 1:3){
  vr.trial <- tsvreq_classic(compiled[[i]])
  aggresShort <- aggts(vr.trial, vr.trial$ts[vr.trial$ts<4])
  aggresLong <- aggts(vr.trial, vr.trial$ts[vr.trial$ts>=4])
  vr.trial.short[[i]] <- aggresShort[[3]]
  vr.trial.long[[i]] <- aggresLong[[3]]
  temp <-  vreq_classic(compiled[[i]])
  vr.class[[i]] <- temp[[3]]
}

# Bar graph of variance ratios
plot(c(1,2,3), c(vr.trial.short[[1]], vr.trial.long[[1]], vr.class[[1]]), col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=TRUE, cex=1.25, tck=-.03, padj=.25)
axis(side=2, at=c(0,1,2), labels=TRUE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
#mtext("Variance Ratio", side=3, outer=FALSE, line=1)
abline(h=1, col="grey", lty=2, lwd=2)
mtext("Variance Ratio", side=2, outer=FALSE, line=2)
text(x=0.72, y= 1.98, "f)", cex=1.25)

plot(c(1,2,3), c(vr.trial.short[[2]], vr.trial.long[[2]], vr.class[[2]]), col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=TRUE, cex=1.25, tck=-.03, padj=.25)
axis(side=2, at=c(0,1,2), labels=FALSE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
abline(h=1, col="grey", lty=2, lwd=2)
text(x=0.72, y= 1.98, "g)", cex=1.25)

plot(c(1,2,3), c(vr.trial.short[[3]], vr.trial.long[[3]], vr.class[[3]]), col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=TRUE, cex=1.25, tck=-.03, padj=.25)
axis(side=2, at=c(0,1,2), labels=FALSE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
abline(h=1, col="grey", lty=2, lwd=2)
text(x=0.72, y= 1.98, "h)", cex=1.25)

plot.new()
legend_order <- matrix(1:8, ncol=4, byrow = TRUE)


legend("bottomleft", c("Short-term Driver ", "Long-term Driver ","Species 1 (r = 1.00)", "Species 2 (r = 0.15)", 
                       "Species 3 (r = 1.8)", "Total Biomass"), 
       col=c("blue", "darkgreen",  "black", "orange3", "darkorchid1", "red"), 
       lwd=2, bty="n", cex=1, horiz=TRUE)

dev.off()