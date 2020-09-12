# Case study 3
# Comparing two patches in a metacommunity, each with different environmental drivers. 
# Comparison of two different dispersal rates 
# Code to accompany Figure 5 in the manuscript
# The long and the short of it: 
# Decomposing synchrony and compensation across temporal scales

rm(list=ls())
library(tsvr) # library for the timescale specific variance ratio
library(here)

source(here("/Theory/Conceptual_fig_functions.R"))

time <- 100        # run model for 100 timesteps; use timesteps 50-100 for calculations
time_graph <- 100  # use timesteps 50-100 for visualizations

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

# graph the driver
#quartz(width=8, height=6)
pdf(here("Figures/Case_study_3.pdf"), width=8, height=6)
par(mfcol=c(4,3), mar=c(.5,2,.5,.5), oma=c(2,2,2,.5))


plot(x_fine, env1_fine, type="l", ylim=c(-1.5, 1.5), ylab="", xlab="", xaxt="n", yaxt="n", col="blue")
text(x=50.25, y= 1.42, "a)", cex=1.25)
mtext("Patch 1", side=2, outer=FALSE, line=.5)
mtext("Drivers", side=3, outer=FALSE, line=.5)

plot(x_fine, env2_fine, type="l", ylim=c(-1.5, 1.5), ylab="", xlab="", xaxt="n", yaxt="n", col="darkgreen")
text(x=50.25, y= 1.42, "d)", cex=1.25)
mtext("Patch 2", side=2, outer=FALSE, line=.5)

plot.new()
plot.new()
legend("topleft", c("Short-term Driver", "Long-term Driver", "", "Species 1 Biomass", "Species 2 Biomass", "Total Biomass"),
       col=c("blue", "darkgreen", NA, "black", "darkgrey", "red"), lwd=2, bty="n")

# Run abundances for each species --------------------------------------------------------------------------------------------

# No dispersal scenario
# Parameters
# patches are denotes by p1 versus p2 for parameters
r1p1 <- .5
r2p1 <- .5
r1p2 <- .5
r2p2 <- .5
K1P1 <- 1100
K2P1 <- 1000
K1P2 <- 1100
K2P2 <- 1000
env1_sigma1 <- 0.5  # environmental signature in patch 1
env1_sigma2 <- 0.5
env2_sigma1 <- .1   # environmental signature in patch 2
env2_sigma2 <- -.1
beta12 <- 0.5
beta21 <- 0.5

disp <- 0

# Starting conditions
N1P1 <- N2P1 <- N1P2 <- N2P2 <- rep(NA, time)
N1P1[1] <- K1P1
N2P1[1] <- K2P1
N1P2[1] <- K1P2
N2P2[1] <- K2P2

results_no <- run.model.one.driver.dispersal(N1P1, N2P1, N1P2, N2P2, r1p1, r2p1, r1p2, r2p2, beta12, beta21, 
                                             K1P1, K2P1, K1P2, K2P2, env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2,
                                             time, env1, env2, disp)

patch1 <- matrix(c(results_no[1,], results_no[2,]), nrow=2, byrow=TRUE)
patch2 <-  matrix(c(results_no[3,], results_no[4,]), nrow=2, byrow=TRUE)
metacommunity <- matrix(c(results_no[1,]+results_no[3,], results_no[2,]+results_no[4,]), nrow=2, byrow=TRUE)

# graph abundances
graph_min <- 300
graph_max <- 3500

graph_min_land <- 500
graph_max_land <- 3500

plot.model.spatial(time_graph, patch1, "black", "darkgrey", "red", graph_min, graph_max)
text(x=1.25, y= graph_max-95, "b)", cex=1.25)
mtext("No Dispersal", side=3, outer=FALSE, line=.5)
mtext("Abundance", side=2, outer=FALSE, line=.5)

plot.model.spatial(time_graph, patch2, "black", "darkgrey", "red", graph_min, graph_max)
text(x=1.25, y= graph_max-95, "e)", cex=1.25)
mtext("Abundance", side=2, outer=FALSE, line=.5)

plot.model.spatial(time_graph, metacommunity, "black", "darkgrey", "red", graph_min_land, graph_max_land)
text(x=1.25, y= graph_max_land-95, "g)", cex=1.25)
mtext("  Landscape \n Abundance", side=2, outer=FALSE, line=1)

# Calculate variance ratio for each patch and for the metacommunity
compiled <- list(patch1[,50:time], patch2[,50:time], metacommunity[,50:time])

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

# Plot variance ratio

# Bar graph of variance ratios
plot(c(2,3,4,8,9,10,14,15,16), c(vr.trial.short[[1]], vr.trial.long[[1]], vr.class[[1]], 
                                 vr.trial.short[[2]], vr.trial.long[[2]], vr.class[[2]], 
                                 vr.trial.short[[3]], vr.trial.long[[3]], vr.class[[3]] ), 
     col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", ylab="", yaxt="n", 
     ylim=c(-.05, 2.05), xlim=c(0,18))
axis(side=1, at=c(3,9,15), labels=c("Patch 1", "Patch 2", "Landscape"), tick=FALSE, cex=.75)
axis(side=2, at=c(0,1,2), labels=TRUE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
mtext("Variance Ratio", side=2, outer=FALSE, line=1.25)
abline(h=1, col="grey", lty=2, lwd=2)
abline(v=6, col="grey", lty=2, lwd=2)
abline(v=12, col="grey", lty=2, lwd=2)
text(x=0.25, y= 1.98, "i)", cex=1.25)

# Moderate Dispersal --------------------------------------------------------------------------------------------------------
disp <- 0.4

results_mod <- run.model.one.driver.dispersal(N1P1, N2P1, N1P2, N2P2, r1p1, r2p1, r1p2, r2p2, beta12, beta21, 
                                             K1P1, K2P1, K1P2, K2P2, env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2,
                                             time, env1, env2, disp)

patch1 <- matrix(c(results_mod[1,], results_mod[2,]), nrow=2, byrow=TRUE)
patch2 <-  matrix(c(results_mod[3,], results_mod[4,]), nrow=2, byrow=TRUE)
metacommunity <- matrix(c(results_mod[1,]+results_mod[3,], results_mod[2,]+results_mod[4,]), nrow=2, byrow=TRUE)

plot.model.spatial(time_graph, patch1, "black", "darkgrey", "red", graph_min, graph_max)
text(x=1.25, y= graph_max-95, "c)", cex=1.25)
mtext("With Dispersal", side=3, outer=FALSE, line=.5)
#mtext("Abundace", side=2, outer=FALSE, line=2.5, at=2200)

plot.model.spatial(time_graph, patch2, "black", "darkgrey", "red", graph_min, graph_max)
text(x=1.25, y= graph_max-95, "f)", cex=1.25)

plot.model.spatial(time_graph, metacommunity, "black", "darkgrey", "red", graph_min_land, graph_max_land)
text(x=1.25, y= graph_max_land-95, "h)", cex=1.25)

# Calculate variance ratio for each patch and for the metacommunity
compiled <- list(patch1[,50:time], patch2[,50:time], metacommunity[,50:time])

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


# Plot variance ratio

# Bar graph of variance ratios
plot(c(2,3,4,8,9,10,14,15,16), c(vr.trial.short[[1]], vr.trial.long[[1]], vr.class[[1]], 
                                 vr.trial.short[[2]], vr.trial.long[[2]], vr.class[[2]], 
                                 vr.trial.short[[3]], vr.trial.long[[3]], vr.class[[3]] ), 
     col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", ylab="", yaxt="n", 
     ylim=c(-.05, 2.05), xlim=c(0,18))
axis(side=1, at=c(3,9,15), labels=c("Patch 1", "Patch 2", "Landscape"), tick=FALSE, cex=.75)
axis(side=2, at=c(0,1,2), labels=TRUE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
abline(h=1, col="grey", lty=2, lwd=2)
abline(v=6, col="grey", lty=2, lwd=2)
abline(v=12, col="grey", lty=2, lwd=2)
text(x=0.25, y= 1.98, "j)", cex=1.25)

dev.off()

