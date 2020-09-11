# Case study 4
# Comparing synchrony with global change altering the environmental driver
# Code to accompany Figure 6 in the manuscript
# The long and the short of it: 
# Decomposing synchrony and compensation across temporal scales

rm(list=ls())
library(tsvr) # library for the timescale specific variance ratio
library(here)

source(here("/Theory/Conceptual_fig_functions.R"))

time_graph <- 100
time <- 100

# environmental drivers
x_fine <- seq(from=1, to=time_graph, by=0.1)
x <- seq(from=1, to=time)
A1 <- .25
B1 <- 2*pi/3
C1 <- 2
D1 <- 0
env1 <- A1*sin(B1*x+C1)+D1
env1_fine <- A1*sin(B1*x_fine+C1)+D1

A2 <- .5
B2 <- 2*pi/20
C2 <- 0
D2 <- 0
env2 <- A2*sin(B2*x+C2)+D2
env2_fine <- A2*sin(B2*x_fine+C2)+D2

# graph the driver
#quartz(width=9, height=6)
pdf(here("Figures/Case_study_4.pdf"), width=9.5, height=5)
par(mfrow=c(3,4), mar=c(1,1,1,.5), oma=c(2,2,2,.5))

plot(x_fine, env1_fine, type="l", ylim=c(-1.5, 1.5), ylab="Environment", xlab="", xaxt="n", yaxt="n", col="blue")
lines(x_fine, env2_fine, col="darkgreen")
text(x=0.25, y= 1.42, "a)", cex=1.25)
mtext("Individual Drivers", side=3, outer=FALSE, line=1)

# composite effect of both drivers --------------------------------------------------------------------------------------------
threshold  <- .5
env_sum <- env1_fine + env2_fine
env_sp2 <- ifelse(env_sum < threshold, NA, env_sum) # if the environment is below the threshold, 
                                                    # species 2 does not respond

# plot it
plot(x_fine, env_sum, type="l", ylim=c(-1.5, 1.5), ylab="", xlab="", xaxt="n", yaxt="n", col="darkorchid")
lines(x_fine, env_sp2, col="palevioletred", lwd=1.5)
text(x=0.25, y= 1.42, "b)", cex=1.25)
abline(h=threshold, col="grey", lty=2)
mtext("Stable Climate", side=3, outer=FALSE, line=1)

# changing climate with both drivers --------------------------------------------------------------------------------------------
temp <- (.5/time_graph)*x_fine
env_sum_change <- temp + env_sum
env_sp2_change <- ifelse(env_sum_change < threshold, NA, env_sum_change)

# plot it
plot(x_fine, env_sum_change, type="l", ylim=c(-1.5, 1.5), ylab="", xlab="", xaxt="n", yaxt="n", col="darkorchid")
lines(x_fine, env_sp2_change, col="palevioletred", lwd=1.5)
text(x=0.25, y= 1.42, "c)", cex=1.25)
abline(h=threshold, col="grey", lty=2)
mtext("Changing Climate", side=3, outer=FALSE, line=1)

# new climate with both drivers --------------------------------------------------------------------------------------------
env_sum_new <- .5 + env_sum
env_sp2_new <- ifelse(env_sum_new < threshold, NA, env_sum_new)

# plot it
plot(x_fine, env_sum_new, type="l", ylim=c(-1.5, 1.5), ylab="", xlab="", xaxt="n", yaxt="n", col="darkorchid")
lines(x_fine, env_sp2_new, col="palevioletred", lwd=1.5)
text(x=0.25, y= 1.42, "d)", cex=1.25)
abline(h=threshold, col="grey", lty=2)
mtext("New Stable Climate", side=3, outer=FALSE, line=1)

# Run abundances for each species --------------------------------------------------------------------------------------------
# Pull out annual timeseries
vals <- which (x_fine %% 1 == 0)

# Parameters
r1 <- .65
r2 <- .65
K1 <- 1100
K2 <- 1000
env1_sigma1 <- 0.5
env1_sigma2 <- 0
env2_sigma1 <- 0
env2_sigma2 <- .5
beta12 <- 0.5
beta21 <- 0.5

N1 <- N2 <- rep(NA, time)
N1[1] <- K1
N2[1] <- K2

# baseline
s1 <- env_sum[vals]
s2 <- env_sp2[vals]
s2 <- ifelse(is.na(s2), 0, s2) 

results1 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, s1, s2)

# transient case
N1 <- N2 <- rep(NA, time)
N1[1] <- results1[1,time]
N2[1] <- results1[2,time]

s1 <- env_sum_change[vals]
s2 <- env_sp2_change[vals]
s2 <- ifelse(is.na(s2), 0, s2) 
results2 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, s1, s2)

# new baseline
N1 <- N2 <- rep(NA, time)
N1[1] <- results2[1,time]
N2[1] <-results2[2,time]

s1 <- env_sum_new[vals]
s2 <- env_sp2_new[vals]
s2 <- ifelse(is.na(s2), 0, s2) 
results3 <- run.model(N1, N2, r1, r2, beta12, beta21, K1, K2, 
                      env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2, time, s1, s2)

# graph abundance results ----------------------------------------------------------------------------------------------
graph_min <- 350
graph_max <- 2750

plot.new()

plot.model.climate(time, results1, "black", "darkgrey", "red", graph_min, graph_max)
text(x=2.25, y= graph_max-50, "e)", cex=1.25)
mtext("Abundance", side=2, outer=FALSE, line=1)

plot.model.climate(time, results2, "black", "darkgrey", "red", graph_min, graph_max)
text(x=2.25, y= graph_max-50, "f)", cex=1.25)

plot.model.climate(time, results3, "black", "darkgrey", "red", graph_min, graph_max)
text(x=2.25, y= graph_max-50, "g)", cex=1.25)

# Calculate variance ratio ----------------------------------------------------------------------------------------------

# Combine results
compiled <- list(results1[,1:time], results2[,1:time], results3[,1:time])

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

plot.new()
legend("bottomleft", c("Short-term Driver", "Long-term Driver", "", "Species 1 Response", 
                       "Species 2 Response", "", "Species 1 Biomass", "Species 2 Biomass", "Total Biomass"),
       col=c("blue", "darkgreen", NA, "darkorchid", "palevioletred", NA, "black", "darkgrey", "red"), lwd=2, bty="n")

plot(c(1,2,3), c(vr.trial.short[[1]], vr.trial.long[[1]], vr.class[[1]]), col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xlab="", xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=TRUE, cex=1.25, tck=-.03, padj=.25)
axis(side=2, at=c(0,1,2), labels=TRUE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
#mtext("Variance Ratio", side=3, outer=FALSE, line=1)
#abline(h=1, col="grey", lty=2, lwd=2)
lines(c(0.7, 3.3), c(1,1), lty=2, lwd=2, col="grey")
mtext("Variance Ratio", side=2, outer=FALSE, line=1)
text(x=0.72, y= 1.98, "h)", cex=1.25)

plot(c(1,2,3), c(vr.trial.short[[2]], vr.trial.long[[2]], vr.class[[2]]), ylab="", xlab="", col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=TRUE, cex=1.25, tck=-.03, padj=.25)
axis(side=2, at=c(0,1,2), labels=FALSE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
#abline(h=1, col="grey", lty=2, lwd=2)
lines(c(0.7, 3.3), c(1,1), lty=2, lwd=2, col="grey")
text(x=0.72, y= 1.98, "i)", cex=1.25)

plot(c(1,2,3), c(vr.trial.short[[3]], vr.trial.long[[3]], vr.class[[3]]), ylab="", xlab="", col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=TRUE, cex=1.25, tck=-.03, padj=.25)
axis(side=2, at=c(0,1,2), labels=FALSE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
#abline(h=1, col="grey", lty=2, lwd=2)
lines(c(0.7, 3.3), c(1,1), lty=2, lwd=2, col="grey")
text(x=0.72, y= 1.98, "j)", cex=1.25)

dev.off()
