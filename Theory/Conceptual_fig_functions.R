
# Conceptual Model of Species Response to Environmental Drivers
run.model <- function(N1, N2, r1, r2, beta12, beta21, K1, K2, env1_sigma1, 
                      env1_sigma2, env2_sigma1, env2_sigma2, time, env1, env2) {
  
  for (t in 1:(time-1)) {
    N1[t+1] <- N1[t]*exp(r1*(1-N1[t]/K1-beta12*N2[t]/K2)+env1_sigma1*env1[t]+env2_sigma1*env2[t])
    N2[t+1] <- N2[t]*exp(r2*(1-N2[t]/K2-beta21*N1[t]/K1)+env1_sigma2*env1[t]+env2_sigma2*env2[t])
  }
  
  results <- matrix(NA, nrow=2, ncol=time)
  results[1,] <- N1
  results[2,] <- N2
  
  return(results)
}

run.model.one.driver <- function(N1, N2, r1, r2, beta12, beta21, K1, K2, env1_sigma1, 
                                 env1_sigma2, time, env1) {
  
  for (t in 1:(time-1)) {
    N1[t+1] <- N1[t]*exp(r1*(1-N1[t]/K1-beta12*N2[t]/K2)+env1_sigma1*env1[t])
    N2[t+1] <- N2[t]*exp(r2*(1-N2[t]/K2-beta21*N1[t]/K1)+env1_sigma2*env1[t])
  }
  
  results <- matrix(NA, nrow=2, ncol=time)
  results[1,] <- N1
  results[2,] <- N2
  
  return(results)
}

run.model.one.driver.dispersal <- function(N1P1, N2P1, N1P2, N2P2, r1p1, r2p1, r1p2, r2p2, beta12, beta21, 
                                           K1P1, K2P1, K1P2, K2P2, env1_sigma1, env1_sigma2, env2_sigma1, env2_sigma2,
                                           time, env1, env2, disp) {
  
  for (t in 1:(time-1)) {
    N1P1[t+1] <- N1P1[t]*exp(r1p1*(1-N1P1[t]/K1P1-beta12*N2P1[t]/K2P1)+env1_sigma1*env1[t]) + disp*(N1P2[t]-N1P1[t])
    N2P1[t+1] <- N2P1[t]*exp(r2p1*(1-N2P1[t]/K2P1-beta21*N1P1[t]/K1P1)+env1_sigma2*env1[t]) + disp*(N2P2[t]-N2P1[t])
    
    N1P2[t+1] <- N1P2[t]*exp(r1p2*(1-N1P2[t]/K1P2-beta12*N2P2[t]/K2P1)+env2_sigma1*env2[t]) + disp*(N1P1[t]-N1P2[t])
    N2P2[t+1] <- N2P2[t]*exp(r2p2*(1-N2P2[t]/K2P2-beta21*N1P2[t]/K1P1)+env2_sigma2*env2[t]) + disp*(N2P1[t]-N2P2[t])
    
    
  }
  
  results <- matrix(NA, nrow=4, ncol=time)
  results[1,] <- N1P1
  results[2,] <- N2P1
  results[3,] <- N1P2
  results[4,] <- N2P2
  
  return(results)
}

plot.model <- function(time, results, col1, col2, col3, graph_min, graph_max) {
  
  plot(x[50:time], results[1,50:time], type="l", xaxt="n", yaxt="n", ylim=c(graph_min, graph_max), ylab="", xlab="", col=col1)
  lines(x[50:time], results[2,50:time], col=col2)
  lines(x[50:time], results[1,50:time]+results[2,50:time], col=col3)
}

plot.model.climate <- function(time, results, col1, col2, col3, graph_min, graph_max) {
  
  plot(x[1:time], results[1,1:time], type="l", xaxt="n", yaxt="n", ylim=c(graph_min, graph_max), ylab="", xlab="", col=col1)
  lines(x[1:time], results[2,1:time], col=col2)
  lines(x[1:time], results[1,1:time]+results[2,1:time], col=col3)
}

plot.model.spatial <- function(time, results, col1, col2, col3, graph_min, graph_max) {
  
  plot(seq(50:time), results[1,50:time], type="l", xaxt="n", yaxt="n", ylim=c(graph_min, graph_max), ylab="", xlab="", col=col1)
  lines(seq(50:time), results[2,50:time], col=col2)
  lines(seq(50:time), results[1,50:time]+results[2,50:time], col=col3)
}

run.model.growth <- function(N1, N2, N3, r1, r2, r3, K1, K2, K3, time) {
  
  for (t in 1:(time-1)) {
    N1[t+1] <- N1[t]*exp(r1*(1-N1[t]/K1))
    N2[t+1] <- N2[t]*exp(r2*(1-N2[t]/K2))
    N3[t+1] <- N3[t]*exp(r3*(1-N3[t]/K3))
  }
  
  results <- matrix(NA, nrow=3, ncol=time)
  results[1,] <- N1
  results[2,] <- N2
  results[3,] <- N3
  
  return(results)
}
