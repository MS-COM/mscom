
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Turn off sci notation
options(scipen=999)

# Packages
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(optimx)

# Helper functions
################################################################################

# Plot catch time series
plot_catch <- function(C_mat, yrs){
  nstocks <- ncol(C_mat)
  cmax <- max(C_mat)
  colors <- brewer.pal(nstocks, "Set1")
  plot(x=yrs, y=C_mat[,1], bty="n", type="l", las=1, col=colors[1],
       ylim=c(0,cmax), xlab="", ylab="Catch", main="Catch") 
  for(i in 2:nstocks){
    lines(x=yrs, y=C_mat[,i], col=colors[i])
  }
  legend("topleft", bty="n", legend=paste("Stock", 1:nstocks), col=colors, lty=1)
}

# Plot biomass time series
plot_biomass <- function(B_mat, yrs){
  nstocks <- ncol(C_mat)
  bmax <- max(B_mat)
  colors <- brewer.pal(nstocks, "Set1")
  plot(x=yrs, y=B_mat[,1], bty="n", type="l", las=1, col=colors[1],
       ylim=c(0,bmax), xlab="", ylab="Biomass", main="Biomass") 
  for(i in 2:nstocks){
    lines(x=yrs, y=B_mat[,i], col=colors[i])
  }
  legend("topleft", bty="n", legend=paste("Stock", 1:nstocks), col=colors, lty=1)
}

# Plot fishing mortality time series
plot_f <- function(F_mat, yrs){
  nstocks <- ncol(F_mat)
  fmax <- max(F_mat)
  colors <- brewer.pal(nstocks, "Set1")
  plot(x=yrs, y=F_mat[,1], bty="n", type="l", las=1, col=colors[1],
       ylim=c(0,fmax), xlab="", ylab="Catch / Biomass", main="Fishing mortality") 
  for(i in 2:nstocks){
    lines(x=yrs, y=F_mat[,i], col=colors[i])
  }
  legend("topleft", bty="n", legend=paste("Stock", 1:nstocks), col=colors, lty=1)
}

# Categorical status
bbmsy2status <- function(bbmsy){
  as.character(cut(bbmsy, breaks=c(0,0.5,1.5,999), labels=c("over", "fully", "under")))
}

# Calculate r prior
calc_r_priors <- function(res){
  for(i in 1:length(res)){
    if(res[i]=="High"){r_prior <- c(0.6,1.5)}
    if(res[i]=="Medium"){r_prior <- c(0.2,0.8)}
    if(res[i]=="Low"){r_prior <- c(0.05,0.5)}
    if(res[i]=="Very low"){r_prior <- c(0.015,0.1)}
    if(i==1){r_priors <- r_prior}else{r_priors <- rbind(r_priors, r_prior)}
  }
  colnames(r_priors) <- c("r_lo", "r_hi")
  rownames(r_priors) <- NULL
  return(r_priors)
}

# Calculate k prior
calc_k_priors <- function(C_mat, r_priors){
  for(i in 1:ncol(C_mat)){
    catch <- C_mat[,i]
    r_prior <- r_priors[i,]
    r_lo <- r_prior[1]
    r_hi <- r_prior[2]
    c_max <- max(catch)
    k_lo <- c_max / r_hi
    k_hi <- 12 * c_max / r_lo
    k_prior <- c(k_lo, k_hi)
    if(i==1){k_priors <- k_prior}else{k_priors <- rbind(k_priors, k_prior)}
  }
  colnames(k_priors) <- c("k_lo", "k_hi")
  rownames(k_priors) <- NULL
  return(k_priors)
}

# Simulate data
################################################################################

# Params
yrs <- 1970:2010
nyrs <- length(yrs)
nstocks <- 3
stocks <- paste("Stock", 1:nstocks)
ri_true <- c(0.1, 0.3, 0.5)
ki_true <- c(500, 300, 100)
res_true <- c("Low", "Medium", "Medium")
e_true <- seq(0.02,0.25,length.out=nyrs)
e_scalar_true <- c(1, 0.8, 0.6)

# Simulate data
# B[t+1] = bt[t] + r*bt[t]*(1-bt[t]/k) - ct[t]
B_mat_true <- matrix(data=NA, ncol=nstocks, nrow=nyrs)
C_mat <- matrix(data=NA, ncol=nstocks, nrow=nyrs)
B_mat_true[1,] <- ki_true
C_mat[1,] <- ki_true * e_true[1] * e_scalar_true
for(i in 1:nstocks){
  for(j in 2:nyrs){
    B_mat_true[j,i] <- B_mat_true[j-1,i] + ri_true[i] * B_mat_true[j-1,i] * (1-B_mat_true[j-1,i]/ki_true[i]) - C_mat[j-1,i]
    C_mat[j,i] <-  B_mat_true[j,i] * e_true[j] * e_scalar_true[i]
  }
}

# Calculate fishing mortality time series
U_mat_true <- C_mat / B_mat_true

# Plot simulated data
plot_catch(C_mat, yrs)
plot_biomass(B_mat_true, yrs)
plot_f(U_mat_true, yrs)

# Management quantities
msy_true <- ri_true * ki_true / 4
bmsy_true <- ki_true / 2
fmsy_true <- ri_true / 2
bbmsy_mat_true <- t(t(B_mat_true) / bmsy_true)
bbmsy_end_true <- bbmsy_mat_true[nyrs,]
status_end_true <- bbmsy2status(bbmsy_end_true)

# Model functions
################################################################################

# Parameters
params <- c("r", "k")
nparam <- length(params)

# Starting values
par_init <- c(0.15, 0.29, 0.55, 470, 320, 110)

# Multi-species catch-only model
# For testing: par <- par_init
obj_func <- function(par){
  
  # Parameters
  ri <- par[1:3]
  ki <- par[4:6] 

  # Build biomass time series
  # B[t+1] = bt[t] + r*bt[t]*(1-bt[t]/k) - ct[t]
  B_mat <- matrix(data=NA, ncol=nstocks, nrow=nyrs)
  B_mat[1,] <- ki
  for(i in 1:nstocks){
    for(j in 2:nyrs){
      B_mat[j,i] <- B_mat[j-1,i] + ri[i] * B_mat[j-1,i] * (1-B_mat[j-1,i]/ki[i]) - C_mat[j-1,i]
    }
  }
  
  # Are all biomasses positive?
  # If biomass remains positive, calculate NLL
  # If biomass goes negative, report a large NLL
  pos_check <- sum(B_mat < 0)==0
  if(pos_check==T){
    
    # Derive fishing mortality time series
    F_mat <- C_mat / B_mat
    
    # Calculate log-ratio of secondary stocks to reference stock (first stock)
    F_ref <- F_mat[,2]
    F_other <- F_mat[,c(1,3)]
    z <- log(F_other / F_ref)
    z_avg <- apply(z, 2, mean)
    
    # Calculate SSQ of log-ratio and nll
    sq <- (z-z_avg)^2
    ssq <- colSums(sq)
    nll_E <- sum(-(nyrs/2)*log(ssq))
    
  }else{
    
    nll_E <- 9999
    
  }
  
  # Return NLL
  nll <- nll_E
  return(nll)
  
}

# Calculate r and k priors
r_priors <- calc_r_priors(res_true)
k_priors <- calc_k_priors(C_mat, r_priors)

# Create lower/upper bounds to pass to optim()
lower_bounds <- c(r_priors[,"r_lo"],k_priors[,"k_lo"])
upper_bounds <- c(r_priors[,"r_hi"],k_priors[,"k_hi"])

# Estimate parameters using optim()
# optfit <- optim(par=par_init, fn=obj_func)
optfit <- optim(par=par_init, fn=obj_func,
                lower=lower_bounds, upper=upper_bounds, method="L-BFGS-B")
optfit$convergence




