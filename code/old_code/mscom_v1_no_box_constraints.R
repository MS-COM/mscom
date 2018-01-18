
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


# Plotting functions
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

# Plot B/BMSY time series
plot_bbmsy <- function(mat, yrs){
  nstocks <- ncol(mat)
  mmax <- max(mat)
  colors <- brewer.pal(nstocks, "Set1")
  plot(x=yrs, y=mat[,1], bty="n", type="l", las=1, col=colors[1],
       ylim=c(0,mmax), xlab="", ylab="B/BMSY", main="Stock status") 
  for(i in 2:nstocks){
    lines(x=yrs, y=mat[,i], col=colors[i])
  }
  abline(h=0.5, lty=2)
  legend("bottomleft", bty="n", legend=paste("Stock", 1:nstocks), col=colors, lty=1)
}

# Categorical status
bbmsy2status <- function(bbmsy){
  as.character(cut(bbmsy, breaks=c(0,0.5,1.5,999), labels=c("over", "fully", "under")))
}

# Simulate data
################################################################################

# Params
yrs <- 1970:2010
nyrs <- length(yrs)
nstocks <- 3
stocks <- paste("Stock", 1:nstocks)
r_true <- c(0.1, 0.3, 0.5)
k_true <- c(500, 300, 100)
msy_true <- r_true*k_true/4
bmsy_true <- k_true/2
fmsy_true <- r_true/2
res_true <- c("Low", "Medium", "Medium")
e_true <- seq(0.02,0.25,length.out=nyrs)
e_scalar_true <- c(1, 0.8, 0.6)

# Simulate data
# B[t+1] = bt[t] + r*bt[t]*(1-bt[t]/k) - ct[t]
B_mat_true <- matrix(data=NA, ncol=nstocks, nrow=nyrs)
C_mat <- matrix(data=NA, ncol=nstocks, nrow=nyrs)
B_mat_true[1,] <- k_true
C_mat[1,] <- k_true * e_true[1] * e_scalar_true
for(i in 1:nstocks){
  for(j in 2:nyrs){
    B_mat_true[j,i] <- B_mat_true[j-1,i] + r_true[i] * B_mat_true[j-1,i] * (1-B_mat_true[j-1,i]/k_true[i]) - C_mat[j-1,i]
    C_mat[j,i] <-  B_mat_true[j,i] * e_true[j] * e_scalar_true[i]
  }
}

# Calculate other time series
U_mat_true <- C_mat / B_mat_true
BBMSY_mat_true <- t(t(B_mat_true) / bmsy_true)

# Plot simulated data
par(mfrow=c(2,2))
plot_catch(C_mat, yrs)
plot_biomass(B_mat_true, yrs)
plot_f(U_mat_true, yrs)
plot_bbmsy(BBMSY_mat_true, yrs)



# Model functions
################################################################################

# Parameters
params <- c("r", "k")
nparam <- length(params)

# Starting values
par_init <- matrix(data=c(0.15,550,
                          0.35,250,
                          0.55,150), byrow=F, nrow=nparam, ncol=nstocks, dimnames=list(params, stocks))

# Multi-species catch-only model
# For testing: par <- par_init
mscom_nll <- function(par){
  
  # Parameters
  par <- matrix(par, nrow=nparam, ncol=nstocks, byrow=F, dimnames=list(params, stocks))
  ri <- par[1,]
  ki <- par[2,] 

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

# Estimate parameters using optim()
optfit <- optim(par=par_init, fn=mscom_nll, control=list(maxit=10000))


# Extract parameters and build time series
################################################################################

# Extract results
get_results <- function(msfit){
  
  # Parameters
  par_est <- exp(msfit$par)
  r_est <- par_est[1:3]
  k_est <- par_est[4:6]
  msy_est <- r_est*k_est/4
  bmsy_est <- k_est/2
  fmsy_est <- r_est/2
  params <- data.frame(stocks, r_est, k_est, msy_est, bmsy_est, fmsy_est)
  
  # Build data
  B_mat <- matrix(data=NA, ncol=nstocks, nrow=nyrs)
  B_mat[1,] <- k_est
  for(i in 1:nstocks){
    # Extract coefs
    r <- r_est[i]
    k <- k_est[i]
    # Build trajectories
    for(j in 2:nyrs){
      B_mat[j,i] <- B_mat[j-1,i] + r * B_mat[j-1,i] * (1-B_mat[j-1,i]/k) - C_mat[j-1,i]
    }
  }
  
  # Other time series
  U_mat <- C_mat / B_mat
  BBMSY_mat <- t(t(B_mat) / bmsy_est)
  
  # Return output
  output <- list(params=params, B_mat=B_mat, U_mat=U_mat, BBMSY_mat=BBMSY_mat)
  return(output)
  
}


# Extract parameters
par <- optfit$par
r_pred <- par["r",]
k_pred <- par["k",]
msy_pred <- r_pred * k_pred / 4
bmsy_pred <- k_pred / 2
fmsy_pred <- r_pred / 2

# Build data
B_mat <- matrix(data=NA, ncol=nstocks, nrow=nyrs)
B_mat[1,] <- par["k",]
for(i in 1:nstocks){
  # Extract coefs
  r <- par["r",i]
  k <- par["k",i]
  # Build trajectories
  for(j in 2:nyrs){
    B_mat[j,i] <- B_mat[j-1,i] + r * B_mat[j-1,i] * (1-B_mat[j-1,i]/k) - C_mat[j-1,i]
  }
}
F_mat <- C_mat / B_mat
bbmsy_mat <- t(t(B_mat) / bmsy_pred)
bbmsy_end <- bbmsy_mat[nyrs,]
status_end <- bbmsy2status(bbmsy_end)

# Compare predictions
pmat <- data.frame(stocks, 
                   r_true=ri_true, r_pred,
                   k_true=ki_true, k_pred,
                   msy_true, msy_pred,
                   bmsy_true, bmsy_pred,
                   fmsy_true, fmsy_pred,
                   bbmsy_end_true, bbmsy_end,
                   status_end_true, status_end)

# Plot predictions
plot_biomass(B_mat, yrs)
plot_f(F_mat, yrs)













  
