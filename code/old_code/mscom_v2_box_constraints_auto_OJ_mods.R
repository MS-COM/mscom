
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

# Multi-species catch-only model
# For testing: par <- start_vals_log
mscom_nll <- function(par, C_mat, nstocks, nyrs, r_priors, k_priors){
  
  # Parameters
  ri_log <- par[1:3]
  ki_log <- par[4:6]
  ri <- exp(ri_log)
  ki <- exp(ki_log)
  
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
  pos_check <- sum(B_mat < 0)==0 & sum(is.na(B_mat))==0
  if(pos_check==F){
    nll_E <- 9999
  }else{
    
    # Derive fishing mortality time series
    F_mat <- C_mat / B_mat
    
    # Calculate log-ratio of secondary stocks to reference stock (first stock)
    F_ref <- F_mat[,3]
    F_other <- F_mat[,1:2]
    z <- log(F_other / F_ref)
    z_avg <- apply(z, 2, mean)
    q_other = exp(z_avg)
    E_other=F_other/q_other
    
    # Calculate SSQ of log-ratio and nll
    #sq <- (z-z_avg)^2
    #ssq <- colSums(sq)
    #nll_E <- sum(-(nyrs/2)*log(ssq))
    
    #Calculate SSQ and likelihood component related to shared effort
    sq=(E_other-F_other)^2
    nll_E = -sum(dnorm(sq, mean=0, sd=0.5, log=T))
    
    # Add normal priors on r and K
    # nll_r1 <- -dnorm(ri[1], mean = 0.3, sd = 1, log=T)
    # nll_r2 <- -dnorm(ri[2], mean = 0.3, sd = 1, log=T)
    # nll_r3 <- -dnorm(ri[3], mean = 0.3, sd = 1, log=T)
    # nll_k1 <- -dnorm(ki[1], mean = 300, sd = 1000, log=T)
    # nll_k2 <- -dnorm(ki[2], mean = 300, sd = 1000, log=T)
    # nll_k3 <- -dnorm(ki[3], mean = 300, sd = 1000, log=T)
    
    # Add log-normal priors on r and K
    nll_r1 <- -dlnorm(ri[1], meanlog=0.3, sdlog=1, log=T)
    nll_r2 <- -dlnorm(ri[2], meanlog=0.3, sdlog=1, log=T)
    nll_r3 <- -dlnorm(ri[3], meanlog=0.3, sdlog=1, log=T)
    nll_k1 <- -dlnorm(ki[1], meanlog=300, sdlog=1000, log=T)
    nll_k2 <- -dlnorm(ki[2], meanlog=300, sdlog=1000, log=T)
    nll_k3 <- -dlnorm(ki[3], meanlog=300, sdlog=1000, log=T)
    
    # Cumulative NLL
    nll_E = nll_r1 + nll_r2 + nll_r3 + nll_k1 + nll_k2 + nll_k3 + nll_E
    
  }
  
  # Return NLL
  nll <- nll_E
  return(nll)
  
}

# Multi-species catch-only model
# For testing: C_mat<-C_mat; years<-yrs; stocks<-stocks; res<-res_true
mscom <- function(C_mat, years, stocks, res){
  
  # Time series info
  nyrs <- length(years)
  nstocks <- length(stocks)
  
  # Calculate r and k priors
  # R prior based on resilience; K prior based on max catch and r prior
  r_priors <- calc_r_priors(res_true)
  k_priors <- calc_k_priors(C_mat, r_priors)
  
  # Create lower/upper bounds to pass to optim()
  lower_bounds_log <- log(c(r_priors[,"r_lo"],k_priors[,"k_lo"]))
  upper_bounds_log <- log(c(r_priors[,"r_hi"],k_priors[,"k_hi"]))
  
  # Specify initial values
  # Median of r and k priors for each stock
  r_starts <- apply(r_priors, 1, median)
  k_starts <- apply(k_priors, 1, median)
  start_vals_log <- log(c(r_starts, k_starts))
  
  # Fit model using optim() --- original version
  #optfit <- optim(par=start_vals_log, fn=mscom_nll,
  #                C_mat=C_mat, nstocks=nstocks, nyrs=nyrs, 
  #                lower=lower_bounds_log, upper=upper_bounds_log, method="L-BFGS-B",
  #                control=list(trace=1, maxit=10000))
  
  # Fit model using optim() ---- play version
  optfit <- optim(par=start_vals_log, fn=mscom_nll,
                  C_mat=C_mat, nstocks=nstocks, nyrs=nyrs, 
                  method="Nelder-Mead",
                  control=list(trace=1, maxit=10000))
  
  # Return fit
  return(optfit)
  
}

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


# Test model
################################################################################

# Fit MSCOM
msfit <- mscom(C_mat=C_mat, years=yrs, stocks=stocks, res=res_true)

# Get MSCOM results
msout <- get_results(msfit)

# Plot results
par(mfrow=c(1,1))
plot_bbmsy(msout[["BBMSY_mat"]], yrs)


# Plot r/k priors
###################################

# Get estimated parameters
params <- msout$params

par(mfcol=c(2,3))
for(i in 1:nstocks){
  
  # Plot R prior
  r_obs <- r_true[i]
  r_est <- params$r_est[i]
  x <- seq(0,1,0.01)
  # y <- dnorm(x, mean=0.3, sd=1)
  y <- dlnorm(x, meanlog=0.3, sdlog=1)
  plot(1:10, 1:10, type="n", bty="n", las=1, xlab="r", ylab="", 
       xlim=c(0,1), ylim=c(0,max(y)))
  polygon(x=c(x, rev(x)), 
          y=c(rep(0, length(x)), rev(y)), col="grey80", border=F)
  abline(v=r_obs, lty=1, lwd=2)
  abline(v=r_est, lty=2)
  
  # Plot K prior
  k_obs <- k_true[i]
  k_est <- params$k_est[i]
  x <- seq(0,1000,1)
  y <- dnorm(x, mean=300, sd=1000)
  # y <- dlnorm(x, meanlog=300, sdlog=1000)
  plot(1:10, 1:10, type="n", bty="n", las=1, xlab="k", ylab="", 
       xlim=c(0,1000), ylim=c(0,max(y)))
  polygon(x=c(x, rev(x)), 
          y=c(rep(0, length(x)), rev(y)), col="grey80", border=F)
  abline(v=k_obs, lty=1, lwd=2)
  abline(v=k_est, lty=2)
  
}




x <- 0:1000
y <- dlnorm(x, meanlog=log(500, sdlog=1000)
plot(x, y, type="l")

x <- seq(0,1,0.01)
# y <- dnorm(x, mean=0.3, sd=1)
y <- stats::dlnorm(x, meanlog=0.3, sdlog=2)
plot(x, y)

