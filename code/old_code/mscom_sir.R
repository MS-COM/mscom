
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
# For testing: catch<-C_mat; years<-yrs; stocks<-stocks; res<-res_true
mscom <- function(catch, years, stocks, res){
  
  # Time series info
  nyrs <- length(years)
  nstocks <- length(stocks)
  
  # Calculate r and k priors
  # R prior based on resilience; K prior based on max catch and r prior
  r_priors <- calc_r_priors(res_true)
  k_priors <- calc_k_priors(catch, r_priors)
  r_priors_ln <- log(r_priors)
  k_priors_ln <- log(k_priors)
  
  # Randomly sample r-k pairs in log-space
  npairs <- 10000
  ri <- sapply(1:nstocks, function(x) exp(runif(npairs, r_priors_ln[x,1], r_priors_ln[x,2])))
  ki <- sapply(1:nstocks, function(x) exp(runif(npairs, k_priors_ln[x,1], k_priors_ln[x,2])))
  
  # Loop through stocks
  for(i in 1:nstocks){
    
    # Get catch
    c_vec <- catch[,i]
    
    # Get r/k pairs to evaluate
    ri1 <- ri[,i]
    ki1 <- ki[,i]
    rk_pairs <- cbind(r=ri1, k=ki1, viable=rep(NA, npairs))
    
    # Plot r/k pairs to evaluate
    plot(x=ri1, y=ki1, log="xy", bty="n", las=1,
         xlim=r_priors[i,], ylim=k_priors[i,], xlab="r", ylab="k", col="gray95")
    
    # Loop through r/k pairs to see if viable
    for(j in 1:npairs){
      r <- ri1[j]
      k <- ki1[j]
      b_vec <- vector()
      b_vec[1] <- k
      for(yr in 2:nyrs){
        b_vec[yr] <- b_vec[yr-1] +  r*b_vec[yr-1] * (1-b_vec[yr-1]/k) - c_vec[yr-1]
        # if(b_vec[yr]<0){
        #   rk_pairs[j,"viable"] <- 0
        #   break
        # }
      }
      plot(b_vec)
    }
  
    
  }
  
  
}
  
plot(x=ri[,1], y=ki[,2], log="xy", xlim=r_prior[1,], ylim=k_prior[1,])



# Identify viable r/k pairs with cMSY
out1 <- cmsy2(year=years, catch=C_mat[,1], resilience=res_true[1])
out2 <- cmsy2(year=years, catch=C_mat[,2], resilience=res_true[2])
out3 <- cmsy2(year=years, catch=C_mat[,3], resilience=res_true[3])






