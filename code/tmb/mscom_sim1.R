
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
# e_scalar_true <- c(0.6, 0.6, 0.6)
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


input_data <- NULL
input_data$years <- yrs
input_data$catch <- C_mat

true_vals <- NULL
true_vals$r <- r_true
true_vals$K <- k_true
true_vals$B_ts <- B_mat_true
true_vals$U_ts <- U_mat_true
true_vals$Bmsy <- bmsy_true
true_vals$Umsy <- fmsy_true

wd <- "C:\\merrill\\MS-COM\\mscom\\code\\tmb"
saveRDS(input_data, file.path(wd, "simtest.rds"))
saveRDS(true_vals, file.path(wd, "simtest_true.rds"))

