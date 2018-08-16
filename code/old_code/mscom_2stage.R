
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

# Model
################################################################################

# Load cMSY function
source("/Users/cfree/Dropbox/Chris/Rutgers/rpackages/datalimited2/R/cmsy2.R")

# Run cMSY on all three time series
out1 <- cmsy2(year=yrs, catch=C_mat[,1], resilience=res_true[1])
out2 <- cmsy2(year=yrs, catch=C_mat[,2], resilience=res_true[2])
out3 <- cmsy2(year=yrs, catch=C_mat[,3], resilience=res_true[3])

# Extract viable rk values
rv1 <- out1[["rv.all"]]
rv2 <- out2[["rv.all"]]
rv3 <- out3[["rv.all"]]
kv1 <- out1[["kv.all"]]*1000
kv2 <- out2[["kv.all"]]*1000
kv3 <- out3[["kv.all"]]*1000

# Extract viable biomass trajectories
btv1 <- out1$btv.all * 1000
btv2 <- out2$btv.all * 1000
btv3 <- out3$btv.all * 1000

# Calculate fishing mortality (U) trajectories
utv1 <- t(C_mat[,1] / t(btv1[,1:(ncol(btv1)-1)]))
utv2 <- t(C_mat[,2] / t(btv2[,1:(ncol(btv1)-1)]))
utv3 <- t(C_mat[,3] / t(btv3[,1:(ncol(btv1)-1)]))

# Calculate correlation matrix
utv <- rbind(utv1, utv2, utv3)
cov_mat <- utv %*% t(utv)
corr_mat <- cov2cor(cov_mat)


####

a<-length(rv1)
b<-length(rv2)
c<-length(rv3)

Max1<-NULL
for (i in a+1:a+b){
  for(j in 1:a){
    max_cor<-max(cov_mat[j,i])
    Max1<-rbind(max_cor)
  }
}

Max2<-NULL
for (i in a+b+1:a+b+c){
  for(j in 1:a){
    max_cor<-max(cov_mat[j,i])
    Max2<-rbind(max_cor)
  }
}

Max3<-NULL
for (i in 1:a){
  for(j in a+1:a+b){
    max_cor<-max(cov_mat[j,i])
    Max3<-rbind(max_cor)
  }
}

Max4<-NULL
for (i in a+b+1:a+b+c){
  for(j in a+1:a+b){
    max_cor<-max(cov_mat[j,i])
    Max4<-rbind(max_cor)
  }
}

Max5<-NULL
for (i in 1:a){
  for(j in a+b+1:a+b+c){
    max_cor<-max(cov_mat[j,i])
    Max5<-rbind(max_cor)
  }
}

Max6<-NULL
for (i in a+1:a+b){
  for(j in a+b+1:a+b+c){
    max_cor<-max(cov_mat[j,i])
    Max6<-rbind(max_cor)
  }
}



##################################

# Build index
ind <- c(rep(1, length(rv1)), 
         rep(2, length(rv2)),
         rep(3, length(rv3)))

# Identify break points
a <- length(rv1) + 1
b <- a + length(rv2)

# Stock indices
stock1ind <- 1:length(rv1)
stock2ind <- 1:length(rv2)
stock3ind <- 1:length(rv3)

# Build matrix to hold results
results <- matrix(data=NA, ncol=nstocks, nrow=length(ind))

# Loop through rows
for(i in 1:length(ind)){
  
  # Which stock?
  print(i)
  stock.now <- ind[i]
  
  # If stock==1, search over stocks 2+3
  if(stock.now==1){
    # Best stock in 2
    stock2 <- corr_mat[i,a:(b-1)]
    stock2max <- which.max(stock2)
    # Best stock in 3
    stock3 <- corr_mat[i,b:length(ind)]
    stock3max <- which.max(stock3)
    # Record results
    results[i,1] <- stock1ind[i] # stock 1
    results[i,2] <- stock2max # stock 2
    results[i,3] <- stock3max # stock 3
  }
  
  # If stock==2, search over stocks 1+3
  if(stock.now==2){
    # Best stock in 1
    stock1 <- corr_mat[i,1:(a-1)]
    stock1max <- which.max(stock1)
    # Best stock in 3
    stock3 <- corr_mat[i,b:length(ind)]
    stock3max <- which.max(stock3)
    # Record results
    results[i,1] <- stock1max # stock 1
    results[i,2] <- stock2ind[i] # stock 2
    results[i,3] <- stock3max # stock 3
  }
  
  # If stock==3, search over stocks 1+2
  if(stock.now==2){
    # Best stock in 2
    stock2 <- corr_mat[i,a:(b-1)]
    stock2max <- which.max(stock2)
    # Best stock in 3
    stock3 <- corr_mat[i,b:length(ind)]
    stock3max <- which.max(stock3)
    # Record results
    results[i,1] <- stock1max # stock 1
    results[i,2] <- stock2max # stock 2
    results[i,3] <- stock3ind[i] # stock 3
  }
  
}

sum(!is.na(results))
  
  
  









# Make stock 1 the reference stock
par(mfrow=c(1,3))
plot(rv1, kv1, log="xy", col="grey20", cex.axis=1.2, cex.lab=1.2)
points(r_true[1], k_true[1], pch=16, col="red")

plot(rv2, kv2, log="xy", col="grey20", cex.axis=1.2, cex.lab=1.2)
points(r_true[2], k_true[2], pch=16, col="red")

plot(rv3, kv3, log="xy", col="grey20", cex.axis=1.2, cex.lab=1.2)
points(r_true[3], k_true[3], pch=16, col="red")



