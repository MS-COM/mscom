
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Turn off sci notation
options(scipen=999)

# Packages
library(rio)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "data/test_cases/data"
plotdir <- "figures"

# Read data
catch <- read.csv(file.path(datadir, "gulf_of_thailand_trawl_catch_data.csv"), as.is=T)
load(file.path(datadir, "got_ms_cmsy_predictions.Rdata"))

# Format catch data
catch_w <- na.omit(dcast(catch, year ~ comm_name, value.var="catch"))


# Visualize data
################################################################################

# Species
species <- out$stocks

# Setup figure
figname <- "Fig6_GOT_predictions.png"
png(file.path(plotdir, figname), width=6.5, height=6.5, units="in", res=600)
par(mfcol=c(4, length(species)), mar=c(3,3,0.5,0.5), oma=c(0,2,1,0), mgp=c(3,0.8,0), xpd=NA)

# Loop through species
for(i in 1:length(species)){
  
  # Species
  spp <- species[i]
  yrs <- out$yrs
  yr1 <- freeR::floor1(min(yrs),5)
  yr2 <- freeR::ceiling1(max(yrs),5)
  
  # Extract stuff
  ids <- out$id_try
  rs <- out$r_try[,i]
  ks <- out$k_try[,i]
  id_rk_viable <- out$id_rk_v[[i]]
  b_viable <- out$b_v[[i]]
  bbmsy_viable <- out$bbmsy_v[[i]]
  er_viable <- out$er_v[[i]]
  s1_priors <- out$s1_priors
  s2_priors <-out$s2_priors
  r_priors <- out$r_priors
  k_priors <- out$k_priors
  er_vv <- out$er_vv[[i]]
  bbmsy_vv <- out$bbmsy_vv[[i]]
  bbmsy_vv_median <- out$bbmsy_vv_median[[i]]
  top_corr <- out$top_corr
  
  # Plot catch
  #############################################
  
  # Plot catch
  ylabel <- ifelse(i==1, "Catch (1000s mt)", "")
  c <- catch_w[,1+i]/1000
  cmax <- freeR::ceiling1(max(c),5)
  plot(x=yrs, y=c, bty="n", type="l", las=2, cex.lab=1.2, main=spp,
       xlim=c(yr1, yr2), ylim=c(0, cmax), xlab="", ylab=ylabel)
  
  
  # Plot r/k pairs
  #############################################
  
  # Plot r/k pairs
  ylabel <- ifelse(i==1, "Carring capacity, K", "")
  plot(ks/1000 ~ rs, log="xy", type="n", bty="n", las=1, pch=15, col="gray80", cex.lab=1.2,
       xlim=r_priors[i,], ylim=k_priors[i,]/1000,
       # xlim=c(0.01,1.2), ylim=c(0.1, k_priors[i,2]/1000),
       xlab="", ylab=ylabel)
  if(i==2){mtext("Intrinsic growth rate, r", side=1, adj=0.5, cex=0.75, line=2)}
  
  # Add viable r/k pairs
  # Potentially reduce this to unique r/k pairs
  # There could be redundancy when evaluating multiple IDs
  points(id_rk_viable$r, id_rk_viable$k/1000, pch=15, col="grey70")
  
  # Add most highly correlated pairs
  id_rk_v_ind <- unlist(top_corr[,paste0("index", i)])
  rk_corr <- subset(id_rk_viable, index %in% id_rk_v_ind)
  points(x=rk_corr$r, y=rk_corr$k/1000, pch=15, col=freeR::tcolor("darkorange", 0.6))

  # Plot status estimate
  #############################################
  
  # Plot BBMSY trajectories
  ylabel <- ifelse(i==1, expression("B/B"["MSY"]), "")
  ymax <- freeR::ceiling1(max(bbmsy_viable, na.rm=T), 0.5)
  plot(bbmsy_viable[,1] ~ yrs, type="n", bty="n", las=2, cex.lab=1.2,
       xlim=c(yr1, yr2), ylim=c(0, ymax), xlab="", ylab=ylabel)
  for(k in 1:ncol(bbmsy_viable)){lines(x=yrs, y=bbmsy_viable[,k], col="grey70")}
  for(k in 1:ncol(bbmsy_vv)){lines(x=yrs, y=bbmsy_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
  lines(x=yrs, y=bbmsy_vv_median, lwd=1.5, col="black")
  lines(x=c(yr1, yr2), y=c(0.5, 0.5), lty=3)
  lines(x=yrs, y=apply(bbmsy_viable,1,median), lwd=1.5, lty=3, col="green")
  
  # Plot exploitation rate
  #############################################
  
  # Plot exploitation trajectories
  ylabel <- ifelse(i==1, "Exploitation rate", "")
  ermax <- freeR::ceiling1(max(er_viable, na.rm=T), 0.5)
  plot(er_viable[,1] ~ yrs, type="n", bty="n", las=2, cex.lab=1.2,
       xlim=c(yr1, yr2), ylim=c(0, ymax), xlab="", ylab=ylabel)
  for(k in 1:ncol(er_viable)){lines(x=yrs, y=er_viable[,k], col="grey80")}
  for(k in 1:ncol(er_vv)){lines(x=yrs, y=er_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
  
  # Add observed ER
  er <- subset(catch, comm_name==spp & year%in%1995:2015)
  lines(x=er$year, y=er$er, lwd=1.3, col="red")
  
  
}

# Off
dev.off()



