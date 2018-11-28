
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
outdir <- "data/output"
simdir <- "data/simulated/data"
plotdir <- "figures"

# Read output
results <- read.csv(file.path(outdir, "simtest_results.csv"))

# Read sim parameters
params <- read.csv(file.path(simdir, "demersal_trawl_fishery_10spp.csv"), as.is=T)

# Look for example iteration
################################################################################

# Look for example
ex <- results %>% 
  filter(nstocks==5 & year==2015) %>% 
  filter(abs(mscmsy_diff)<10)

# Do a lot of these good ones occur in a single scenario?
sort(table(ex$scenario)) # 4 in here: unfished_one-way_0.15_0_5_25

# Load scenario
iter_id <- as.character(unique(ex$iter[ex$scenario=="unfished_one-way_0.15_0_5_25"]))
load(file.path(outdir, paste0(iter_id, ".Rdata")))

# Plot data
################################################################################

# Extract info
output <- out
key <- output$key
preds <- output$preds
stocks <- key$stock
nstocks <- length(stocks)

# Setup figure
figname <- "Fig3_mscmsy_viz.png"
png(file.path(plotdir, figname), width=6.5, height=5.5, units="in", res=600)
par(mfcol=c(4,nstocks), mar=c(3,2,2,1), oma=c(0,2,0,0), mgp=c(2.5,0.7,0), xpd=NA)

# Loop through species
for(i in 1:nstocks){
  
  # Stock
  stock_i <- stocks[i]
  
  # Year info
  #########################################
  
  # Year info
  sdata <- filter(preds, stock==stock_i)
  yrs <- sdata$year
  yr1 <- floor(min(yrs) / 10) * 10
  yr2 <- ceiling(max(yrs) / 10) * 10
  
  # Plot catch
  #########################################
  
  # Plot catch
  ymax <- max(sdata$catch)
  ylabel <- ifelse(i==1, "Catch", "")
  plot(catch ~ year, sdata, type="l", bty="n", las=2, xaxt="n",
       xlim=c(yr1, yr2), ylim=c(0, ymax), xlab="", ylab=ylabel, main=stock_i)
  axis(side=1, at=seq(yr1, yr2, 10), las=2)
  
  # Plot r/k pairs
  #########################################
  
  # Extract r/k info
  id_rk_v <- output$id_rk_v[[i]]
  id_rk_v$k <-   id_rk_v$k
  top_corr <- output$top_corr
  r_prior <- unlist(select(filter(output$r_priors, stock==stock_i), r_lo, r_hi))
  k_prior <- unlist(select(filter(output$k_priors, stock==stock_i), k_lo, k_hi))
  
  # Plot viable r/k pairs
  # Potentially reduce this to unique r/k pairs
  # There could be redundancy when evaluating multiple IDs
  rmin <- floor(min(id_rk_v$r) / 0.1) * 0.1
  rmax <- ceiling(max(id_rk_v$r) / 0.1) * 0.1
  kmin <- min(id_rk_v$k) 
  kmax <- max(id_rk_v$k)
  ylabel <- ifelse(i==1, "Carrying capacity, K", "")
  plot(k ~ r, id_rk_v, bty="n", las=1, pch=15, col="gray70",
       xlim=c(rmin, rmax), ylim=c(kmin, kmax), 
       xlab="Intrinsic growth rate, r", ylab=ylabel)
  
  # Add most highly correlated pairs
  id_rk_v_ind <- unlist(top_corr[,paste0("index", i)])
  rk_corr <- subset(id_rk_v, index %in% id_rk_v_ind)
  points(x=rk_corr$r, y=rk_corr$k, pch=15, col=freeR::tcolor("darkorange", 0.6))
  
  # # Add most common highly correlated pair
  # rk_mode <- mode(rk_v_ind)
  # rk_corr <- subset(rk_viable, index==rk_mode)
  # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="green")
  
  # # Add repeatedly highly correlated r/k pairs
  # rk_corr <- subset(rk_viable, ncorr>=5)
  # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="darkorange")
  
  # Add truth
  rk <- filter(params, comm_name==stock_i)
  points(x=rk$r, y=rk$k, col="red", pch=15, cex=1.2)
  
  # 
  # # Add legend
  # if(i==1){
  #   legend("topright", bty="n", pch=15, pt.cex=1.3, cex=0.9,
  #          col=c("grey70", "darkorange"),
  #          legend=c("Viable", "Effort constrained"))
  # }
  
  # Plot exploitation trajectories
  #########################################
  
  # Extract ER trajectories
  er_v <- output$er_v[[i]]
  er_vv <- output$er_vv[[i]]
  
  # Plot exploitation trajectories
  # if(!missing(true)){
  #   ymax <- ceiling(max(er_v, true$ts$er, na.rm=T)/0.5) * 0.5
  # }else{
  #   ymax <- ceiling(max(er_v, na.rm=T)/0.5) * 0.5
  # }
  ylabel <- ifelse(i==1, "Exploitation rate", "")
  plot(er_v[,1] ~ yrs, type="n", bty="n", las=2, xaxt="n",
       xlim=c(yr1, yr2), ylim=c(0, 1.0), xlab="", ylab=ylabel)
  axis(side=1, at=seq(yr1, yr2, 10), las=2)
  for(k in 1:ncol(er_v)){lines(x=yrs, y=er_v[,k], col="grey80")}
  for(k in 1:ncol(er_vv)){lines(x=yrs, y=er_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
  lines(x=sdata$year, y=sdata$er, lwd=1.5, col="black")
  
  # Add truth
  ts <- filter(results, stock==stock_i & iter==iter_id)
  lines(x=ts$year, y=ts$er, col="red", lwd=1, lty=2)

  
  # Plot B/BMSY trajectories
  #########################################
  
  # Extract B/BMSY trajectories
  bbmsy_v <- output$bbmsy_v[[i]]
  bbmsy_vv <- output$bbmsy_vv[[i]]
  
  # Plot B/BMSY trajectories
  # if(!missing(true)){
  #   ymax <- ceiling(max(bbmsy_v, true$ts$bbmsy, na.rm=T)/0.5) * 0.5
  # }else{
  #   ymax <- ceiling(max(bbmsy_v, na.rm=T)/0.5) * 0.5
  # }
  ylabel <- ifelse(i==1, expression("B/B"["MSY"]), "")
  plot(bbmsy_v[,1] ~ yrs, type="n", bty="n", las=2, xaxt="n",
       xlim=c(yr1, yr2), ylim=c(0, 2.5), xlab="", ylab=ylabel)
  axis(side=1, at=seq(yr1, yr2, 10), las=2)
  # Viable trajectories
  for(k in 1:ncol(bbmsy_v)){lines(x=yrs, y=bbmsy_v[,k], col="grey70")}
  # Effort constrained trajectories
  for(k in 1:ncol(bbmsy_vv)){lines(x=yrs, y=bbmsy_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
  lines(x=sdata$year, y=sdata$bbmsy, lwd=1.5, col="black")
  # cMSY trajectory
  # lines(x=sdata$year, y=sdata$bbmsy_cmsy, lwd=1.2, col="black")
  # Overfishing line
  lines(x=c(yr1, yr2), y=c(0.5, 0.5), lty=3)
  lines(x=c(yr1, yr2), y=c(1, 1), lty=2)
  
  # Add truth
  lines(x=ts$year, y=ts$bbmsy, col="red", lwd=1, lty=2)
  
  
}

# Off
dev.off()

