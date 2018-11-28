

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

# Directories
datadir <- "data/simulated/data"
codedir <- "code/helper_functions"
outdir <- "data/output"
plotdir <- "figures"

# Read species
spp <- read.csv(file.path(datadir, "demersal_trawl_fishery_10spp.csv"), as.is=T)

# Read r priors by family
r_vals_fam <- read.csv("data/priors/data/r_priors_by_family.csv", as.is=T)

# Source FishLife r prior code
source("data/priors/fishlife2/fishlife_r_priors.R")

# Helper functions
################################################################################

# R priors
calc_r_priors <- function(key){
  
  # Get FishLife r priors
  fl_r_list <- lapply(key$species, function(x) r_prior(x))
  fl_r_df <- do.call("rbind", fl_r_list)
  fl_r_df_format <- fl_r_df %>% 
    select(species, ln_r_mu, ln_r_sd) %>% 
    mutate(species=as.character(species))
  
  # Build r prior table
  r_priors <- key %>% 
    select(stock, species, family, resilience) %>% 
    # Add FishLife r priors
    left_join(fl_r_df_format, by="species") %>% 
    # Add Neubauer family-level r priors
    left_join(r_vals_fam, by="family") %>% 
    # Set resilience priors and record final prior source
    mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")),
           r_lo=c(0.01, 0.01, 0.01, 0.50)[resilience],
           r_hi=c(0.15, 0.50, 1.00, 1.25)[resilience],
           use=ifelse(!is.na(ln_r_mu), "FishLife", ifelse(!is.na(r_mu), "family", "resilience")))
  
  # Print message if resilience priors are used
  if(any(r_priors$use=="resilience")){
    stocks_with_res_r <- r_priors$stock[r_priors$use=="resilience"]
    print(paste("r priors for the following stocks are based on resilience instead of family: ", paste( stocks_with_res_r, collapse=", ")))
  }
  return(r_priors)
}

# Plot data
################################################################################

# Calculate priors for species of interest
spp1 <- spp %>% 
  rename(stock=comm_name) %>% 
  arrange(stock)
data <- calc_r_priors(spp1)

# Setup figure
figname <- "FigS2_fishlife_r_priors.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3.5, units="in", res=600)
par(mfrow=c(2,5), mar=c(2.5, 0.5, 1, 0.5), mgp=c(2,0.7,0), oma=c(1,0,0,0), xpd=NA)

# Loop through species
for(i in 1:nrow(data)){
  
  # Species
  spp <- data$stock[i]
  
  # Plot Neubauer et al. (2013) distribution
  if(!is.na(data$r_mu[i])){
    x <- seq(0,1,0.005)
    hx2 <- dlnorm(x, meanlog=log(data$r_mu[i]), sdlog=data$r_sd[i])
    ymax <- max(hx2)
    ymin <- 0 - ymax*0.05
    plot(x, hx2, type="n", bty="n", cex.axis=0.8, lwd=0.8, ylim=c(ymin, ymax),
         xlim=c(0,1), yaxt="n", xlab="", ylab="", main=spp, cex.main=0.8)
    polygon(x=c(x, rev(x)), y=c(rep(0,length(x)), rev(hx2)), col="grey80", border=F)
  }else{
    plot(1, 1, type="n", bty="n", cex.axis=0.8, lwd=0.8,
         xlim=c(0,1), yaxt="n", xlab="", ylab="", main=spp, cex.main=0.8)
  }
  
  # Add FishLife distribution
  x <- seq(0,1,0.005)
  hx <- dlnorm(x, meanlog=data$ln_r_mu[i], sdlog=data$ln_r_sd[i])
  ymax <- max(hx)
  ymin <- 0 - ymax*0.05
  par(new=T)
  plot(x, hx, type="n", bty="n", cex.axis=0.8, lwd=0.8, ylim=c(ymin, ymax),
       xlim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", main="")
  polygon(x=c(x, rev(x)), y=c(rep(0,length(x)), rev(hx)), col="grey50", border=F)
  
  # Add resilience prior
  lines(x=c(data$r_lo[i], data$r_hi[i]), y=c(ymin/2,ymin/2))
  
  # Add legend
  if(i==1){
    legend("topright", bty="n", legend=c("FishLife", "Neubauer", "Resilience"), cex=0.7,
           col=c("grey50", "grey80", "black"), pch=c(15, 15, NA), lty=c(NA, NA, 1), pt.cex=1.5)
  }
  
}

# Add axis label
mtext("Intrinsic growth rate, r", outer=T, side=1, adj=0.5, line=-0.5, cex=0.7, font=2)

# Off
dev.off()
