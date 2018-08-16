
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(quantreg)

# Directories
plotdir <- "figures"

# Read RAM Legacy Database
load("/Users/cfree/Dropbox/Prelim Database Files/Versions/RAM v4.40 (6-4-18)/DB Files With Assessment Data/DBdata.RData")

# Load helper functions
source("code/ram_stocks/helper_functions.R")

# Build data
################################################################################

# Build stock data
colnames(timeseries_values_views )
data <- timeseries_values_views %>% 
  # Select columns of interest
  select(stockid, year, TCbest, BdivBmsypref) %>% 
  # Rename columns
  setNames(tolower(names(.))) %>% 
  rename(tc=tcbest, bbmsy=bdivbmsypref)

# Build useful reference points
refs <- bioparams_values_views %>% 
  # Get BMSY
  select(stockid, TBmsybest) %>% 
  rename(bmsy=TBmsybest) %>% 
  # Add BMSY units
  left_join(select(bioparams_units_views, stockid, TBmsybest), by="stockid") %>% 
  rename(bmsy_units=TBmsybest) %>% 
  # Add catch units
  left_join(select(timeseries_units_views, stockid, TCbest), by="stockid") %>% 
  rename(tc_units=TCbest) %>% 
  # Filter
  filter(!is.na(bmsy))

# Build data for examining cMSY prior performance
results <- data %>% 
  filter(!is.na(tc)) %>% 
  group_by(stockid) %>% 
  summarize(nyr=n(),
            cmax=max(tc),
            clast=tc[n()],
            cratio=clast/cmax,
            bbmsy1=bbmsy[1],
            bbmsy2=bbmsy[n()],
            s1=datalimited2::bbmsy2s(bbmsy1),
            s2=datalimited2::bbmsy2s(bbmsy2)) %>% 
  # Add BMSY to calculate K
  left_join(refs, by="stockid") %>% 
  mutate(k=bmsy*2) %>% 
  # How much larger is K than Cmax? Only if in same units
  mutate(cmax2k=ifelse(bmsy_units==tc_units, k/cmax, NA)) %>% 
  # Add species name
  left_join(select(stock, stockid, scientificname, commonname), by="stockid") %>% 
  rename(species=scientificname, comm_name=commonname) %>% 
  # Rearrange
  select(stockid, species, comm_name, everything())

# Get and add resilience
spp <- sort(unique(results$species))
res <- resilience(spp)
results <- results %>% 
  left_join(res, by="species")

# Helper stuff
################################################################################

# Calculate cMSY saturation prior
calc_cmsy_s_prior <- function(cr){
  if(cr>=0.8){s_prior <- c(0.4,0.8)}
  if(cr<0.8&cr>=0.5){s_prior <- c(0.2,0.6)}
  if(cr<0.5&cr>=0.35){s_prior <- c(0.01,0.4)}
  if(cr<0.35&cr>=0.15){s_prior <- c(0.01,0.3)}
  if(cr<0.15&cr>=0.05){s_prior <- c(0.01,0.2)}
  if(cr<0.05){s_prior <- c(0.01,0.1)}
  return(s_prior)
}

# Data for plotting polygon
cmsy_priors <- data.frame(cr=seq(0,1,0.001), s_min=NA, s_max=NA)
for(i in 1:nrow(cmsy_priors)){
  cmsy_priors[i,2:3]<- calc_cmsy_s_prior(cmsy_priors$cr[i])
}

# Plot data
################################################################################

# Setup figure
figname <- "Fig0_cmsy_prior_check.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6.5, units="in", res=600)
par(mfrow=c(2,2), mar=c(3.5, 3.5, 2, 0.5), mgp=c(2.2,0.8,0))

# A. Initial year saturation
##############################################

# Initial year saturation (cuts off outliers)
plot(s1 ~ nyr, results, type="n", bty="n", las=1, 
     xlim=c(0,150), ylim=c(0,5), main="Initial saturation prior",
     xlab="# of years of catch data", ylab="Initial year saturation")
polygon(x=c(0,50,50,150,150,50,50,0), y=c(0.2,0.2,0.5,0.5,0.9,0.9,0.6,0.6),
        col="coral", border=NA)
points(results$nyr, results$s1, col="grey30")
n <- nrow(filter(results, !is.na(s1)))
text(x=150, y=5, label=paste(n, "stocks"), pos=2, xpd=NA)

# B. Final year saturation
##############################################

# Final year saturation (cuts off outliers)
plot(s2 ~ cratio, results, type="n", bty="n", las=1, 
     xlim=c(0,1), ylim=c(0,3),  main="Final saturation prior",
     xlab=expression("C"["final"]*" / C"["max"]), ylab="Final year saturation")
polygon(x=c(cmsy_priors$cr, rev(cmsy_priors$cr)),
        y=c(cmsy_priors$s_min, rev(cmsy_priors$s_max)), col="coral", border=NA)
points(results$cratio, results$s2, col="grey30")
n <- nrow(filter(results, !is.na(s2)))
text(x=1, y=3, label=paste(n, "stocks"), pos=2, xpd=NA)

# C. Carrying capacity overall
##############################################

# Carrying capacity overall
boxplot(results$cmax2k, log="y", las=1, frame=F, lty=1,
        ylab=expression("K / C"["max"]),  main="K prior (Catch-MSY)")
n <- nrow(filter(results, !is.na(cmax2k)))
text(x=1.5, y=100, label=paste(n, "stocks"), pos=2, xpd=NA)

# Add lines marking resilience priors
lines(x=c(0.65,0.65), y=c(1,100), lwd=1.5)

# D. Carrying capacity by resilience
##############################################

# Carrying capacity by resilience
res_colors <- rev(brewer.pal(nlevels(results$resilience), "Spectral"))
boxplot(cmax2k ~ resilience, results, log="y", las=2, frame=F, lty=1, at=2:5, xlim=c(0,5.5),
        ylab=expression("K / C"["max"]),  main="K prior by resilience (cMSY)", col=res_colors)
n <- nrow(filter(results, !is.na(cmax2k) & !is.na(resilience)))
text(x=5.5, y=100, label=paste(n, "stocks"), pos=2, xpd=NA)

# Add prior lines
res_rs <- matrix(data=c(0.6,1.5,
                        0.2,0.8,
                        0.05, 0.5,
                        0.015,0.1), ncol=2, byrow=T)
for(i in 1:nrow(res_rs)){
  lines(x=c(0.5+0.2*i,0.5+0.2*i), y=c(1/res_rs[i,2], 4/res_rs[i,1]), lwd=1.5, col=res_colors[i])
  lines(x=c(0.6+0.2*i,0.6+0.2*i), y=c(2/res_rs[i,2], 12/res_rs[i,1]), lwd=1.5, col=res_colors[i])
}



# Off
dev.off()
graphics.off()



