
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

# Load scenarios
load(file.path(datadir, "simulated_multispecies_fisheries.Rdata"))


# Build data
################################################################################

# Calculate stats
stats <- sims %>% 
  # Derive values used to set priors
  group_by(iter, stock) %>% 
  summarize(cmax=max(catch),
            cfinal=catch[year==max(year)],
            cratio=cfinal/cmax,
            yr1=min(year),
            b1=biomass[year==min(year)],
            b2=biomass[year==max(year)]) %>% 
  # Add carrying capacity
  left_join(select(fdata, comm_name, k), by=c("stock"="comm_name")) %>% 
  # Calculate stats used to set priors
  mutate(s1=b1/k,
         s2=b2/k,
         kdivcmax=k/cmax) %>% 
  # Add scenario info
  left_join(select(scenarios, iter_id, ed, id, nyrs), by=c("iter"="iter_id"))

# For visualizing cMSY prior
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
figname <- "Fig6_priors_on_simstocks.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3, units="in", res=600)
layout(matrix(data=1:3, ncol=3), widths=c(0.2, 0.4, 0.4))
par(mar=c(3.5, 3.7, 2, 0.5), mgp=c(2.5,0.8,0), xpd=NA)

# Carrying capacity
##################################

# Carrying capacity
boxplot(stats$kdivcmax, log="y", las=1, frame=F, lty=1, col="grey95", lwd=0.9,
        cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex=1.2, pch=16,
        ylim=c(1,100), ylab=expression("K / C"["max"]), main="Carrying capacity, K")
# mtext("C", side=3, adj=0.05, line=-1.5, font=2)

# Add priors
lines(x=c(0.7, 0.7), y=c(1,100), col=tcolor("orange", 0.5), lwd=3)
lines(x=c(0.75, 0.75), y=c(2,25), col=tcolor("darkgreen", 0.3), lwd=3)

# Initial saturation
##################################

# Initial saturation
plot(s1 ~ yr1, stats, bty="n", type="n", las=2, ylim=c(0,2),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2, xlim=c(1880, 2000),
     xlab="", ylab="Saturation", main="Initial saturation")
# mtext("A", side=3, adj=0.05, line=-1.5, font=2)

# Add MS-cMSY prior
polygon(x=c(1880, 1945, 1980, 2000, 2000, 1980, 1945, 1880),
        y=c(0.75, 0.75, 0.2, 0.2, 1, 1, 1, 1),
        col=tcolor("darkgreen", 0.3), border=F)

# Add cMSY prior
polygon(x=c(1880,1960,1960,2000,2000, 1960,1960, 1880),
        y=c(0.5, 0.5, 0.2, 0.2, 0.6, 0.6, 0.9, 0.9),
        col=tcolor("orange", 0.5), border=F)

# Add points
points(x=stats$yr1, y=stats$s1, pch=16, col=c("grey20", "grey50")[factor(stats$nyrs)], cex=1.2)

# Add legend
legend("topright", bty="n", legend=c("25 yr", "50 yr"), pch=16, col=c("grey20", "grey50"))

# Full saturation line
lines(x=c(1880,2000), y=c(1,1), lwd=1.2, lty=3, col="black")
text("Fully saturated", x=1880, y=2000, pos=2, col="black", cex=0.9, xpd=NA)

# Final saturation
##################################

# Final saturation
plot(s2 ~ cratio, stats, bty="n", type="n", las=1, xlim=c(0,1), ylim=c(0,2),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2,
     xlab=expression("C"["final"]*" / C"["max"]), ylab="Saturation", main="Final saturation")
# mtext("B", side=3, adj=0.05, line=-1.5, font=2)

# Add MS-cMSY prior
x <- seq(0,1,0.1)
y1 <- 0 + 0.2*x
# y2 <- 0.8 + 0.2*x
y2 <- 0.5 + 0.4*x
polygon(x=c(x, rev(x)), y=c(y2, rev(y1)), col=tcolor("darkgreen", 0.3), border=F)

# Add cMSY prior
x <- cmsy_priors$cr
y2 <- cmsy_priors$s_min
y1 <- cmsy_priors$s_max
polygon(x=c(x, rev(x)), y=c(y1, rev(y2)), col=tcolor("orange", 0.5), border=F)

# Add points
points(x=stats$cratio, y=stats$s2, pch=16, col=c("grey20", "grey50")[factor(stats$ed)], cex=1.2)

# Add legend
legend("topright", bty="n", legend=c("one-way", "two-way"), pch=16, col=c("grey20", "grey50"))

# Full saturation line
lines(x=c(0,1), y=c(1,1), lwd=1.2, lty=3, col="black")
text("Fully saturated", x=1.0, y=1.08, pos=2, col="black", cex=1.1, xpd=NA)

###
dev.off()
