
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
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "data/test_cases/data"
plotdir <- "figures"

# Read data
data <- import(file.path(datadir, "gulf_of_thailand_trawl_catch_data.csv"))


# Visualize data
################################################################################

# Setup figure
figname <- "FigX_GOT_predictions.png"
png(file.path(plotdir, figname), width=6.5, height=5, units="in", res=600)
species <- c("Blue swimming crab", "Bigeye scad", "Largehead hairtail")

# Loop through species
par(mfcol=c(3, length(species)), mar=c(3,3,0.5,0.5), oma=c(0,1,1,0), mgp=c(2.5,0.8,0), xpd=NA)
for(i in 1:length(species)){
  
  # Subset data
  spp <- species[i]
  sdata <- data %>% 
    filter(comm_name==spp) %>% 
    mutate(catch=catch/1000)
  
  # Plot catch
  xmin <- freeR::floor1(min(sdata$year), 10)
  xmax <- freeR::ceiling1(max(sdata$year), 10)
  ymax <- freeR::ceiling1(max(sdata$catch), 10)
  ylabel <- ifelse(i==1, "Catch (1000s mt)", "")
  plot(catch ~ year, sdata, bty="n", type="l", las=2, cex.lab=1.1, cex.axis=0.9,
       xlab="", ylab=ylabel, xlim=c(xmin, xmax), ylim=c(0, ymax))
  title(main=spp, line=0.5)
  
  # Plot exploitation rate
  ylabel <- ifelse(i==1, "Exploitation rate", "")
  plot(er ~ year, sdata, bty="n", type="l", las=2, cex.lab=1.1, cex.axis=0.9,
       xlab="", ylab=ylabel, ylim=c(0,1), xlim=c(xmin, xmax))
  
  # Plot status estimate
  ylabel <- ifelse(i==1, expression("B/B"["MSY"]), "")
  plot(er ~ year, sdata, bty="n", type="n", las=2, cex.lab=1.1, cex.axis=0.9,
       xlab="", ylab=ylabel, ylim=c(0,2.5), xlim=c(xmin, xmax))
  lines(x=c(xmin, xmax), y=rep(0.5, 2), lty=3)
  
  
}

# Off
dev.off()



