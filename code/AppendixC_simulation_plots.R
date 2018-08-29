
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
datadir <- "data/simulated"
plotdir <- "figures"


# Format data
################################################################################

# Load scenarios
load(file.path(datadir, "scenarios_pelagic_longline.Rdata"))

# Load simulatons
load(file.path(datadir, "sim1_pelagic_longline_byScenario.Rdata"))

# Reorder scenarios
scenarios1 <- scenarios %>% 
  mutate(num=1:n()) %>% 
  arrange(EffDyn, PopDyn)
byScen1 <- byScen[scenarios1$num]

# Setup figure
figname <- "AppendixC_simulation_plots.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(4, 3), mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,3,3,3), lwd=0.8)

# Loop through scenarios
for(i in 1:nrow(scenarios1)){
  
  # Subset data
  sdata <- byScen1[[i]]
  
  # Params
  species <- unique(sdata$Species)
  colors <- brewer.pal(length(species), "Set1")
  
  # Plot catch
  ymax <- freeR::ceiling1(max(sdata$Catch), 20)
  plot(Catch ~ Year, sdata, bty="n", las=1, type="n", ylim=c(0,ymax))
  for(j in 1:length(species)){
    sdata1 <- filter(sdata, Species==species[j])
    lines(sdata1$Year, sdata1$Catch, col=colors[j])
  }
  
  # Plot B/BMSY
  ymax <- freeR::ceiling1(max(sdata$BBmsy), 0.5)
  plot(BBmsy ~ Year, sdata, bty="n", las=1, type="n", ylim=c(0,ymax))
  for(j in 1:length(species)){
    sdata1 <- filter(sdata, Species==species[j])
    lines(sdata1$Year, sdata1$BBmsy, col=colors[j])
  }
  abline(h=0.5, lty=3)
  
  # Plot exploitation rate
  ymax <- freeR::ceiling1(max(sdata$ExploitRate), 0.2)
  plot(ExploitRate ~ Year, sdata, bty="n", las=1, type="n", ylim=c(0,ymax))
  for(j in 1:length(species)){
    sdata1 <- filter(sdata, Species==species[j])
    lines(sdata1$Year, sdata1$ExploitRate, col=colors[j])
  }

  
}

# Off
dev.off()
