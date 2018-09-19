

# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Turn off sci notation
options(scipen=999)

# Packages
library(freeR)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# Directories
plotdir <- "figures"
outdir <- "data/output"

# Read data
data <- read.csv(file.path(outdir, "pelagic_longline_simtest_results.csv"), as.is=T)

# Format data
################################################################################

# Setup figure
figname <- "FigX_performance_scatterplot.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=3, units="in", res=600)
par(mfrow=c(1,2), mar=c(3.5, 4.5, 0.5, 0.5), mgp=c(2.5,0.8,0))

# Plot all
plot(mscmsy ~ true, data, bty="n", las=1, col="grey60",
     xlim=c(0,2), ylim=c(0,2),
     xlab=expression("B/B"["MSY"]*" true"), 
     ylab=expression("B/B"["MSY"]*" estimated"))
lines(x=c(0.5, 0.5), y=c(0,2), lty=3)
lines(x=c(0, 2), y=c(0.5,0.5), lty=3)

# Plot by ED scenario
plot(mscmsy ~ true, data, bty="n", las=1, col=c("blue", "red")[as.factor(data$ed)],
     xlim=c(0,2), ylim=c(0,2),
     xlab=expression("B/B"["MSY"]*" true"), 
     ylab=expression("B/B"["MSY"]*" estimated"))
lines(x=c(0.5, 0.5), y=c(0,2), lty=3)
lines(x=c(0, 2), y=c(0.5,0.5), lty=3)

# Add legend
legend("bottomleft", bty="n", fill=c("blue", "red"), legend=c("one-way", "two-way"))

# Off
dev.off()

