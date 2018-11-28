

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
data <- read.csv(file.path(outdir, "simtest_results.csv"), as.is=T)

# Filter data
data1 <- filter(data, year==max(year))

# Format data
################################################################################

# Setup figure
figname <- "Fig4_performance_overall.png"
png(paste(plotdir, figname, sep="/"), width=4.5, height=4.5, units="in", res=600)
par(mar=c(3.5, 4.5, 0.5, 0.5), mgp=c(3.3,1,0))

# Plot
boxplot(data1[,c("cmsy_diff", "mscmsy_diff")], frame=F, las=1, lty=1, ylim=c(-100, 400),
        col=c(tcolor("orange", 0.5), tcolor("darkgreen", 0.3)),
        names=c("cMSY", "MS-cMSY"), ylab=expression("Percent error in B/B"["MSY"]*" estimate"))
lines(x=c(0.5,3.5), y=c(0,0), lty=2, col="black")

# Off
dev.off()

