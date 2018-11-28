

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

# Scenarios
scenarios <- data %>% 
  select(iter, ed, id, sigmaP, sigmaF, nstocks, nyrs) %>% 
  unique()

# Reshape for boxplots
data1 <- data %>% 
  filter(year==2015) %>% 
  select(iter, cmsy_diff, mscmsy_diff) %>% 
  rename(cmsy=cmsy_diff, mscmsy=mscmsy_diff)
data1_long <- melt(data1, id.vars=c("iter"), variable.name="model", value.name="diff")
data2 <- data1_long %>% 
  left_join(scenarios, by="iter") %>% 
  mutate(model=factor(model, levels=c("cmsy", "mscmsy")),
         ed=factor(ed, levels=c("one-way", "two-way")),
         id=factor(id, levels=c("unfished", "fished")),
         sigmaP=factor(sigmaP),
         sigmaF=factor(sigmaF),
         nstocks=factor(nstocks),
         nyrs=factor(nyrs))


# Format data
################################################################################

# Colors
cmsy_col <- tcolor("orange", 0.5)
mscmsy_col <- tcolor("darkgreen", 0.3)
cols <- c(cmsy_col, mscmsy_col)

# Factor info  
sfs <- c("ed", "id", "sigmaP", "sigmaF", "nstocks", "nyrs")
sf_titles <- c("Effort dynamics", "Initial saturation", "Process variability", "Effort correlation", "# of stocks", "# of years")


# Setup figure
figname <- "Fig5_performance_by_scenario.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=5, units="in", res=600)
par(mfrow=c(2,3), mar=c(2.5, 2.5, 1.5, 0.5), mgp=c(3.3,1,0), oma=c(0,3,0,0))

# Loop through factors
for(i in 1:length(sfs)){
  
  # Factor info
  sf <- sfs[i]
  n_levels <- nlevels(data2[,sf])
  
  # Plotting position
  ats1 <- c(0.5, 0.5 + 2.5 * 1:n_levels)[1:n_levels]
  ats2 <- ats1 + 1
  ats <- sort(c(ats1, ats2))
  
  # Calculate midpoints
  mpoints <- apply(matrix(ats, ncol=2, byrow=T), 1, mean)
  
  # Plot data
  boxplot(data2$diff ~ data2$model + data2[,sf], data2, las=2, frame=F, lty=1, at=ats,
          ylim=c(-100, 400), col=cols, main=sf_titles[i], xaxt="n")
  lines(x=c(0,max(ats)+0.5), y=c(0,0), lty=2, col="black")
  axis(side=1, at=mpoints, labels=levels(data2[,sf]))
  
}

# Add y-axis label
ylab <- expression("Percent error in B/B"["MSY"]*" estimate")
mtext(ylab, outer=T, side=2, adj=0.5, line=0.5)

# Off
dev.off()
