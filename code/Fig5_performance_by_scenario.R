

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

# Scenarios
scenarios <- data %>% 
  select(sim_id, dyn, id, ed, maxf, effsd, iter) %>% 
  unique()

# Reshape for boxplots
data1 <- data %>% 
  select(stock_id, sim_id, cmsy_diff, mscmsy_diff) %>% 
  rename(cmsy=cmsy_diff, mscmsy=mscmsy_diff)
data1_long <- melt(data1, id.vars=c("sim_id", "stock_id"), variable.name="model", value.name="diff")
data2 <- data1_long %>% 
  left_join(scenarios, by="sim_id") %>% 
  select(sim_id, dyn:iter, everything()) %>% 
  mutate(model=factor(model, levels=c("cmsy", "mscmsy")),
         dyn=factor(dyn, levels=c("deterministic", "variable")),
         ed=factor(ed, levels=c("one-way", "two-way")),
         id=factor(id, levels=c("unfished", "fished")))


# Format data
################################################################################

# Colors
cmsy_col <- tcolor("orange", 0.5)
mscmsy_col <- tcolor("darkgreen", 0.3)
cols <- rep(c(cmsy_col, mscmsy_col), 2)
  
# Setup figure
figname <- "Fig5_performance_by_scenario.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6.5, units="in", res=600)
par(mfrow=c(2,2), mar=c(2.5, 2.5, 1.5, 0.5), mgp=c(3.3,1,0), oma=c(0,3,0,0))

# Initial saturation
boxplot(diff ~ model + dyn, data2, las=2, frame=F, lty=1, at=c(0.5, 1.5, 3, 4), 
        col=cols, ylim=c(-200, 400), xaxt="n", main="Initial saturation")
lines(x=c(0,4.5), y=c(0,0), lty=2, col="black")
axis(1, at=c(1,3.5), labels=c("Unfished", "Fished"), las=1)

# Exploitation dynamics
boxplot(diff ~ model + ed, data2, las=2, frame=F, lty=1, at=c(0.5, 1.5, 3, 4),
        col=cols, ylim=c(-200, 400), xaxt="n", main="Exploitation dynamics")
lines(x=c(0,4.5), y=c(0,0), lty=2, col="black")
axis(1, at=c(1,3.5), labels=c("One-way", "Two-way"), las=1)
legend("topright", fill=c(cmsy_col, mscmsy_col), legend=c("cMSY", "MS-cMSY"), bty="n")

# Process error
boxplot(diff ~ model + dyn, data2, las=2, frame=F, lty=1, at=c(0.5, 1.5, 3, 4),
        col=cols, ylim=c(-200, 400), xaxt="n", main="Process variability")
lines(x=c(0,4.5), y=c(0,0), lty=2, col="black")
axis(1, at=c(1,3.5), labels=c("Deterministic", "Variable"), las=1)

# Effort correlation
boxplot(diff ~ model + effsd, data2, las=2, frame=F, lty=1, at=c(0.5, 1.5, 3, 4),
        col=cols, ylim=c(-200, 400), xaxt="n", main="Effort correlation")
lines(x=c(0,4.5), y=c(0,0), lty=2, col="black")
axis(1, at=c(1,3.5), labels=c("Perfect", "Imperfect"), las=1)

# Add y-axis label
ylab <- expression("Percent error in B/B"["MSY"]*" estimate")
mtext(ylab, outer=T, side=2, adj=0.5, line=0.5)

# Off
dev.off()
