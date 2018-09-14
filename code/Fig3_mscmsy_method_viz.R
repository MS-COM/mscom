

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

# Directories
outdir <- "data/output"
simdir <- "data/simulated"
plotdir <- "figures"
tabledir <- "tables"

# Read RAM Legacy Database
load(file.path(outdir, "pelagic_longline_simtest.Rdata"))

# Load scenarios/simulations
load(file.path(simdir, "scenarios_pelagic_longline.Rdata"))
load(file.path(simdir, "sim1_pelagic_longline_byScenario.Rdata"))


# Pick scenario and build "true" data
################################################################################

# Which scenario?
s <- 4
sdata <- byScen[[s]]
out <- output_all[[1]][s]

# True values
true <- NULL
true$r <- arrange(unique(sdata[,c("Species", "r")]), Species)$r
true$k <- arrange(unique(sdata[,c("Species", "K")]), Species)$K
true$b_ts <- dcast(sdata, Year ~ Species, value.var="Biomass")
true$er_ts <- dcast(sdata, Year ~ Species, value.var="ExploitRate")
true$bbmsy_ts <- dcast(sdata, Year ~ Species, value.var="BBmsy")

# Plot data
################################################################################

# Setup figure
figname <- "Fig1_RAM_MSF_correlations.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=4, units="in", res=600)

# Extract predictions
results <- out[[1]]
b <- results$B_ts
u <- results$Upred_ts
bbmsy <- results$BBmsy_ts
uumsy <- results$UUmsy_qE_ts

# Params
xmin <- freeR::floor1(min(years), 10)
xmax <- freeR::ceiling1(max(years), 10)

# Extract true values
true_yn <- T
if(true_yn){
  b_true <- true$B_ts
  u_true <- true$U_ts
  bbmsy_true <- true$BBMSY_ts
  uumsy_true <- true$UUMSY_ts
}

# Species and colors
nyr  <- nrow(b)
nspp <- ncol(b)
colors <- RColorBrewer::brewer.pal(nspp, "Set1")
par(mfrow=c(2,2), mar=c(2,4,0.5,0.5), xpd=F)

# A. Plot biomass
ymax <- ifelse(true_yn, max(c(b, b_true[,2:ncol(b_true)]), na.rm=T), max(b))
plot(1, 1, type="n", bty="n", las=1, 
     xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab="", ylab="Biomass")
for(i in 1:nspp){lines(x=years, y=b[,i], col=colors[i])}
if(true_yn){
  for(i in 1:nspp){
    lines(x=b_true[,1], y=b_true[,i+1], col=colors[i], lty=2)
  }
}

# B. Plot exploitation rate
ymax <- ifelse(true_yn, max(c(u, u_true[,2:ncol(u_true)]), na.rm=T), max(u))
plot(1, 1, type="n", bty="n", las=1, 
     xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab="", ylab="Exploitation rate")
for(i in 1:nspp){lines(x=years, y=u[,i], col=colors[i])}
if(true_yn){
  for(i in 1:nspp){
    lines(x=u_true[,1], y=u_true[,i+1], col=colors[i], lty=2)
  }
}

# C. Plot B/BMSY
ymax <- ifelse(true_yn, max(c(bbmsy, bbmsy_true[,2:ncol(bbmsy_true)]), na.rm=T), max(bbmsy))
plot(1, 1, type="n", bty="n", las=1, 
     xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab="", ylab=expression("B/B"["MSY"]))
for(i in 1:nspp){lines(x=years, y=bbmsy[,i], col=colors[i])}
if(true_yn){
  for(i in 1:nspp){lines(x=bbmsy_true[,1], y=bbmsy_true[,i+1], col=colors[i], lty=2)}
}

# D. Plot U/UMSY
ymax <- ifelse(true_yn, max(c(uumsy, uumsy_true[,2:ncol(uumsy_true)]), na.rm=T), max(uumsy))
plot(1, 1, type="n", bty="n", las=1, 
     xlim=c(xmin, xmax), ylim=c(0, ymax),
     xlab="", ylab=expression("U/U"["MSY"]))
for(i in 1:nspp){lines(x=years, y=uumsy[,i], col=colors[i])}
if(true_yn){
  for(i in 1:nspp){lines(x=uumsy_true[,1], y=uumsy_true[,i+1], col=colors[i], lty=2)}
}
legend("topleft", bty="n", col=colors, lwd=1.2, legend=stocks)



# Off
dev.off()


