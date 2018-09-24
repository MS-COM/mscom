

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
datadir <- "data/ram/data"
plotdir <- "data/ram/figures"
tabledir <- "data/ram/tables"

# Read time series data
load(file.path(datadir, "ram_multispecies_fisheries.Rdata"))

# Read correlation table
corrs <- read.csv(file.path(tabledir, "Table1_ram_ms_fisheries_er_correlation.csv"), as.is=T)


# Plot correlations
################################################################################

# Setup figure
figname <- "Fig1_RAM_MSF_correlations.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=4, units="in", res=600)
layout(matrix(data=c(1,2,
                     1,3,
                     1,4), ncol=2, byrow=T), width=c(0.55, 0.45))
par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.2,0.8,0))

# Plot example fisheries
ex_fisheries <- c("ICCAT East Atlantic tunas", "USA West Coast DTS", "USA GB groundfish")
corr_ex <- filter(corrs, fishery %in% ex_fisheries)

# Plot correlatations
plot(vr_scaled ~ corr_avg, corrs, bty="n", las=1, bg="grey80", pch=21, col="grey40",
     cex=corrs$nstocks/min(corrs$nstocks), cex.axis=1.1, cex.lab=1.2,
     xlab="Mean correlation", ylab="Variance ratio", xlim=c(-0.2,1), ylim=c(0,4))
lines(x=c(0,0), y=c(0,4), lty=3, col="grey30")
lines(x=c(-0.2,1), y=c(1,1), lty=3, col="grey30")
text(x=corr_ex$corr_avg, y=corr_ex$vr_scaled, labels=corr_ex$fishery, pos=2, offset=0.8)
text(x=1, y=0, labels=paste0(nrow(corrs), " fisheries"), adj=1, cex=1.1)

# New par
par(mar=c(3.5, 1.5, 1.0, 6), mgp=c(1,0.8,0), xpd=NA)

# Loop through fisheries
for(i in 1:length(ex_fisheries)){
  
  # Subset fishery data
  fishery1 <- ex_fisheries[i]
  fdata <- filter(data, fishery==fishery1 & !is.na(er))
  
  # Format ER data
  erdata <- dcast(fdata, year ~ stockid, value.var="er")
  
  # Parameters
  fstocks <- colnames(erdata)[2:ncol(erdata)] # I get stock ids here b/c sometime stocks in $info don't have catch
  nstocks <- length(fstocks)
  species <- stocks %>% 
    select(stockid, comm_name) %>% 
    filter(stockid %in% fstocks) %>% 
    mutate(comm_name=revalue(comm_name, c("Longspine thornyhead"="LS thornyhead",
                                          "Shortspine thornyhead"="SS thornyhead"))) %>% 
    arrange(stockid)
  
  # Plot initial stock
  xmin <- freeR::floor1(min(erdata$year),10)
  xmax <- freeR::ceiling1(max(erdata$year),10)
  colors <- brewer.pal(pmax(nstocks,3), "Set1")
  ylabel <- ifelse(i==2, "Exploitation rate", "")
  plot(erdata$year, erdata[,2], type="l", bty="n", 
       xlab="", ylab=ylabel, main="",
       xaxt="n", yaxt="n", col=colors[1], lwd=1.0, cex.main=0.9, cex.lab=1.2,
       xlim=c(xmin, xmax), ylim=c(min(erdata[,2], na.rm=T), max(erdata[,2], na.rm=T)))
  title(main=fishery1, adj=0, cex.main=0.9)
  axis(1, at=seq(xmin, xmax, 10), las=2, cex.axis=1)
  
  # Plot additional stocks
  for(j in 3:ncol(erdata)){
    par(new=T)
    ymin <- min(erdata[,j], na.rm=T)
    ymax <- max(erdata[,j], na.rm=T)
    plot(erdata$year, erdata[,j], type="l", bty="n", 
         xaxt="n", yaxt="n", xlab="", ylab="", col=colors[j-1], lwd=1.0,
         xlim=c(xmin, xmax), ylim=c(ymin, ymax))
  }
  
  # Add legend
  legend(x=xmax-(xmax-xmin)*0.1, y=ymax+(ymax-ymin)*0.2,
         bty="n", legend=species$comm_name, col=colors, lwd=1.1, xpd=NA, cex=0.8)
  
}

# Off
dev.off()


