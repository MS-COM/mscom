

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

# Read RAM Legacy Database
load(file.path(datadir, "ram_multispecies_fisheries.Rdata"))


# Plot overlaid effort time series
################################################################################

# Fisheries
fisheries <- sort(unique(stocks$fishery))

# Setup figure
figname <- "AppendixA_ram_msf_effort_time_series.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(4, 2), mar=c(3.5, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,3,3,3), lwd=0.8)

# Loop through fisheries
for(i in 1:length(fisheries)){
  
  # Subset fishery data
  fishery1 <- fisheries[i]
  fishery2 <- tolower(gsub("/", "_", gsub(" ", "_", fishery1)))
  
  # Build ER fata
  fdata <- data %>%
    filter(fishery==fishery1 & !is.na(er)) 
  
  # If there is ER data available...
  if(nrow(fdata)>0){
    
    # Format ER data
    erdata <- dcast(fdata, year ~ stockid, value.var="er")
    erdata_use <- na.omit(erdata)
    
    # Stocks to plot
    fstocks <- sort(unique(fdata$stockid))
    nstocks <- length(fstocks)
      
    # # Setup figure
    # figname <- paste0("SFig", i, "_effort_ts_", fishery1, ".png")
    # png(paste(plotdir, figname, sep="/"), width=5, height=4, units="in", res=600)
    # par(mfrow=c(1,1), mar=c(4,2,2,0.5))
    
    # Plot initial stock
    xmin <- freeR::floor1(min(erdata$year),10)
    xmax <- freeR::ceiling1(max(erdata$year),10)
    colors <- brewer.pal(pmax(nstocks,3), "Set1")
    plot(erdata$year, erdata[,2], type="l", bty="n", 
         xlab="", ylab="Effort", main=fishery1,
         xaxt="n", yaxt="n", col=colors[1], lwd=1.3, cex.main=1.5,
         xlim=c(xmin, xmax), ylim=c(min(erdata[,2], na.rm=T), max(erdata[,2], na.rm=T)))
    axis(1, at=seq(xmin, xmax, 10), las=2, cex.axis=1.3)
    
    # Plot additional stocks
    for(j in 3:ncol(erdata)){
      par(new=T)
      plot(erdata$year, erdata[,j], type="l", bty="n", 
           xaxt="n", yaxt="n", xlab="", ylab="", col=colors[j-1], lwd=1.3,
           xlim=c(xmin, xmax), ylim=c(min(erdata[,j], na.rm=T), max(erdata[,j], na.rm=T)))
    }
    
    # Add line marking shared time series
    ymin <- min(erdata[,j], na.rm=T)
    ymax <- max(erdata[,j], na.rm=T)
    yline <- ymin - (ymax-ymin)*0.02
    lines(x=erdata_use$year, y=rep(yline, nrow(erdata_use)), lwd=2)
    
    # Add legend
    legend("topleft", bty="n", legend=fstocks, col=colors, lwd=1.3)
    
    # # Off
    # dev.off()
    # graphics.off()
      
  }
    
  
}

# Off
dev.off()
graphics.off()


# Plot pairwise correlations
################################################################################

# Setup figure
figname <- "AppendixB_ram_msf_effort_corr_matrix.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(2, 1), mar=c(3.5, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,3,3,3), lwd=0.8)

# Loop through fisheries
for(i in 1:length(fisheries)){

  # Subset fishery data
  fishery1 <- fisheries[i]
  fishery2 <- tolower(gsub("/", "_", gsub(" ", "_", fishery1)))
  
  # Build ER fata
  fdata <- data %>%
    filter(fishery==fishery1 & !is.na(er)) 
  
  # If there is ER data available...
  if(nrow(fdata)>0){
    
    # Format ER data
    erdata <- dcast(fdata, year ~ stockid, value.var="er")
    erdata_use <- na.omit(erdata)
    
    # Isolate time series
    edata <- select(erdata_use, -year)
  
    # Plot pairwise correlations
    pairs(edata, upper.panel=NULL, las=1, main=fishery1)
    
  }
  
}

# Off
dev.off()
graphics.off()


# Plot correlations
################################################################################

# Loop through fisheries
for(i in 1:length(fisheries)){
  
  # Subset fishery data
  fishery1 <- fisheries[i]
  fishery2 <- tolower(gsub("/", "_", gsub(" ", "_", fishery1)))
  
  # Build ER fata
  fdata <- data %>%
    filter(fishery==fishery1 & !is.na(er)) 
  
  # If there is ER data available...
  if(nrow(fdata)>0){
    
    # Format ER data
    erdata <- dcast(fdata, year ~ stockid, value.var="er")
    erdata_use <- na.omit(erdata)
    
    # Isolate time series
    edata <- select(erdata_use, -year)
  
    # Params
    stocks <- colnames(edata)
    nstocks <- length(stocks)
    
    # Setup figure
    figname <- paste0("SFig", i, "_effort_corr_", fishery2, ".png")
    png(file.path(plotdir, "ram_correlations", figname), width=5, height=5, units="in", res=600)
    par(mfrow=rep(nstocks-1,2), mar=c(2,2,0.5,0.5), oma=c(2,3,1,0))
    
    # Combos
    cells_draw <- matrix(t(combn(1:nstocks,2))[,2:1], ncol=2)
    cell_codes_draw <- apply(cells_draw, 1, function(x) paste(x, collapse=""))
    cell_codes_skip <- unique(c(paste0(1, rep(1:nstocks)), paste0(rep(1:nstocks), nstocks)))  # start with 1, end with nstocks
  
    # Loop through stocks
    for(j in 1:nstocks){ # columns
      for(k in 1:nstocks){ # rows
        code <- paste0(j,k)
        if(code%in%cell_codes_draw){
          xmax <- freeR::ceiling1(max(edata[,k]),0.1)
          ymax <- freeR::ceiling1(max(edata[,j]),0.1)
          xlabel <- ifelse(j==nstocks, stocks[k], "") #"" # stocks[k]
          ylabel <- ifelse(k==1, stocks[j], "")
          plot(edata[,j] ~ edata[,k], bty="n", las=1, 
               xlim=c(0, xmax), ylim=c(0,ymax), cex.axis=0.9,
               xlab=xlabel, ylab=ylabel, main="", col="grey50", xpd=NA)
          lmfit <- lm(edata[,j] ~ edata[,k])
          curve(coef(lmfit)[1]+coef(lmfit)[2]*x, from=0, to=xmax, n=50, add=T)
          r2 <- freeR::r2(lmfit)
          corr <- cor(edata[,j], edata[,k])
          text(x=xmax, y=0+ymax*0.05, label=paste0("r=",freeR::roundf(corr, 2)), offset=0.01, pos=2, xpd=NA)
        }else{
          if(!code%in%cell_codes_skip){
            plot(1:10, 1:10, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", main="")
          }
        }
      }
    }
    
    # Add title
    mtext(text=fishery, outer=T, side=3, adj=0.5, line=-0.5, font=2, xpd=NA)
    
    # Off
    dev.off()
    graphics.off()
      
  }
  
}


# Calculate mean correlation
################################################################################

# Container
corr_df <- data.frame(fishery=fisheries, nstocks=NA, stocks=NA, 
                      nyr=NA, corr_avg=NA, corr_se=NA, vr=NA, vr_scaled=NA)

# Loop through fisheries
for(i in 1:length(fisheries)){

  # Subset fishery data
  fishery1 <- fisheries[i]
  fdata <- filter(data, fishery==fishery1)
  
  # Build ER fata
  fdata <- data %>%
    filter(fishery==fishery1 & !is.na(er)) 
  
  # If there is ER data available...
  if(nrow(fdata)>0){
    
    # Format ER data
    erdata_all <- dcast(fdata, year ~ stockid, value.var="er")
    erdata_use <- na.omit(erdata_all)
    erdata <- select(erdata_use, -year)
  
    # Scale ER data
    erdata_scaled <- t(t(erdata) / apply(erdata, 2, mean))
    # Make sure scaled properly: apply(erdata_scaled, 2, mean)
    
    # Params
    stocks <- colnames(erdata)
    nstocks <- length(stocks)

    # Identify unique combinations of variables
    inds <- 1:nstocks
    combos <- t(combn(inds, 2))
  
    # Loop through combos and calculate correlation
    corrs <- NA
    for(j in 1:nrow(combos)){
  
      # Subset data
      cols <- combos[j,]
      sdata <- as.matrix(na.omit(erdata[,cols]))
  
      # Calculate and record correlation
      corrs[j] <- cor(x=sdata[,1], y=sdata[,2])
  
    }
    
    # Record info
    corr_df$fishery[i] <- fishery1
    corr_df$nstocks[i] <- nstocks
    corr_df$stocks[i] <- paste(stocks, collapse=", ")
    corr_df$nyr[i] <- nrow(erdata)
    corr_df$corr_avg[i] <- mean(corrs)
    corr_df$corr_se[i] <- plotrix::std.error(corrs)
    
    # Calculate variance ratio
    corr_df$vr[i] <- var(apply(erdata, 1, sum)) / sum(apply(erdata, 2, var))
    corr_df$vr_scaled[i] <- var(apply(erdata_scaled, 1, sum)) / sum(apply(erdata_scaled, 2, var))
    
  }

}

# Format table
corr_df1 <- corr_df %>%
  filter(!is.na(nstocks)) %>% 
  arrange(desc(corr_avg)) %>%
  mutate(stocks1=paste(nstocks, stocks, sep=" - ")) %>%
  select(fishery, nstocks, stocks1, nyr, corr_avg, vr, vr_scaled)

# Export table
write.csv(corr_df1, file=file.path(tabledir, "Table1_ram_ms_fisheries_er_correlation.csv"), row.names=F)

