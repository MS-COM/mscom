

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
datadir <- "code/ram_stocks/data"
plotdir <- "figures"
tabledir <- "tables"

# Read RAM Legacy Database
load(file.path(datadir, "ram_multispecies_fisheries.Rdata"))


# Plot overlaid effort time series
################################################################################

# Setup figure
figname <- "AppendixA_ram_msf_effort_time_series.pdf"
pdf(paste(plotdir, figname, sep="/"), width=8.5, height=11)
par(mfrow=c(4, 2), mar=c(3.5, 2.5, 2.5, 0.5), mgp=c(2.5,0.8,0), oma=c(3,3,3,3), lwd=0.8)

# Loop through fisheries
for(i in 1:length(data)){
  
  # Subset fishery data
  fishery <- names(data)[i]
  fishery1 <- tolower(gsub("/", "_", gsub(" ", "_", fishery)))
  fdata <- data[[i]]
  erdata_orig <- data[[i]]$er
  
  # Format effort data
  col_use <- apply(erdata_orig, 2, function(x) sum(!is.na(x))) > 0
  erdata <- erdata_orig %>% 
    # Reduce to years where at least 1 stock has effort data
    mutate(any_data=apply(erdata_orig, 1, function(x) sum(!is.na(x)))) %>% 
    filter(any_data>1)
  erdata <- erdata[,col_use] # only use cols with data
  erdata <- select(erdata, -any_data)
  
  # Parameters
  stocks <- colnames(erdata)[2:ncol(erdata)] # I get stock ids here b/c sometime stocks in $info don't have catch
  nstocks <- length(stocks)
  
  # Setup figure
  figname <- paste0("SFig", i, "_effort_ts_", fishery1, ".png")
  png(paste(plotdir, figname, sep="/"), width=5, height=4, units="in", res=600)
  par(mfrow=c(1,1), mar=c(4,2,2,0.5))

  # Plot initial stock
  xmin <- freeR::floor1(min(erdata$year),10)
  xmax <- freeR::ceiling1(max(erdata$year),10)
  colors <- brewer.pal(nstocks, "Set1")
  plot(erdata$year, erdata[,2], type="l", bty="n", 
       xlab="", ylab="Effort", main=fishery,
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
  
  # Add legend
  legend("topleft", bty="n", legend=stocks, col=colors, lwd=1.3)
  
  # Off
  dev.off()
  graphics.off()
  
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
for(i in 1:length(data)){

  # Subset fishery data
  fishery <- names(data)[i]
  fdata <- data[[i]]
  
  # Format ER data
  edata_orig <- fdata$er
  col_use <- apply(edata_orig, 2, function(x) sum(!is.na(x))) > 0
  edata <- select(na.omit(edata_orig[,col_use]), -year)

  # Plot pairwise correlations
  pairs(edata, upper.panel=NULL, las=1, main=fishery)
  
}

# Off
dev.off()
graphics.off()


# Plot correlations
################################################################################

# Loop through fisheries
for(i in 1:length(data)){
  
  # Subset fishery data
  fishery <- names(data)[i]
  fishery1 <- tolower(gsub("/", "_", gsub(" ", "_", fishery)))
  fdata <- data[[i]]
  
  # Format ER data
  edata_orig <- fdata$er
  col_use <- apply(edata_orig, 2, function(x) sum(!is.na(x))) > 0
  edata <- select(na.omit(edata_orig[,col_use]), -year)
  
  # Params
  stocks <- colnames(edata)
  nstocks <- length(stocks)
  
  # Setup figure
  figname <- paste0("SFig", i, "_effort_corr_", fishery1, ".png")
  png(paste(plotdir, figname, sep="/"), width=5, height=5, units="in", res=600)
  par(mfrow=rep(nstocks-1,2), mar=c(2,2,0.5,0.5), oma=c(2,3,1,0))
  
  # Combos
  cells_draw <- t(combn(1:nstocks,2))[,2:1]
  cell_codes_draw <- apply(cells_draw, 1, function(x) paste(x, collapse=""))
  cell_codes_skip <- unique(c(paste0(1, rep(1:nstocks)), paste0(rep(1:nstocks), nstocks)))  # start with 1, end with nstocks

  # Loop through stocks
  for(j in 1:nstocks){ # columns
    for(k in 1:nstocks){ # rows
      code <- paste0(j,k)
      if(code%in%cell_codes){
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

# Calculate mean correlation
################################################################################

# Container
corr_df <- data.frame(fishery=rep(NA, length(data)), nstocks=NA, stocks=NA, nyr=NA, corr_avg=NA)

# Loop through fisheries
for(i in 1:length(data)){

  # Subset fishery data
  fishery <- names(data)[i]
  fdata <- data[[i]]
  erdata_orig <- data[[i]]$er
  tc_use <- data[[i]]$catch_use
  fstocks <- data[[i]]$info

  # Format effort data
  col_use <- apply(erdata_orig, 2, function(x) sum(!is.na(x))) > 0
  erdata <- erdata_orig %>%
    # Reduce to years where at least 1 stock has effort data
    mutate(any_data=apply(erdata_orig, 1, function(x) sum(!is.na(x)))) %>%
    filter(any_data>1)
  erdata <- erdata[,col_use] # only use cols with data
  erdata <- select(erdata, -c(year, any_data))

  # Identify unique combinations of variables
  inds <- 1:ncol(erdata)
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
  corr_df$fishery[i] <- fishery
  corr_df$nstocks[i] <- ncol(erdata)
  corr_df$stocks[i] <- paste(fstocks$comm_name, collapse=", ")
  corr_df$nyr[i] <- nrow(tc_use)
  corr_df$corr_avg[i] <- mean(corrs)

}

# Format table
corr_df1 <- corr_df %>%
  arrange(desc(corr_avg)) %>%
  mutate(stocks1=tolower(paste(nstocks, stocks, sep=" - "))) %>%
  select(fishery, stocks1, nyr, corr_avg)

# Export table
write.csv(corr_df1, file=file.path(tabledir, "Table5_ram_ms_fisheries_overview.csv"), row.names=F)



