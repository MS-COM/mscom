

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
  er_use <- na.omit(erdata)
  
  # Plot if there is data
  if(ncol(erdata)>2){
  
    # Parameters
    stocks <- colnames(erdata)[2:ncol(erdata)] # I get stock ids here b/c sometime stocks in $info don't have catch
    nstocks <- length(stocks)
    
    # # Setup figure
    # figname <- paste0("SFig", i, "_effort_ts_", fishery1, ".png")
    # png(paste(plotdir, figname, sep="/"), width=5, height=4, units="in", res=600)
    # par(mfrow=c(1,1), mar=c(4,2,2,0.5))
  
    # Plot initial stock
    xmin <- freeR::floor1(min(erdata$year),10)
    xmax <- freeR::ceiling1(max(erdata$year),10)
    colors <- brewer.pal(pmax(nstocks,3), "Set1")
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
    
    # Add line marking shared time series
    ymin <- min(erdata[,j], na.rm=T)
    ymax <- max(erdata[,j], na.rm=T)
    yline <- ymin - (ymax-ymin)*0.02
    lines(x=er_use$year, y=rep(yline, nrow(er_use)), lwd=2)
    
    # Add legend
    legend("topleft", bty="n", legend=stocks, col=colors, lwd=1.3)
  
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
for(i in 1:length(data)){

  # Subset fishery data
  fishery <- names(data)[i]
  fdata <- data[[i]]
  
  # Format ER data
  edata_orig <- fdata$er
  col_use <- apply(edata_orig, 2, function(x) sum(!is.na(x))) > 0
  
  # Plot if there is data
  if(sum(col_use)>1){
    
    # Isolate time series
    edata <- select(na.omit(edata_orig[,col_use]), -year)
  
    # Plot pairwise correlations
    pairs(edata, upper.panel=NULL, las=1, main=fishery)
    
  }
  
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
  
  # Plot if there is ER data
  if(sum(col_use)>1){
    
   # Subset data
    edata <- select(na.omit(edata_orig[,col_use]), -year)
  
    # Params
    stocks <- colnames(edata)
    nstocks <- length(stocks)
    
    # Setup figure
    figname <- paste0("SFig", i, "_effort_corr_", fishery1, ".png")
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
corr_df <- data.frame(fishery=rep(NA, length(data)), nstocks=NA, stocks=NA, 
                      nyr=NA, corr_avg=NA, corr_se=NA, vr=NA, vr_scaled=NA)

# Loop through fisheries
for(i in 1:length(data)){

  # Subset fishery data
  fishery <- names(data)[i]
  fdata <- data[[i]]
  erdata_orig <- data[[i]]$er
  tc_use <- data[[i]]$catch_use
  fstocks <- data[[i]]$info

  # Usable data
  col_use <- apply(erdata_orig, 2, function(x) sum(!is.na(x))) > 0
  if(sum(col_use)>=3){
  
    # Format effort data
    erdata <- erdata_orig %>%
      # Reduce to years where at least 1 stock has effort data
      mutate(any_data=apply(erdata_orig, 1, function(x) sum(!is.na(x)))) %>%
      filter(any_data>1)
    erdata <- erdata[,col_use] # only use cols with data
    erdata <- na.omit(select(erdata, -c(year, any_data)))
    erdata_scaled <- t(t(erdata) / apply(erdata, 2, mean))
    # Make sure scaled properly: apply(erdata_scaled, 2, mean)

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
    corr_df$corr_se[i] <- plotrix::std.error(corrs)
    
    # Calculate variance ratio
    corr_df$vr[i] <- var(apply(erdata, 1, sum)) / sum(apply(erdata, 2, var))
    corr_df$vr_scaled[i] <- var(apply(erdata_scaled, 1, sum)) / sum(apply(erdata_scaled, 2, var))
    
  }

}

# Format table
corr_df1 <- corr_df %>%
  filter(!is.na(fishery)) %>% 
  arrange(desc(corr_avg)) %>%
  mutate(stocks1=tolower(paste(nstocks, stocks, sep=" - "))) %>%
  select(fishery, nstocks, stocks1, nyr, corr_avg, vr, vr_scaled)

# Export table
write.csv(corr_df1, file=file.path(tabledir, "Table5_ram_ms_fisheries_overview.csv"), row.names=F)


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
corr_ex <- filter(corr_df1, fishery %in% ex_fisheries)

# Plot correlatations
plot(vr_scaled ~ corr_avg, corr_df1, bty="n", las=1, bg="grey80", pch=21, col="grey40",
     cex=corr_df1$nstocks/min(corr_df1$nstocks), cex.axis=1.1, cex.lab=1.2,
     xlab="Mean correlation", ylab="Variance ratio", xlim=c(-0.2,1), ylim=c(0,4))
lines(x=c(0,0), y=c(0,4), lty=3, col="grey30")
lines(x=c(-0.2,1), y=c(1,1), lty=3, col="grey30")
text(x=corr_ex$corr_avg, y=corr_ex$vr_scaled, labels=corr_ex$fishery, pos=2, offset=0.8)
text(x=1, y=0, labels=paste0(nrow(corr_df1), " fisheries"), adj=1, cex=1.1)

# New par
par(mar=c(3.5, 1.5, 1.0, 6), mgp=c(1,0.8,0), xpd=NA)

# Loop through fisheries
for(i in 1:length(ex_fisheries)){
  
  # Subset fishery data
  fishery <- ex_fisheries[i]
  fdata <- data[[fishery]]
  erdata_orig <- fdata$er
  
  # Format effort data
  col_use <- apply(erdata_orig, 2, function(x) sum(!is.na(x))) > 0
  erdata <- erdata_orig %>% 
    # Reduce to years where at least 1 stock has effort data
    mutate(any_data=apply(erdata_orig, 1, function(x) sum(!is.na(x)))) %>% 
    filter(any_data>1)
  erdata <- erdata[,col_use] # only use cols with data
  erdata <- select(erdata, -any_data)
  er_use <- na.omit(erdata)

  # Parameters
  stocks <- colnames(erdata)[2:ncol(erdata)] # I get stock ids here b/c sometime stocks in $info don't have catch
  nstocks <- length(stocks)
  species <- fdata$info %>% 
    select(stockid, comm_name) %>% 
    filter(stockid %in% stocks) %>% 
    mutate(comm_name=revalue(comm_name, c("Longspine thornyhead"="LS thornyhead",
                                          "Shortspine thornyhead"="SS thornyhead"))) %>% 
    arrange(stocks)
  
  # Plot initial stock
  xmin <- freeR::floor1(min(erdata$year),10)
  xmax <- freeR::ceiling1(max(erdata$year),10)
  colors <- brewer.pal(pmax(nstocks,3), "Set1")
  ylabel <- ifelse(i==2, "Exploitation rate", "")
  plot(erdata$year, erdata[,2], type="l", bty="n", 
       xlab="", ylab=ylabel, main="",
       xaxt="n", yaxt="n", col=colors[1], lwd=1.0, cex.main=0.9, cex.lab=1.2,
       xlim=c(xmin, xmax), ylim=c(min(erdata[,2], na.rm=T), max(erdata[,2], na.rm=T)))
  title(main=fishery, adj=0, cex.main=0.9)
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


