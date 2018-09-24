

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
datadir <- "data/ram/data"
tabledir <- "data/ram/tables"
codedir <- "code/helper_functions"
outdir <- "data/ram/output"
plotdir <- "data/ram/figures"

# Read data
load(file.path(outdir, "ram_msf_mscmsy_results.Rdata"))

# Read correlations
corrs <- read.csv(file.path(tabledir, "Table1_ram_ms_fisheries_er_correlation.csv"), as.is=T)


# Build data
################################################################################

# Loop through analyzed fisheries
i <- 1
for(i in 1:length(fits)){
  
  # Get output
  fishery <- names(fits[i])
  f_out <- fits[[i]]
  
  # Get predictions
  preds <- f_out$preds
  true <- f_out$data %>% 
    select(stock, year, bbmsy) %>% 
    rename(true=bbmsy)
  
  # Build data
  f_data <- preds %>% 
    left_join(true, by=c("stock", "year")) %>% 
    rename(cmsy=bbmsy_cmsy,
           mscmsy=bbmsy) %>% 
    mutate(fishery=fishery) %>% 
    select(fishery, stock, year, catch, true, cmsy, mscmsy)
  
  # Merge data
  if(i==1){data <- f_data}else{data <- rbind(data, f_data)}
  
}

# Reduce to final year
data1 <- data %>% 
  filter(!is.na(true)) %>% 
  group_by(stock) %>% 
  filter(year==max(year)) %>% 
  mutate(cmsy_diff=(cmsy-true)/true*100,
         mscmsy_diff=(mscmsy-true)/true*100)


# Plot data
################################################################################

# Setup figure
figname <- "Fig2_ram_msf_performance.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=2.25, units="in", res=600)
par(mfrow=c(1,3), mar=c(3.5, 4, 0.5, 0.5), mgp=c(2.5,0.8,0))

# Boxplot percent differenc
boxplot(data1[,c("cmsy_diff", "mscmsy_diff")], frame=F, las=1, lty=1, ylim=c(-100, 400),
        names=c("cMSY", "MS-cMSY"), ylab=expression("Percent error in B/B"["MSY"]*" estimate"), col=c("orange", "darkgreen"))
lines(x=c(0.5,3.5), y=c(0,0), lty=3, col="black")
mtext(paste(nrow(data1), "stocks"), side=3, adj=0.95, line=-1.5, cex=0.7)

# Reset par
par(mgp=c(2,0.8,0))

# cMSY vs. true
plot(cmsy ~ true, data1, bty="n", las=1, pch=16, col="orange", type="n",
     xlim=c(0,4), ylim=c(0,4), xlab=expression("Data-rich B/B"["MSY"]), ylab=expression("Estimated B/B"["MSY"]))
lines(x=c(0,4), y=c(1,1), lty=3)
lines(x=c(1,1), y=c(0,3.5), lty=3)
lines(x=c(0,4), y=c(0,4), lty=2)
lmfit <- lm(cmsy ~ true, data1)
curve(coef(lmfit)[1]+coef(lmfit)[2]*x, from=0, to=4, n=100, add=T, col="orange")
points(x=data1$true, y=data1$cmsy, pch=16, col="orange")
text(x=0, y=3.9, pos=4, labels="cMSY", col="orange", font=2)

# MS-cMSY vs. true
plot(mscmsy ~ true, data1, bty="n", las=1, pch=16, col="darkgreen", type="n",
     xlim=c(0,4), ylim=c(0,4), xlab=expression("Data-rich B/B"["MSY"]), ylab=expression("Estimated B/B"["MSY"]))
lines(x=c(0,4), y=c(1,1), lty=3)
lines(x=c(1,1), y=c(0,3.5), lty=3)
lines(x=c(0,4), y=c(0,4), lty=2)
lmfit <- lm(mscmsy ~ true, data1)
curve(coef(lmfit)[1]+coef(lmfit)[2]*x, from=0, to=4, n=100, add=T, col="darkgreen")
points(x=data1$true, y=data1$mscmsy, pch=16, col="darkgreen")
text(x=0, y=3.9, pos=4, labels="MS-cMSY", col="darkgreen", font=2)

# Off
dev.off()

