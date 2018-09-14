
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(freeR)
library(datalimited2)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(quantreg)

# Directories
plotdir <- "figures"

# Read RAM Legacy Database
load("/Users/cfree/Dropbox/Prelim Database Files/Versions/RAM v4.41 (8-20-18)/DB Files With Assessment Data/DBdata (assessment data only).RData")

# Helper functions
################################################################################

# Load helper functions
source("code/ram_stocks/helper_functions.R")

# Calculate cMSY saturation prior
calc_cmsy_s_prior <- function(cr){
  if(cr>=0.8){s_prior <- c(0.4,0.8)}
  if(cr<0.8&cr>=0.5){s_prior <- c(0.2,0.6)}
  if(cr<0.5&cr>=0.35){s_prior <- c(0.01,0.4)}
  if(cr<0.35&cr>=0.15){s_prior <- c(0.01,0.3)}
  if(cr<0.15&cr>=0.05){s_prior <- c(0.01,0.2)}
  if(cr<0.05){s_prior <- c(0.01,0.1)}
  return(s_prior)
}

# Data for plotting polygon
cmsy_priors <- data.frame(cr=seq(0,1,0.001), s_min=NA, s_max=NA)
for(i in 1:nrow(cmsy_priors)){
  cmsy_priors[i,2:3]<- calc_cmsy_s_prior(cmsy_priors$cr[i])
}


# Build data
################################################################################

# Get TB0 and SSB0 values and units
refs <- bioparams_values_views %>% 
  # Get TB0 and SSB0
  select(stockid, TB0, SSB0) %>% 
  filter(!is.na(TB0) | !is.na(SSB0)) %>% 
  rename(tb0=TB0, ssb0=SSB0) %>% 
  # Add BMSY TB0 and SSB0 units
  left_join(select(bioparams_units_views, stockid, TB0, SSB0), by="stockid") %>% 
  rename(tb0_units=TB0, ssb0_units=SSB0)
  
# Get stock time series data
data <- timeseries_values_views %>% 
  # Reduce to stocks with either TB0 or SSB0
  filter(stockid %in% refs$stockid) %>% 
  # Select columns of interest
  select(stockid, year, TB, SSB, TC, TL) %>% 
  setNames(., tolower(colnames(.))) %>% 
  # Add units
  left_join(select(timeseries_units_views, stockid, TB, SSB, TC, TL), by="stockid") %>% 
  rename(tb_units=TB, ssb_units=SSB, tc_units=TC, tl_units=TL) %>% 
  # Add TB0 and SSB0 values
  left_join(refs, by="stockid")  %>% 
  # Create B0 reference points
  mutate(tbdivtb0=tb/tb0,
         ssbdivssb0=ssb/ssb0)

# Inspect data availability
data_n <- data %>% 
  group_by(stockid) %>% 
  summarize(n_tb=sum(!is.na(tb)),
            n_ssb=sum(!is.na(ssb)),
            n_tc=sum(!is.na(tc)),
            n_tl=sum(!is.na(tl))) %>% 
  mutate(b_use=ifelse(n_tb>=n_ssb, "TB", "SSB"),
         c_use=ifelse(n_tc>=n_tl, "TC", "TL"))

# Create catch-to-use and biomass-to-use columns
data1 <- data %>%
  left_join(select(data_n, stockid, c_use, b_use), by="stockid") %>% 
  group_by(stockid) %>% 
  mutate(c=ifelse(c_use=="TC", tc, tl),
         c_units=ifelse(c_use=="TC", tc_units, tl_units),
         b=ifelse(b_use=="TB", tb, ssb),
         b_units=ifelse(b_use=="TB", tb_units, ssb_units), 
         bdivb0=ifelse(b_use=="TB", tbdivtb0, ssbdivssb0))

# Initial saturation
results <- data1 %>% 
  # Only years with both B/B0 and catch
  filter(!is.na(bdivb0) & !is.na(c)) %>%
  # Calculate statistics
  group_by(stockid) %>% 
  summarize(n=n(),
            cmax=max(c),
            cfinal=c[year==max(year)],
            cratio=cfinal/cmax,
            yr1=min(year),
            sat1=bdivb0[year==min(year)],
            sat2=bdivb0[year==max(year)]) %>% 
  # Add TB / SSB
  left_join(refs, by="stockid") %>% 
  mutate(b0=ifelse(!is.na(tb0), tb0, ssb0),
         b0_units=ifelse(!is.na(tb0), tb0_units, ssb0_units), 
         b0divcmax=b0/cmax) %>% 
  # Add species names
  left_join(select(stock, stockid, scientificname, commonname), by="stockid") %>% 
  rename(sci_name=scientificname, comm_name=commonname) %>% 
  mutate(sci_name=revalue(sci_name, c("Chrysophrys auratus"="Pagrus auratus",
                                      "Neoplatycephalus richardsoni"="Platycephalus richardsoni",
                                      "Tetrapturus albidus"="Kajikia albida")),
         comm_name=sentcase(comm_name))
  
# Check species names
check_species(results$sci_name)

# Add resilience values
res <- resilience(results$sci_name)
results <- results %>% 
  left_join(res, by=c("sci_name"="species")) %>% 
  mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")))


# Build data
################################################################################

# Read RAM SP-fits
load("/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output/ramldb_v3.8_sp.Rdata")
lhdata <- read.csv("/Users/cfree/Dropbox/Chris/Rutgers/projects/productivity/models/sst_productivity/output/ramldb_v3.8_spsst_pella_cobe_lme.csv")

# Format r values
rs <- results.wide %>% 
  select(stockid, r) %>% 
  left_join(select(lhdata, stockid, species, comm_name), by="stockid")

# Add resilience
res <- resilience(rs$species)
rs1 <- rs %>% 
  left_join(res, by="species") %>% 
  mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")))


# Plot data
################################################################################

# Setup figure
figname <- "Fig2_prior_development.png"
png(paste(plotdir, figname, sep="/"), width=6.5, height=6.5, units="in", res=600)
# layout(matrix(c(1,1,1,2,2,2,
#                 3,3,3,3,4,4), byrow=T, ncol=6))
layout(matrix(c(1,1,1,1,2,2,
                3,3,3,4,4,4), byrow=T, ncol=6))
par(mar=c(3.5, 3.7, 2, 0.5), mgp=c(2.5,0.8,0))

# C. Intrinsic growth rate
##################################

# Intrinsic rate of growth
boxplot(r ~ resilience, rs1, frame=F, las=1, lty=1, ylim=c(0,2.5), col="grey95", lwd=0.9,
        cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex=1.2, pch=16,
        xlab="Resilience", ylab="r", main="Intrinsic growth rate, r")

# Build r priors
# cMSY
r_priors1 <- matrix(c(0.015, 0.1,
                      0.05, 0.5,
                      0.2, 0.8,
                      0.6, 1.5), ncol=2, byrow=T)
# MS-cMSY
r_priors2 <- matrix(c(0.01, 0.15,
                      0.01, 0.5,
                      0.01, 1.0,
                      0.5, 1.25), ncol=2, byrow=T)

# Add r priors
for(i in 1:nrow(r_priors1)){
  lines(x=rep(-0.5+i, 2), y=c(r_priors1[i,1], r_priors1[i,2]), col=tcolor("orange", 0.5), lwd=3)
  lines(x=rep(-0.45+i, 2), y=c(r_priors2[i,1], r_priors2[i,2]), col=tcolor("darkgreen", 0.3), lwd=3)
}

# Add legend
legend("topleft", legend=c("cMSY", "This study"), bty="n", border=F,
       fill=c(tcolor("orange", 0.5), tcolor("darkgreen", 0.3)), cex=1.3)

# Carrying capacity
##################################

# Carrying capacity
boxplot(results$b0divcmax, log="y", las=1, frame=F, lty=1, col="grey95", lwd=0.9,
        cex.axis=1.2, cex.lab=1.2, cex.main=1.2, cex=1.2, pch=16,
        ylim=c(1,100), ylab=expression("B"["0"]*" / C"["max"]), main="Carrying capacity, K")
# mtext("C", side=3, adj=0.05, line=-1.5, font=2)

# Add priors
lines(x=c(0.7, 0.7), y=c(1,100), col=tcolor("orange", 0.5), lwd=3)
lines(x=c(0.75, 0.75), y=c(2,25), col=tcolor("darkgreen", 0.3), lwd=3)

# Initial saturation
##################################

# Initial saturation
plot(sat1 ~ yr1, results, bty="n", type="n", las=2, ylim=c(0,2),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2,
     xlab="", ylab="Saturation", main="Initial saturation")
# mtext("A", side=3, adj=0.05, line=-1.5, font=2)

# Add MS-cMSY prior
polygon(x=c(1880, 1945, 1980, 2000, 2000, 1980, 1945, 1880),
        y=c(0.8, 0.8, 0.1, 0.1, 1, 1, 1, 1),
        col=tcolor("darkgreen", 0.3), border=F)

# Add cMSY prior
polygon(x=c(1880,1960,1960,2000,2000, 1960,1960, 1880),
        y=c(0.5, 0.5, 0.2, 0.2, 0.6, 0.6, 0.9, 0.9),
        col=tcolor("orange", 0.5), border=F)

# Add points
points(x=results$yr1, y=results$sat1, pch=16, col="grey30", cex=1.2)

# Full saturation line
lines(x=c(1880,2000), y=c(1,1), lwd=1.2, lty=3, col="black")
text("Fully saturated", x=1880, y=2000, pos=2, col="black", cex=0.9, xpd=NA)

# Final saturation
##################################

# Final saturation
plot(sat2 ~ cratio, results, bty="n", type="n", las=1, ylim=c(0,2),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2,
     xlab=expression("C"["final"]*" / C"["max"]), ylab="Saturation", main="Final saturation")
# mtext("B", side=3, adj=0.05, line=-1.5, font=2)

# Add MS-cMSY prior
x <- seq(0,1,0.1)
y1 <- 0 + 0.4*x
# y2 <- 0.8 + 0.2*x
y2 <- 0.5 + 0.4*x
polygon(x=c(x, rev(x)), y=c(y2, rev(y1)), col=tcolor("darkgreen", 0.3), border=F)

# Add cMSY prior
x <- cmsy_priors$cr
y2 <- cmsy_priors$s_min
y1 <- cmsy_priors$s_max
polygon(x=c(x, rev(x)), y=c(y1, rev(y2)), col=tcolor("orange", 0.5), border=F)

# Add points
points(x=results$cratio, y=results$sat2, pch=16, col="grey30", cex=1.2)

# Full saturation line
lines(x=c(0,1), y=c(1,1), lwd=1.2, lty=3, col="black")
text("Fully saturated", x=1.0, y=1.08, pos=2, col="black", cex=1.1, xpd=NA)

# Off
dev.off()
graphics.off()



