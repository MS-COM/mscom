
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
datadir <- "data/simulated"
codedir <- "code/helper_functions"

# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Format data
################################################################################

# Load scenarios
load(file.path(datadir, "scenarios_pelagic_longline.Rdata"))

# Files containing each iteration
files <- list.files(datadir, pattern="byScenario")

# Loop through files: i <- 1; j <- 1
for(i in 1:length(files)){
  
  # Load file
  file_i <- files[i]
  load(file.path(datadir, file_i))
  
  # Loop through scenarios
  nscenarios <- length(byScen)
  for(j in 1:nscenarios){
    
    # Scenario data
    sdata <- byScen[[j]]
    
    # Format data for MS-cMSY
    catch <- dcast(sdata, Year ~ Species, value.var="Catch")
    yrs <- catch$Year
    catch <- as.matrix(select(catch, -Year))
    species <- colnames(catch)
    res <- c("Medium", "Low", "Medium")
    
    # billfish, shark, tuna
    species <- colnames(catch)
    r_true <- c(0.67, 0.11, 0.60)
    k_true <- c(227, 738, 304)
    res <- c("Medium", "Low", "Medium")
    yrs <- 1:nrow(catch)
    
    # True
    true <- NULL
    true$r <- arrange(unique(sdata[,c("Species", "r")]), Species)$r
    true$k <- arrange(unique(sdata[,c("Species", "K")]), Species)$K
    true$b_ts <- dcast(sdata, Year ~ Species, value.var="Biomass")
    true$er_ts <- dcast(sdata, Year ~ Species, value.var="ExploitRate")
    true$bbmsy_ts <- dcast(sdata, Year ~ Species, value.var="BBmsy")
    
    # Fit MS-cMSY
    out <- fit_mssra(catch=catch, years=yrs, stocks=species, res=res, id_fixed=F)
    
    # Plot MS-cMSY
    plot_mssra(out, true)
    
    
    
  }
  
  
}





