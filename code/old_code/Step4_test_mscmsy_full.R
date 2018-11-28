
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
outdir <- "data/output"

# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Format data
################################################################################

# Load scenarios
load(file.path(datadir, "scenarios_pelagic_longline.Rdata"))

# Files containing each iteration
files <- list.files(datadir, pattern="byScenario")

# Loop through files: 
i <- 1; j <- 1
output_all <- NULL
for(i in 1:length(files)){
  
  # Load file
  file_i <- files[i]
  load(file.path(datadir, file_i))
  
  # Loop through scenarios
  output_i <- NULL
  nscenarios <- length(byScen)
  for(j in 1:nscenarios){
    
    # Scenario data
    sdata <- byScen[[j]]
    print(paste(file_i, j))
    
    # Format data for MS-cMSY
    catch <- dcast(sdata, Year ~ Species, value.var="Catch") %>% 
      rename(year=Year) %>% 
      mutate(year=1960:2009)
    species <- colnames(catch)
    species <- species[species!="year"]
    res <- c("Medium", "Low", "Medium")
    
    # True values
    true <- NULL
    true$r <- arrange(unique(sdata[,c("Species", "r")]), Species)$r
    true$k <- arrange(unique(sdata[,c("Species", "K")]), Species)$K
    true$b_ts <- dcast(sdata, Year ~ Species, value.var="Biomass")
    true$er_ts <- dcast(sdata, Year ~ Species, value.var="ExploitRate")
    true$bbmsy_ts <- dcast(sdata, Year ~ Species, value.var="BBmsy")
    
    # Fit MS-cMSY
    out <- ms_cmsy(catch=catch, stocks=species, res=res, id_fixed=F, npairs=10000)
    output_i[[j]] <- out
    
    # Plot MS-cMSY
    # plot_ms_cmsy(out, true)
    
    
  }
  
  # Combine outputs
  output_all[[i]] <- output_i
  
  
}

# Export data
save(output_all, file=file.path(outdir, "pelagic_longline_simtest.Rdata"))





