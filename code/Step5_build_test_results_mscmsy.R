
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


# Format data
################################################################################

# Load MS-cMSY results
load(file.path(outdir, "pelagic_longline_simtest.Rdata"))

# Load simulation scenarios
load(file.path(datadir, "scenarios_pelagic_longline.Rdata"))
scenarios <- scenarios %>% 
  rename(dyn=PopDyn, id=DataStart, ed=EffDyn, maxf=MaxF, nyr=Nyears, effsd=EffSD)

# Simulation output: 1 file per iteration
# These were run through MS-cMSY in the same order
files <- list.files(datadir, pattern="byScenario")

# Loop through iterations: i <- 1; j <- 1
for(i in 1:length(files)){
  
  # MS-cMSY predictions
  preds_i <- output_all[[1]]
  
  # Load simulation output (to get true B/BMSY values)
  file_i <- files[i]
  load(file.path(datadir, file_i))
  
  # Loop through scenarios
  nscenarios <- length(byScen) # 1 element per scenario
  for(j in 1:nscenarios){
    
    # Observed status
    true <- byScen[[j]] %>% 
      filter(Year==max(Year)) %>% 
      select(Species, BBmsy) %>% 
      rename(stock=Species, true=BBmsy)
    
    # MS-cMSY status
    mscmsy <- sapply(preds_i[[j]]$bbmsy_vv_median, cbind)
    colnames(mscmsy) <- preds_i[[j]]$stocks
    mscmsy_end <- mscmsy[nrow(mscmsy),]
    mscmsy_end <- data.frame(stock=colnames(mscmsy), mscmsy=mscmsy_end)
    
    # Add scenario meta-data
    sdata <- true %>%
      left_join(mscmsy_end, by="stock")
    sdata1 <- cbind(scenarios[j,], sdata) # Raises WARNINGS that don't matter
    
    # Merge data
    if(j==1){data_i <- sdata1}else{data_i <- rbind(data_i, sdata1)}
    
  }
  
  # Format data 
  fdata <- data_i %>% 
    mutate(simfile=file_i, iter=i) %>% 
    select(simfile, dyn:effsd, iter, everything())
  
  # Merge
  if(i==1){data <- fdata}else{data <- rbind(data, fdata)}
  
}

# Format
data <- data %>% 
  arrange(dyn, id, ed, maxf, effsd, iter)


