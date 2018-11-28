
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
datadir <- "data/simulated/data"
codedir <- "code/helper_functions"
outdir <- "data/output"

# Read MS-cMSY
source(file.path(codedir, "ms_cmsy.R"))


# Format data
################################################################################

# Load scenarios
load(file.path(datadir, "simulated_multispecies_fisheries.Rdata"))

# Loop through files:
i <- 1
# outputs <- NULL
for(i in 1:nrow(scenarios)){
# for(i 1:nrow(scenarios)){
  
  # Container
  print(i)
  # output_i <- NULL

  # Get data
  iter_i <- scenarios$iter_id[i]
  data_i <- filter(sims, iter==iter_i)
  
  # Proceed if data is available
  # Sometimes simulations resulted in negative biomass and were stripped out
  if(nrow(data_i)>0){
    
    # Format data for MS-cMSY
    data <- data_i %>% 
      select(stock, year, catch) %>% 
      mutate(stock=as.character(stock))
    key <- fdata %>% 
      select(comm_name, species, family, resilience) %>% 
      rename(stock=comm_name) %>% 
      mutate(id_fixed=F) %>% 
      filter(stock %in% data$stock)
    
    # Fit MS-cMSY
    out <- ms_cmsy(data, key, npairs=2000)
    # outputs[[i]] <- out$preds
    
    # Setup true
    true <- NULL
    true$ts <- data_i %>% 
      select(stock, year, catch, bbmsy, er)
    true$rk <- fdata %>% 
      filter(comm_name %in% data$stock)
      
    # # Plot MS-cMSY
    # plot_ms_cmsy(out, true)
    
    # Export results
    save(out, file=file.path(outdir, paste0(iter_i, ".Rdata")))
    rm(out)
    
  }
    
}
  
# Export data
# save(outputs, file=file.path(outdir, "simtest_results.Rdata"))

# Identify scenarios which failed
################################################################################

# Which didn't work?
iter_done <- gsub(".Rdata", "", list.files(outdir))
iter_fail <- scenarios$iter_id[!scenarios$iter_id%in%iter_done]
iter_fail_num <- which(scenarios$iter_id%in%iter_fail)



