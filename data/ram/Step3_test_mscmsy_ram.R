
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

# Read correlations
corrs <- read.csv(file.path(tabledir, "Table1_ram_ms_fisheries_er_correlation.csv"), as.is=T)

# Read data
load(file.path(datadir, "ram_multispecies_fisheries.Rdata"))

# Read helper functions
source(file.path(codedir, "ms_cmsy.R"))


# Format data
################################################################################

# Build key
key_all <- stocks %>% 
  select(fishery, stockid, comm_name, species, resilience) %>% 
  mutate(id_fixed=F)
freeR::complete(key_all)

# Add family to key
families <- get_family(key_all$species)
key_all <- key_all %>% 
  left_join(families, by="species") %>% 
  select(fishery, stockid, comm_name, species, family, resilience, id_fixed) %>%
  rename(stock=stockid)


# Analyze data
################################################################################

# Fisheries to analyze
# fisheries <- sort(unique(stocks$fishery))
fisheries <- sort(unique(corrs$fishery[corrs$corr_avg>0.4])) # only highly correlated fisheries
fisheries <- fisheries[!fisheries%in%c("ICCAT East Atlantic tunas", "ICCAT West Atlantic tunas")] # failure to find viable trajectories

# Loop through fisheries and analyse
i <- 1
fits <- list()
for(i in 1:length(fisheries)){
  
  # Fishery
  fishery1 <- fisheries[i]
  print(fishery1)
  
  # Format key
  fkey <- filter(key_all, fishery==fishery1)
  
  # Format data
  fdata <- data %>% 
    filter(fishery==fishery1 & !is.na(catch)) %>% 
    select(stockid, year, catch, bbmsy, er) %>% 
    rename(stock=stockid)
  fkey <- filter(fkey, stock %in% unique(fdata$stock))
  
  # Build true data
  true <- list()
  true_ts <- data %>% 
    filter(fishery==fishery1) %>% 
    select(stockid, year, bbmsy, er) %>% 
    filter(!is.na(bbmsy) | !is.na(er)) %>% 
    rename(stock=stockid)
  true$ts <- true_ts
  
  # Fit model
  output <- ms_cmsy(data=fdata, key=fkey, npairs=1000)
  
  # Plot
  plot_ms_cmsy(output, true)
  
  # Records results
  fits[[i]] <- output
  
}

# Name fits
names(fits) <- fisheries

# Save output
save(fits, stocks, data, file=file.path(outdir, "ram_msf_mscmsy_results.Rdata"))

