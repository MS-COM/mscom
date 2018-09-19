
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
datadir <- "code/ram_stocks/data"
tabledir <- "tables"
codedir <- "code/helper_functions"
outdir <- "data/output"

# Read data
load(file.path(datadir, "ram_multispecies_fisheries.Rdata"))

# Read correlation key
corr <- read.csv(file.path(tabledir, "Table5_ram_ms_fisheries_overview.csv"), as.is=T)

# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Format data
################################################################################

# Fisheries to analyse
fisheries <- filter(corr, corr_avg > 0.4)$fishery

# Build key
key_all <- ms_stocks %>% 
  filter(fishery %in% fisheries) %>% 
  select(fishery, stockid, comm_name, species, resilience) %>% 
  mutate(id_fixed=F)

# Add family to key
families <- get_family(key_all$species)
key_all <- key_all %>% 
  left_join(families, by="species") %>% 
  select(fishery, stockid, comm_name, species, family, resilience, id_fixed) %>%
  rename(stock=stockid)


# Analyze data
################################################################################


# Loop through fisheries and analyse 
i <- 3
for(i in 6:length(fisheries)){
  
  # Fishery
  fishery1 <- fisheries[i]
  print(fishery1)
  
  # Format key
  fkey <- filter(key_all, fishery==fishery1)
  
  # Format data
  fdata <- ms_data %>% 
    filter(fishery==fishery1 & !is.na(catch)) %>% 
    select(stockid, year, catch) %>% 
    rename(stock=stockid)
  fkey <- filter(fkey, stock %in% unique(fdata$stock))
  
  # Build true data
  true <- list()
  true_ts <- ms_data %>% 
    filter(fishery==fishery1) %>% 
    select(stockid, year, bbmsy, er) %>% 
    filter(!is.na(bbmsy) | !is.na(er)) %>% 
    rename(stock=stockid)
  true$ts <- true_ts
  
  # Fit model
  out <- ms_cmsy(data=fdata, key=fkey, npairs=200)
  
  # Plot
  plot_ms_cmsy(out)
  
  
}







