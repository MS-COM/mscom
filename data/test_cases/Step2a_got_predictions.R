
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Turn off sci notation
options(scipen=999)

# Packages
library(rio)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# Datalimited2
# devtools::install_github("cfree14/datalimited2")
# library(datalimited2)

# Directories
datadir <- "data/test_cases/data"
plotdir <- "figures"
codedir <- "code/helper_functions"

# Read data
data_orig <- import(file.path(datadir, "gulf_of_thailand_trawl_catch_data.csv"))

# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Make predictions
################################################################################

# Look up resilience
# spp <- unique(data$sci_name)
# res <- datalimited2::resilience(spp)

# Format data
data <- data_orig %>% 
  rename(stock=comm_name) %>% 
  select(stock, year, catch)

# Build true
true <- list()
true$ts <- data_orig %>% 
  rename(stock=comm_name)

# Build key
key <- data.frame(stock=c("Bigeye scad", "Blue swimming crab", "Largehead hairtail"),
                  family=c("Carangidae", "Portunidae", "Trichiuridae"),
                  resilience=c("High", "High", "Low"),
                  id_fixed=F,
                  stringsAsFactors=F)

# Make predictions
out <- ms_cmsy(data=data, key=key, npairs=5000)

# Plot data
plot_ms_cmsy(out)

# Export predictions
save(out, file=file.path(datadir, "got_ms_cmsy_predictions.Rdata"))



