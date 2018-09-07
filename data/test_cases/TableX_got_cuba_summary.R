

# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Packages
library(rio)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "data/test_cases/data"
tabledir <- "tables"


# Read data
################################################################################

# Read data
got <- import(file.path(datadir, "gulf_of_thailand_trawl_catch_data.csv"), which=2)
cuba <- import(file.path(datadir, "cuban_nearshore_finfish_catch_data.csv"), which=1)

# Merge data
got$fishery <- "Gulf of Thailand trawl fishery"
cuba$fishery <- "Cuban nearshore finfish fishery"
data <- rbind.fill(got, cuba) %>% 
  select(fishery, everything())


# Summarize data
################################################################################

# Summarize trawl fishery
results <- data %>% 
  filter(year>2005) %>% 
  group_by(fishery, comm_name, sci_name) %>% 
  summarize(c_avg=mean(catch, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(c_perc=c_avg/sum(c_avg)*100,
         name_format=paste0(comm_name, " (", sci_name, ")"),
         c_format=paste0(round(c_avg,1), " / ", round(c_perc, 1), "%")) %>%
  select(fishery, name_format, comm_name, sci_name, c_avg, c_perc, c_format) %>% 
  arrange(desc(fishery), desc(c_perc))

# Export table
write.csv(results, file.path(tabledir, "TableX_real_multispecies_fisheries.csv"), row.names=F)

