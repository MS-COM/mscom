
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

# Directories
datadir <- "data/test_cases/original"
outdir <- "data/test_cases/data"
tabledir <- "tables"


# Format data
################################################################################

# Read data
data_orig <- import(file.path(datadir, "Original_A_Catch_1981-2015.xlsx"))

# Format data
data <- data_orig %>% 
  select(-ID) %>% 
  setNames(., gsub(" ", "_", colnames(.))) %>% 
  rename(comm_name=Nombre_Común, sci_name=Scientific_name) %>% 
  melt(id.vars=c("comm_name", "sci_name"), variable.name="year", value.name="catch") %>% 
  mutate(year=as.numeric(as.character(year))) %>% 
  arrange(comm_name, year)

# Check scientific names
freeR::suggest_names(unique(data$sci_name))

# Summarize catch
catch_sum <- data %>% 
  filter(year>2005) %>% 
  group_by(comm_name, sci_name) %>% 
  summarize(c_avg=mean(catch, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(c_perc=c_avg/sum(c_avg)*100,
         name_format=paste0(comm_name, " (", sci_name, ")"),
         c_format=paste0(round(c_avg,1), " / ", round(c_perc, 1), "%")) %>%
  select(name_format, comm_name, sci_name, c_avg, c_perc, c_format) %>% 
  arrange(desc(c_perc))
  
# Export data
write.csv(data, file.path(outdir, "cuban_nearshore_finfish_catch_data.csv"), row.names=F)
write.csv(catch_sum, file.path(tabledir, "TableX_cuba_summary.csv"), row.names=F)



