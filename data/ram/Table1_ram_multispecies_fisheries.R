

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "data/ram/data"
tabledir <- "data/ram/tables" 

# Load data
load(file.path(datadir, "ram_multispecies_fisheries.Rdata"))


# Build table
################################################################################

# Build catch statistics
# Mean catch and proportion of total over 
# the most recent 5 years of overlapping catch time series
fisheries <- sort(unique(stocks$fishery))
for(i in 1:length(fisheries)){
  fishery1 <- fisheries[i]
  cdata <- data %>% 
    filter(fishery==fishery1) %>% 
    dcast(year ~ stockid, value.var="catch") %>% 
    na.omit() %>% 
    arrange(desc(year)) %>% 
    slice(1:5) %>% 
    arrange(year) %>% 
    select(-year)
  mt <- apply(cdata, 2, mean)
  perc <- mt / sum( mt )
  stats <- data.frame(stockid=names(mt), mt=mt, perc=perc)
  if(i==1){cstats <- stats}else{cstats <- rbind(cstats, stats)}
}

# ER stats
erstats <- data %>%
  group_by(stockid) %>% 
  summarize(er=ifelse(sum(!is.na(er)>0), "Yes", "---"))

# Build table
results <- stocks %>% 
  # Add catch stats
  left_join(cstats, by="stockid") %>% 
  # Add ER stats
  left_join(erstats, by="stockid") %>% 
  # Reduce and arrange
  select(fishery, comm_name, species, stockid, mt, perc, er, family, resilience) %>% 
  arrange(fishery, desc(perc)) %>% 
  # Fix some names
  mutate(fishery=paste(fishery, "fishery")) %>% 
  # Create name and catch stat columns
  mutate(name=paste0(comm_name, " (", species, ")"),
         er=ifelse(is.na(er), "---", er),
         catch=paste0(round(mt), " / ", round(perc*100,1), "%"),
         catch=ifelse(catch=="NA / NA%", "no catch data", catch)) %>% 
  # Reduce again
  select(fishery, name, stockid, catch, er, family, resilience) 


# Export table
################################################################################

write.csv(results, file.path(tabledir, "Table1_ram_ms_stocks.csv"), row.names=F)

