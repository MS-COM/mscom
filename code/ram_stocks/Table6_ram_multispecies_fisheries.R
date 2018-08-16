

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "code/ram_stocks/data"
tabledir <- "tables" 

# Load data
load(file.path(datadir, "ram_multispecies_fisheries.Rdata"))


# Build table
################################################################################

# Calculate recent 10 yrs catch
for(i in 1:length(data)){
  
  # Subset data
  sdata <- data[[i]]$catch_use %>% 
    arrange(desc(year)) %>% 
    slice(1:10) %>% 
    arrange(year)
  
  # Stats
  mt <- apply(sdata[,2:ncol(sdata)], 2, mean)
  perc <-  mt / sum( mt )
  stats <- data.frame(stockid=names(mt), mt=mt, perc=perc, stringsAsFactors=F)
  
  # Merge together
  if(i==1){stats1 <- stats}else{stats1 <- rbind(stats1, stats)}
  
}

# Build table
results <- ms_stocks %>% 
  # Add catch stats
  left_join(stats1, by="stockid") %>% 
  # Reduce and arrange
  select(fishery, comm_name, species, stockid, mt, perc, resilience, m) %>% 
  arrange(fishery, desc(perc)) %>% 
  # Fix some names
  mutate(fishery=revalue(fishery, c("Australia SESSF"="Australia SE scalefish and shark",
                                    "NZ Chatham Rise"="NZ Chatham Rise middle-depth")),
         fishery=paste(fishery, "fishery")) %>% 
  # Create name and catch stat columns
  mutate(name=paste0(comm_name, " (", species, ")"),
         catch=paste0(round(mt), " / ", round(perc*100,1), "%"),
         catch=ifelse(catch=="NA / NA%", "no catch data", catch)) %>% 
  # Reduce again
  select(fishery, name, stockid, catch, resilience, m) 


# Export table
################################################################################

write.csv(results, file.path(tabledir, "Table6_ram_multispecies_fisheries.csv"), row.names=F)

