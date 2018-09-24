
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Packages
library(plyr)
library(dplyr)
library(reshape2)

# Directories
datadir <- "data/ram/data"
codedir <- "code/ram_stocks"

# Read RAM multispecies fisheries
ms_key <- read.csv(file.path(datadir, "RAM_MS_fisheries.csv"), as.is=T)

# Read RAM Legacy Database
load("/Users/cfree/Dropbox/Prelim Database Files/Versions/RAM v4.41 (8-20-18)/DB Files With Assessment Data/DBdata (assessment data only).Rdata")


# Build data
################################################################################

# Are all of the stockids in RAM?
ms_stockids <- ms_key$stockid
ms_stockids[!ms_stockids %in% stock$stockid]

# Build stock key
ms_stocks <- stock %>% 
  select(-c(tsn, inmyersdb, myersstockid)) %>% 
  # Add area name
  left_join(select(area, areaid, country, areaname), by="areaid") %>% 
  # Rename columns
  rename(species=scientificname, comm_name=commonname) %>% 
  # Filter to multi-species stocks
  filter(stockid %in% ms_stockids) %>% 
  # Correct some species names
  mutate(species=revalue(species, c("Neoplatycephalus richardsoni"="Platycephalus richardsoni",
                                    "Tetrapturus albidus"="Kajikia albida")),
         comm_name=freeR::sentcase(comm_name)) %>% 
  # Add fishery name
  left_join(select(ms_key, c(stockid, fishery)), by="stockid") %>% 
  # Add catch/biomass units
  left_join(select(timeseries_ids_views, stockid, TCbest, TBbest), by="stockid") %>% 
  rename(tc_units=TCbest, tb_units=TBbest) %>% 
  # Add reference points
  left_join(select(bioparams_values_views, stockid, MSY, TBmsybest, Fmsy), by="stockid") %>% 
  rename(msy=MSY, bmsy=TBmsybest, fmsy=Fmsy) %>% 
  left_join(select(bioparams_ids_views, stockid, MSY, TBmsybest, Fmsy), by="stockid") %>% 
  rename(msy_units=MSY, bmsy_units=TBmsybest, fmsy_units=Fmsy) %>% 
  # Rearrange 
  select(fishery, everything()) %>% 
  arrange(fishery, stockid)

# Check scientific names
freeR::check_names(ms_stocks$species)

# Remove fisheries with fewer than 2 constituent stocks
nstocks <- ms_stocks %>%
  group_by(fishery) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  filter(n>=2)
ms_stocks <- filter(ms_stocks, fishery %in% nstocks$fishery)

# Build stock data
colnames(timeseries_values_views)
ms_data <- timeseries_values_views %>% 
  # Reduce to stocks of interest
  filter(stockid %in% ms_stocks$stockid) %>% 
  # Select columns of interest
  select(stockid, year, TCbest, TBbest, F, ERbest, BdivBmsypref, UdivUmgtpref) %>% 
  # Rename columns
  rename(catch=TCbest, biomass=TBbest, f=F, er=ERbest, bbmsy=BdivBmsypref, ffmsy=UdivUmgtpref) %>% 
  # Add fishery name
  left_join(select(ms_stocks, stockid, fishery), by="stockid") %>% 
  select(fishery, everything())


# Get and add life history info
################################################################################

# Species
spp <- sort(unique(ms_stocks$species))

# Get FishBase resilience values
res_fb <- datalimited2::resilience(spp)
fam <- select(freeR::taxa(spp), sciname, family)

# Add life history to stock data frame
ms_stocks1 <- ms_stocks %>% 
  left_join(fam, by=c("species"="sciname")) %>% 
  left_join(res_fb, by="species") %>% 
  mutate(resilience=as.character(resilience))

# Add missing resilience
ms_stocks1$resilience[ms_stocks1$stockid=="ILLEXNWATLC"] <- "High"
freeR::complete(ms_stocks1)

# Export data
################################################################################

# Export
data <- ms_data
stocks <- ms_stocks1
save(data, stocks, file=file.path(datadir, "ram_multispecies_fisheries.Rdata"))
