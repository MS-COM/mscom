

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(plyr)
library(dplyr)
library(reshape2)

# Directories
outdir <- "code/ram_stocks/data"
codedir <- "code/ram_stocks"

# Read RAM multispecies fisheries
ms_key <- read.csv(file.path(outdir, "RAM_MS_fisheries.csv"), as.is=T)

# Read RAM Legacy Database
load("/Users/cfree/Dropbox/Prelim Database Files/Versions/RAM v4.41 (8-20-18)/DB Files With Assessment Data/DBdata (assessment data only).Rdata")

# Read helper functions
source(file.path(codedir, "helper_functions.R"))



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
         comm_name=revalue(comm_name, c("albacore tuna"="Albacore tuna",
                                        "American Plaice"="American plaice",
                                        "bigeye tuna"="Bigeye tuna",
                                        "blue grenadier"="Blue grenadier",
                                        "blue marlin"="Blue marlin",
                                        "pacific bluefin tuna"="Pacific bluefin tuna",
                                        "Silverfish"="Silver warehou",
                                        "skipjack tuna"="Skipjack tuna",
                                        "striped marlin"="Striped marlin",
                                        "swordfish"="Swordfish",
                                        "white marlin"="White marlin",
                                        "Winter Flounder"="Winter flounder",
                                        "Witch Flounder"="Witch flounder"))) %>% 
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

# Inspect data
sort(unique(ms_stocks$comm_name))
check_species(ms_stocks$species)

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
  filter(stockid %in% ms_stockids) %>% 
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

# Gather FishLife life history
lh_fl <- fishlife(spp)

# Get FishBase resilience values
res_fb <- resilience(spp)

# Add life history to stock data frame
ms_stocks <- ms_stocks %>% 
  left_join(res_fb, by="species") %>% 
  left_join(lh_fl, by="species")


# Package data for analysis
################################################################################

# Setup list
data <- list()

# Loop through fisheries, package data, and put into master list
fisheries <- sort(unique(ms_stocks$fishery))
for(i in 1:length(fisheries)){
  
  # Setup fishery list
  flist <- NULL
  flist$info <- arrange(filter(ms_stocks, fishery==fisheries[i]), stockid)
  
  # Subset data, reshape, and store
  fdata <- filter(ms_data, fishery==fisheries[i])
  flist$catch <- dcast(fdata, year ~ stockid, value.var="catch")
  flist$catch_use <- na.omit(dcast(fdata, year ~ stockid, value.var="catch"))
  flist$biomass <- dcast(fdata, year ~ stockid, value.var="biomass")
  flist$er <- dcast(fdata, year ~ stockid, value.var="er")
  flist$f <- dcast(fdata, year ~ stockid, value.var="f")
  flist$bbmsy <- dcast(fdata, year ~ stockid, value.var="bbmsy")
  flist$ffmsy <- dcast(fdata, year ~ stockid, value.var="ffmsy")
  
  # Add fishery list to master list
  data[[i]] <- flist
  
}

# Add names to list
names(data) <- fisheries


# Export data
################################################################################

# Export
save(data, ms_stocks, ms_data, file=file.path(outdir, "ram_multispecies_fisheries.Rdata"))
