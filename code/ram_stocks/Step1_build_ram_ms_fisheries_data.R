

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

# Read RAM Legacy Database
load("/Users/cfree/Dropbox/Prelim Database Files/Versions/RAM v4.40 (6-4-18)/DB Files With Assessment Data/DBdata.RData")

# Read helper functions
source(file.path(codedir, "helper_functions.R"))

# Identify stocks
################################################################################

# Look for multi-species fisheries
stocks_sub <- filter(stock, region=="US East Coast")

# Australia / New Zealand
##################################

# NZ Chatham Rise middle-depth
nzcr_stocks <- c("HOKIENZ", "SOUTHHAKECR", "NZLINGLIN3-4")

# SE Australia Scalefish/Sharks
# Blue grenadier (Macruronus novaezelandiae), Tiger flathead (Neoplatycephalus richardsoni)
# Silver warehou (Seriolella punctata), Gummy shark (Mustelus antarcticus), Pink ling (Genypterus blacodes)
sessf_stocks <- c("BGRDRSE", "NZLINGESE", "TIGERFLATSE", "SILVERFISHSE")

# United States
##################################

# Mid-Atlantic mackerel, squid, butterfish
msb_stocks <- c("BUTTERGOMCHATT", "ILLEXNWATLC")

# West Coast DTS fishery
dts_stocks <- c("DSOLEPCOAST", "SSTHORNHPCOAST", "LSTHORNHPCOAST", "SABLEFPCOAST")

# BSAI groundfish 
bsai_stocks <- c("YSOLEBSAI", "PCODBSAI", "WPOLLAI", "WPOLLEBS")

# NE groundfish
# NAFO 5Z = Georges Bank
# NAFO 5Y = Gulf of Maine (GOM)
gb_gf_stocks <- c("CODGB", "HADGB", "WINFLOUN5Z", "YELLGB")
gom_gf_stocks <- c("CODGOM", "HAD5Y", "WINFLOUND5Y", "WITFLOUN5Y", "YELLCCODGOM")
gbgom_gf_stocks <- c("AMPL5YZ", "ACADREDGOMGB", "POLL5YZ", "WHAKEGBGOM")

# Merge all
##################################

# All multi-species stocks
ms_stockids <- c(nzcr_stocks, sessf_stocks,
                 dts_stocks, bsai_stocks,
                 gb_gf_stocks, gom_gf_stocks, gbgom_gf_stocks)


# Build data
################################################################################

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
  mutate(species=revalue(species, c("Neoplatycephalus richardsoni"="Platycephalus richardsoni")),
         comm_name=revalue(comm_name, c("American Plaice"="American plaice",
                                        "blue grenadier"="Blue grenadier",
                                        "Silverfish"="Silver warehou",
                                        "Winter Flounder"="Winter flounder",
                                        "Witch Flounder"="Witch flounder"))) %>% 
  # Add fishery name
  mutate(fishery=NA,
         fishery=ifelse(stockid %in% nzcr_stocks, "NZ Chatham Rise", fishery),
         fishery=ifelse(stockid %in% sessf_stocks, "Australia SESSF", fishery),
         fishery=ifelse(stockid %in% msb_stocks, "US Mid-Atlantic MSB", fishery),
         fishery=ifelse(stockid %in% dts_stocks, "US West Coast DTS", fishery),
         fishery=ifelse(stockid %in% bsai_stocks, "US BSAI groundfish", fishery),
         fishery=ifelse(stockid %in% gb_gf_stocks, "US GB groundfish", fishery),
         fishery=ifelse(stockid %in% gom_gf_stocks, "US GOM groundfish", fishery),
         fishery=ifelse(stockid %in% gbgom_gf_stocks, "US GB/GOM groundfish", fishery)) %>% 
  # Add catch/biomass units
  left_join(select(timeseries_ids_views, stockid, TCbest, TBbest), by="stockid") %>% 
  rename(tc_units=TCbest, tb_units=TBbest) %>% 
  # Add reference points
  left_join(select(bioparams_values_views, stockid, MSY, TBmsybest, Fmsy), by="stockid") %>% 
  rename(msy=MSY, bmsy=TBmsybest, fmsy=Fmsy) %>% 
  left_join(select(bioparams_ids_views, stockid, MSY, TBmsybest, Fmsy), by="stockid") %>% 
  rename(msy_units=MSY, bmsy_units=TBmsybest, fmsy_units=Fmsy) %>% 
  select(fishery, everything())

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
