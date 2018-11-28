

# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Turn off sci notation
options(scipen=999)

# Packages
library(plyr)
library(dplyr)
library(freeR)

# Directories
datadir <- "data/simulated/data"
priordir <- "data/priors/data"
tabledir <- "data/simulated/tables"
codedir <- "data/priors/fishlife2"

# Read data
ne_spp <- read.csv(file.path(datadir, "usa_ne_multispecies_complex.csv"), as.is=T)
freeR::check_names(ne_spp$species)

# Read r prior function
source(file.path(codedir, "fishlife_r_priors.R"))

# Read r priors
rvals <- read.csv(file.path(priordir, "neubauer_etal_2013_r_values_clean.csv"), as.is=T)

# Read RAM Legacy Database
load("/Users/cfree/Dropbox/Prelim Database Files/Versions/RAM v4.41 (8-20-18)/DB Files With Assessment Data/DBdata (assessment data only).RData")


# Northeast multi-species fishery
################################################################################

# Stock info
ne_stocks <- stock %>%
  # Select columns
  select(stocklong, stockid, scientificname, commonname, region, areaid) %>% 
  # Add area info
  left_join(select(area, areaid, areaname), by="areaid") %>% 
  # Format columns
  rename(species=scientificname, comm_name=commonname, area=areaname) %>% 
  mutate(comm_name=freeR::sentcase(comm_name)) %>% 
  # Rearrange columns
  select(-areaid) %>% 
  # Filter to Northeast multispecies complex
  filter(region=="US East Coast") %>% # check species names here: check_names(ne_stocks$species)
  filter(species %in% ne_spp$species)

# Stock stats
ne_stats <- timeseries_values_views %>% 
  # Filter stocks of interest
  filter(stockid %in% ne_stocks$stockid) %>% 
  # Select columns of interest
  select(stockid, year, TB, SSB, TC, TL) %>% 
  # Calculate stats
  group_by(stockid) %>% 
  summarize(tc=quantile(TC, probs=0.90, na.rm=T),
            tl=quantile(TL, probs=0.90, na.rm=T),
            tb=max(TB, na.rm=T),
            ssb=max(SSB, na.rm=T)) %>% 
  # Add units
  left_join(select(timeseries_units_views, stockid, TB, SSB, TC, TL), by="stockid") %>% 
  rename(tb_units=TB, ssb_units=SSB, tc_units=TC, tl_units=TL) %>% 
  # Determine final stats
  mutate(c=ifelse(!is.na(tc), tc, tl),
         b=ifelse(!is.infinite(tb), tb, ssb),
         c_units=ifelse(!is.na(tc), tc_units, tl_units),
         b_units=ifelse(!is.infinite(tb), tb_units, ssb_units)) %>% 
  # Final cols
  select(stockid, c, b, c_units, b_units)

# Get FishLife life history
lh <- freeR::fishlife(ne_spp$species)
taxon <- freeR::taxa(ne_spp$species)
res <- datalimited2::resilience(ne_spp$species)

# Get r priors
set.seed(1)
r_vals <- data.frame(species=ne_spp$species, ln_r_mu=NA, ln_r_sd=NA, r_median=NA, r=NA, stringsAsFactors=F)
for(i in 1:nrow(r_vals)){
  r_info <- r_prior(r_vals$species[i])
  r_vals$ln_r_mu[i] <- r_info$ln_r_mu
  r_vals$ln_r_sd[i] <- r_info$ln_r_sd
  r_vals$ln_r_mu[i] <- r_info$ln_r_mu
  r_vals$r_median[i] <- r_info$r_median
  r_vals$r[i] <- rlnorm(1, meanlog=r_info$ln_r_mu, sdlog= r_info$ln_r_sd)
}

# Add stats to meta-data
ne <- ne_stocks %>% 
  left_join(ne_stats, by="stockid") %>% 
  # One stock of each species (one with greatest catch) 
  group_by(species) %>% 
  filter(!area%in%c("Northwestern Atlantic Coast", "Southern New England /Mid Atlantic")) %>% 
  filter(c==max(c)) %>%
  arrange(desc(c)) %>% 
  # Add life history info
  left_join(select(taxon, sciname, family), by=c("species"="sciname")) %>% 
  left_join(res, by="species") %>% 
  left_join(select(lh, species, m), by="species") %>% 
  left_join(select(r_vals, species, r), by="species") %>% 
  # Final K value
  # K is max(biomass) or 10*catch
  mutate(k=ifelse(b_units=="MT" & !is.na(b_units), b, c*10))

# Format for export
ne_final <- ne %>% 
  ungroup() %>% 
  mutate(fishery="Demersal trawl fishery",
         name=paste0(comm_name, " (", species, ")")) %>% 
  select(fishery, family, resilience, name, comm_name, species, r, k) %>% 
  mutate(k=ceiling(k/1000)) %>% 
  slice(1:10)

# Export
write.csv(ne_final, file.path(datadir, "demersal_trawl_fishery_10spp.csv"), row.names=F)

# Tuna RFMO
################################################################################

# Which RFMO has 10 stocks?
# ICCAT has most (International Commission for the Conservation of Atlantic Tunas)
rfmo_n <- assessment %>% 
  filter(mostrecent==999) %>% 
  left_join(select(stock, stockid, scientificname), by="stockid") %>% 
  group_by(assessorid) %>% 
  summarize(nstocks=n(), 
            nspp=n_distinct(scientificname)) %>% 
  filter(assessorid %in% c("IATTC", "IOTC", "ICCAT", "WCPFC", "SPC", "ISC", "CCSBT"))
  
# Stock info
t_stocks <- stock %>%
  # Select columns
  select(stocklong, stockid, scientificname, commonname, region, areaid) %>% 
  # Add area info
  left_join(select(area, areaid, areaname), by="areaid") %>% 
  # Add assessment info
  left_join(select(filter(assessment, mostrecent==999), stockid, assessorid), by="stockid") %>% 
  # Format columns
  rename(species=scientificname, comm_name=commonname, area=areaname, assessor=assessorid) %>% 
  mutate(comm_name=freeR::sentcase(comm_name),
         species=revalue(species, c("Tetrapturus albidus"="Kajikia albida"))) %>% 
  # Rearrange columns
  select(-areaid) %>% 
  # Filter to ICCAT
  filter(assessor=="ICCAT" & region=="Atlantic Ocean")

# Stock stats
t_stats <- timeseries_values_views %>% 
  # Filter stocks of interest
  filter(stockid %in% t_stocks$stockid) %>% 
  # Select columns of interest
  select(stockid, year, TB, SSB, TC, TL) %>% 
  # Calculate stats
  group_by(stockid) %>% 
  summarize(tc=quantile(TC, probs=0.90, na.rm=T),
            tl=quantile(TL, probs=0.90, na.rm=T),
            tb=max(TB, na.rm=T),
            ssb=max(SSB, na.rm=T)) %>% 
  # Add units
  left_join(select(timeseries_units_views, stockid, TB, SSB, TC, TL), by="stockid") %>% 
  rename(tb_units=TB, ssb_units=SSB, tc_units=TC, tl_units=TL) %>% 
  # Determine final stats
  mutate(c=ifelse(!is.na(tc), tc, tl),
         b=ifelse(!is.infinite(tb), tb, ssb),
         c_units=ifelse(!is.na(tc), tc_units, tl_units),
         b_units=ifelse(!is.infinite(tb), tb_units, ssb_units)) %>% 
  # Final cols
  select(stockid, c, b, c_units, b_units)

# Get FishLife life history
spp <- sort(unique(t_stocks$species))
lh <- freeR::fishlife(spp)
taxon <- freeR::taxa(spp)
res <- datalimited2::resilience(spp)

# Add stats to meta-data
t1 <- t_stocks %>% 
  left_join(t_stats, by="stockid") %>% 
  # One stock of each species (one with greatest catch) 
  group_by(species) %>%
  filter(c==max(c)) %>%
  arrange(desc(c)) %>%
  # Add life history info
  left_join(select(taxon, sciname, family), by=c("species"="sciname")) %>% 
  left_join(res, by="species") %>% 
  left_join(select(lh, species, m), by="species") %>% 
  # Add Neubauer r values
  left_join(unique(select(rvals, species, r_spp)), by="species") %>% 
  left_join(unique(select(rvals, family, r_fam)), by="family") %>% 
  mutate(r_nmort=0.87*m*2) %>% 
  # Final r/k values
  # K is max(biomass) or 10*catch
  mutate(k=ifelse(b_units=="MT" & !is.na(b_units), b, c*10),
         r=ifelse(!is.na(r_spp), r_spp, ifelse(!is.na(r_fam), r_fam, r_nmort)))

# Format for export
t_final <- t1 %>% 
  ungroup() %>% 
  mutate(fishery="Pelagic longline fishery") %>% 
  select(fishery, comm_name, species, r, k) %>% 
  mutate(k=ceiling(k/1000)) 




