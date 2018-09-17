
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

# Read data
spp <- import(file.path(datadir, "GulfThaiSpecies.xlsx"), which=1)
c_rw <- import(file.path(datadir, "GulfThaiSpecies.xlsx"), which=2)
c_fao <- import(file.path(datadir, "GulfThaiSpecies.xlsx"), which=3)
c_fao_t <- import(file.path(datadir, "GulfThaiSpecies.xlsx"), which=4)
er <- import(file.path(datadir, "Gulf_of_Thailand_ER.xlsx"))


# Format data
################################################################################

# Format species key
species <- spp %>% 
  rename(sci_name_orig=sciname, comm_name=commname, sci_name=NewSciName) %>% 
  select(comm_name, sci_name, sci_name_orig) %>% 
  mutate(sci_name=ifelse(is.na(sci_name), sci_name_orig, sci_name),
         sci_name=stringr::str_trim(sci_name), 
         sci_name=revalue(sci_name, c("Selar crumenopthalmus"="Selar crumenophthalmus")),
         comm_name=stringr::str_trim(comm_name)) %>% 
  rbind(c("Monocle breams", rep("Scolopsis spp.", 2)))
check_species(species$sci_name)

# Format FAO catch data
fao_long <- c_fao_t %>% 
  setNames(., paste0("col", 1:ncol(c_fao_t))) %>% 
  select(col10:col16) %>% 
  setNames(., .[1,]) %>% 
  slice(-1) %>% 
  rename(year="Year") %>% 
  melt(id.vars="year", variable.name="comm_name", value.name="catch") %>% 
  select(comm_name, year, catch) %>% 
  mutate(year=as.numeric(year),
         comm_name=as.character(comm_name))
unique(fao_long$comm_name)[!unique(fao_long$comm_name) %in% species$comm_name]

# Format ER data
er_long <- er %>% 
  rename(year=Year) %>% 
  melt(id.vars="year", variable.name="sci_name", value.name="er") %>% 
  select(sci_name, year, er) %>% 
  mutate(sci_name=as.character(sci_name),
         sci_name=revalue(sci_name, c("Atul mate"="Atule mate",
                                     "Loligo chinensis"="Uroteuthis chinensis",
                                     "Loligo duvauceli"="Uroteuthis duvaucelii",
                                     "Lutanus lineolatus"="Lutjanus lutjanus",
                                     "Sepia aculenta"="Sepia aculeata",
                                     "Scolopsis taeniopterus"="Scolopsis taenioptera")))
unique(er_long$sci_name)[!unique(er_long$sci_name) %in% species$sci_name]

# Format SAUP data
rw_long <- c_rw %>% 
  rename(year=Year, sci_name=TaxonName, comm_name=CommonName, fleet=FleetGearName, 
         c_rep=IndReported, c_iuu=IndIUU, c_dis=IndDiscards, 
         c_rep1=NIndReported, c_iuu1=NIndIUU, c_dis1=NIndDiscards) %>% 
  select(fleet, comm_name, sci_name, year, c_rep:c_dis1) %>% 
  arrange(fleet, comm_name, year)
sort(unique(rw_long$comm_name))


# Visualize data
################################################################################

fleet_spp <- rw_long %>% 
  group_by(fleet, comm_name) %>% 
  summarize(tc_avg=mean(c_rep)) %>% 
  arrange(fleet, desc(tc_avg))

ggplot(rw_long, aes(x = year, y = c_rep, colour = comm_name)) +
  labs(y="Catch", x="") +
  geom_area(aes(colour = comm_name, fill=comm_name), position = 'stack') +
  facet_wrap(~ fleet, ncol=4)

ggplot(er_long, aes(x = year, y = er)) +
  geom_line() +
  facet_wrap(~ sci_name, ncol=4)


# Build data
################################################################################

# Examine trawl fishery data from Watson data
trawl_spp <- rw_long %>% 
  filter(fleet=="Trawl" & year>2005) %>% 
  group_by(comm_name, sci_name) %>% 
  summarize(c_avg=mean(c_rep, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(c_perc=c_avg/sum(c_avg)*100,
         name_format=paste0(comm_name, " (", sci_name, ")"),
         c_format=paste0(round(c_avg,1), " / ", round(c_perc, 1), "%")) %>%
  select(name_format, comm_name, sci_name, c_avg, c_perc, c_format) %>% 
  arrange(desc(c_perc))

# Target species
spp_do <- c("Blue swimming crab", "Bigeye scad", "Largehead hairtail")

# Merge data
data <- fao_long %>% 
  filter(comm_name%in%spp_do & !is.na(catch)) %>% 
  left_join(select(species, comm_name, sci_name), by="comm_name") %>% 
  select(comm_name, sci_name, year, catch) %>% 
  left_join(er_long, by=c("year", "sci_name")) %>% 
  mutate(catch=as.numeric(catch))

# What percentage of these species are caught outside the trawl fishery?
spp_do_stats <- rw_long %>% 
  filter(comm_name %in% spp_do & year > 2005) %>% 
  group_by(comm_name, fleet) %>% 
  summarize(c_avg=mean(c_rep, na.rm=T))
  

# Export data
################################################################################

# Export data
write.csv(data, file.path(outdir, "gulf_of_thailand_trawl_catch_data.csv"), row.names=F)
write.csv(trawl_spp, file.path(tabledir, "TableX_got_summary.csv"), row.names=F)

