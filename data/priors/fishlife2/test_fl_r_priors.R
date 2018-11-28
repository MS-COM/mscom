

# Clean
rm(list = ls())

# Setup
################################################################################

# Directories
datadir <- "data/priors/fishlife2"

# Read code
source(file.path(datadir, "fishlife_r_priors.R"))


# Test r prior function on Vienna's species
################################################################################

# Species
spp <- c("Anisotremus interruptus", "Balistes polylepis", "Lutjanus argentiventris", 
         "Mycteroperca jordani", "Mycteroperca prionura", "Scarus perrico") #, " Octopus bimaculatus")
freeR::check_names(spp)

# Calc and export r prior
out <- lapply(spp, r_prior)
out <- do.call("rbind", out)
# write.csv(out, file.path("~/Desktop", "midriff_fish_r_priors.csv"), row.names=F)

# Test r prior function on RAM stocks
################################################################################

# Packages
library(plyr)
library(dplyr)
library(freeR)
library(rfishbase)

# Read RAM Legacy Database
load("/Users/cfree/Dropbox/Prelim Database Files/Versions/RAM v4.41 (8-20-18)/DB Files With Assessment Data/DBdata (assessment data only).RData")

# RAM species
ram_species <- stock %>% 
  dplyr::select(scientificname, commonname) %>% 
  rename(species=scientificname, comm_name=commonname) %>% 
  # Fix species names
  mutate(comm_name=freeR::sentcase(comm_name),
         species=revalue(species, c("Chrysophrys auratus"="Pagrus auratus",
                                    "Clupea pallasii"="Clupea pallasii pallasii",
                                    "Epinephelus flavolimbatus"="Hyporthodus flavolimbatus",
                                    "Epinephelus niveatus"="Hyporthodus niveatus",
                                    "Etrumeus teres"="Etrumeus sadina",
                                    "Merluccius gayi"="Merluccius gayi gayi",
                                    "Mullus barbatus"="Mullus barbatus barbatus",
                                    "Neoplatycephalus richardsoni"="Platycephalus richardsoni",
                                    "Psetta maxima"="Scophthalmus maximus",
                                    "Strangomera bentincki"="Clupea bentincki",
                                    "Tetrapturus albidus"="Kajikia albida",
                                    "Sardinops melanostictus"="Sardinops sagax"))) %>% 
  # Eliminate species not in FB/SLB (Loligo spp) and species not identified to species-level (spp)
  filter(!grepl("spp", species) & !grepl("Loligo", species)) %>% 
  unique()

# Check scientific names
freeR::check_names(ram_species$species)

# Format RAM species
ram_taxa <- taxa(ram_species$species)
ram <- ram_species %>% 
  left_join(ram_taxa, by=c("species"="sciname")) %>% 
  filter(type=="fish") %>% 
  arrange(species) %>% 
  mutate(mu=NA, sd=NA)

# Loop through RAM species
pdf("data/priors/fishlife2/ram_r_priors.pdf", height=11, width=8.5)
par(mfrow=c(6,4))
for(i in 1:nrow(ram)){
  
  out <- r_prior(ram$species[i])
  ram$mu[i] <- out$ln_r_mu
  ram$sd[i] <- out$ln_r_sd
  
}

dev.off()

