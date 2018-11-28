

# Setup
################################################################################

# Install packages
# devtools::install_github("james-thorson/FishLife")
# devtools::install_github("cfree14/freeR")
# devtools::install_github("cfree14/datalimited2")
# devtools::install_github("cttedwards/lhm")

# Clear workspace
rm(list = ls())

# Packages
library(plyr)
library(dplyr)
library(freeR)
library(FishLife)
library(lhm)
library(rfishbase)

# Directories
datadir <- "data/simulated/data"

# Read species key
spp <- read.csv(file.path(datadir, "simulated_species.csv"), as.is=T)


# Helper functions
################################################################################

# Calculate maturity-at-age vector
matvec <- function(tmat, tmax){
  
  # Age vector
  ages <- 0:ceiling(tmax)
  
  # Maturity-at-age
  # Logistic growth: y = L / (1 + exp(-K*(x-x0)))
  L <- 1
  x0 <- tmat
  K <- tmat*2
  mats <- L / (1 + exp(-K*(ages-x0)))
  
  # Plot maturity-at-age
  par(mfrow=c(1,1))
  curve(L / (1 + exp(-K*(x-x0))), from=0, to=ceiling(tmax), n=100, las=1, bty="n", 
        xlab="Age (yr)", ylab="p(mature)", ylim=c(0,1))
  x <- c(0,tmat,tmax)
  y <- L / (1 + exp(-K*(x-x0)))
  points(x=x, y=y, pch=16)
  text(x=x[2:3], y=y[2:3], labels=c("Age-at-maturity", "Max age"), pos=c(4,1))
  
  # Return
  return(mats)
  
}



# Build data
################################################################################

# Check sci names
freeR::check_names(spp$species)

# Get taxanomic info
taxon <- freeR::taxa(spp$species)

# Get resilience
res <- datalimited2::resilience(spp$species)

# Get life history data
lh <- freeR::fishlife(spp$species)

# Get length-weight data
# Types: TL=total, FL=fork, SL=standard, OT=other, UNK=unknown
lw_orig <- length_weight(species_list=spp$species)
lw <- lw_orig %>% 
  # Select and rename columns
  select(sciname, Type, a, b) %>% 
  rename(species=sciname, lw_type=Type) %>% 
  # Replace NA type with "unknown" type
  mutate(lw_type=ifelse(is.na(lw_type), "UNK", lw_type)) %>% 
  # Calculate median LW parameters by species and type
  group_by(species, lw_type) %>% 
  summarize(n=n(),
            lw_a=median(a),
            lw_b=median(b)) %>% 
  # Use the type with the most data
  # If there is a tie, use the preferred type: TL > FL > SL > OT > UNK
  filter(n==max(n)) %>% 
  mutate(rank=as.numeric(factor(lw_type, levels=c("TL", "FL", "SL", "OT", "UNK")))) %>% 
  filter(rank==min(rank)) %>%
  # Rearrange
  arrange(species) %>% 
  select(species, lw_type, lw_a, lw_b)

# Build dataset
data <- spp %>% 
  left_join(select(taxon, sciname, family), by=c("species"="sciname")) %>% 
  left_join(res, by="species") %>% 
  left_join(lh, by="species") %>% 
  left_join(lw, by="species") %>% 
  mutate(r=0.87*m*2)


# Estimate intrinsic growth rate, r
for(i in 1:nrow(data)){
  
  # Subset data
  spp1 <- data$species[i]
  sdata <- filter(data, species==spp1)
  
  # Setup LHM object
  # ainf = assumed asymptotic age
  # iter = number of iterations
  rdat <- lhm(ainf = sdata$tmax_yr, iter = 200)
  
  # Fill LHM object
  nmort(rdat) <- list(mu = sdata$m)
  maturity(rdat) <- matvec(tmat=sdata$tmat_yr, tmax=sdata$tmax_yr)
  size(rdat) <- list(mu = list(Linf = sdata$linf_cm, k = sdata$k, t0 = 0))
  mass(rdat)     <- list(mu = list(a=sdata$lw_a, b = sdata$lw_b))
  sr(rdat) <- list(type = 'BH', mu = 0.90)
  
  # Estimate r
  r <- rCalc(rdat)
  
  # Plot r
  rmax <- freeR::ceiling1(max(r),0.5)
  hist(r, breaks=seq(0,rmax,0.05), las=1, col="grey30", main=spp1)
  
  
}











