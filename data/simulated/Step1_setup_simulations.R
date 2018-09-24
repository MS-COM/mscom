

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
# Types: TL=total, FL=fork, SL=standard, OT=other
lw_orig <- length_weight(species_list=spp$species)
lw <- lw_orig %>% 
  select(sciname, Type, a, b) %>% 
  rename(species=sciname, type=Type) %>% 
  group_by(species, type) %>% 
  summarize(n=n(),
            a=median(a),
            b=median(b)) %>% 
  filter(n==max(n))

# Build dataset
data <- spp %>% 
  left_join(select(taxon, sciname, family), by=c("species"="sciname")) %>% 
  left_join(res, by="species") %>% 
  left_join(lh, by="species")


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
  # maturity(rdat) <- matvec(tmat=sdata$tmat_yr, tmax=sdata$tmax_yr)
  size(rdat) <- list(mu = list(Linf = sdata$linf_cm, k = sdata$k, t0 = 0))
  mass(rdat)     <- list(mu = list(a = 1.88e-9, b = 3.305))
  sr(rdat) <- list(type = 'BH', mu = 0.90, cv = 0.10)
  
  # Estimate r
  r <- rCalc(rdat)
  
  # Plot r
  plot(r)
  
  
}


# initialise lhm data object for calculation of r with uncertainty
rdat <- lhm(ainf = 100, iter = 200)

# then life-history vectors can be assigned to each iteration
# with or without uncertainty
nmort(rdat)    <- list(mu = 0.18)
maturity(rdat) <- c(0.0,0.01,0.02,0.06,0.14,0.28,0.50,0.72,0.86,0.94,0.98,0.99,1.00)
size(rdat)     <- list(mu = list(Linf = 106.5, k = 0.229, t0 = 0.01))
mass(rdat)     <- list(mu = list(a = 1.88e-9, b = 3.305))
sr(rdat)       <- list(type = 'BH', mu = 0.90, cv = 0.10)

# calculate r prior and fit log-normal distribution
r <- rCalc(rdat)
plot(r)



x <-c(0.0,0.01,0.02,0.06,0.14,0.28,0.50,0.72,0.86,0.94,0.98,0.99,1.00)
plot(1:length(x), x)


tmax <- sdata$tmax_yr
tmat <- sdata$tmat_yr







