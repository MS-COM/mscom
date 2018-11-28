
# Load FishLife 2.0 database
# devtools::install_github("james-thorson/FishLife", ref="add-recruitment", force=T)
library(FishLife)
load("~/Dropbox/Chris/Rutgers/projects/mscom/data/priors/fishlife2/Return.RData")

# Predicted variables
# -------------------------------------------
# Loo - asymptotic length (Linf, cm)
# K - growth coefficient (K)
# Winfinity - Asymptotic mass (Winf, g)
# tmax - maximum age (Amax, yr)
# tm - age at maturity (Amat, yr)
# M - mortality rate (M, 1/yr)
# Lm - length at maturity (Lmat, cm)
# Temperature - average temperature (T, °C)
# ln_var - marginal standard deviation of recruitment variability (τ)
# rho - autocorrelation of recruitment variability (ρ)
# ln_MASPS - maximum annual spawners per spawner (r)

# Derived variables
# -------------------------------------------
# ln_margsd - standard deviation for recruitment (σ): σ = sqrt(τ^2 / (1-ρ^2))
# h / logitbound_h - steepness (h): h = ρ / (4 + ρ)
# ln_Fmsy - FMSY
# ln_Fmsy_over_m - FMSY/M ratio
# r / ln_r - Intrinsic growth rate (r): dominant eigen value for Leslie matrix w/ assumptions: length-weight b=3.04, VonB t0=-0.1, maturity ogive slope=0.25*tmat
# G / ln_G - Generation time (G, yr)


# Function to get mean estimates
fishlife2 <- function(species){
  
  # Setup container
  fl <- data.frame(species=sort(unique(species)),
                   linf_cm=NA, k=NA, winf_g=NA, tmax_yr=NA, tmat_yr=NA,
                   m=NA, lmat_cm=NA, temp_c=NA, 
                   rho=NA, h=NA, r=NA, fmsy=NA, g_yr=NA, stringsAsFactors=F)
  
  # Loop through species
  for(i in 1:nrow(fl)){
    
    # Match species to FishLife
    sciname <- fl$species[i]
    genus <- stringr::word(sciname, 1)
    nwords_in_spp <- length(strsplit(sciname, " ")[[1]])
    spp <- stringr::word(sciname, start=2, end=nwords_in_spp)
    spp <- ifelse(spp=="spp", "predictive", spp)
    try(taxa_match <- FishLife::Search_species(Genus=genus, Species = spp, add_ancestors=TRUE)$match_taxonomy)
    
    # Get predictions from FishLife (mean and covariance)
    if(inherits(taxa_match, "try-error")){
      # Record blanks
      fl[i,2:ncol(fl)] <- rep(NA, ncol(fl)-1)
    }else{
      # Values are in log-scale except temperature
      params <- colnames(Return$beta_gv)
      mus <- Return$beta_gv[rownames(Return$beta_gv)==taxa_match[[1]], ]
      mus_use <- mus[c("Loo", "K", "Winfinity", "tmax", "tm", "M", "Lm", "Temperature", "rho", "h", "r", "ln_Fmsy", "G")]
      fl[i,2:ncol(fl)] <- mus_use
    }
    
  }
  
  # Exponentiate columns
  # These columns are not log-transformed: "temp_c", "rho", "h", "r", "g_yr", "fmsy"
  log_cols <- c("linf_cm", "k", "winf_g", "tmax_yr", "tmat_yr", "m", "lmat_cm")
  fl[,log_cols] <- exp(fl[,log_cols])

  # Return
  return(fl)
  
}





