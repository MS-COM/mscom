
# Clear workspace
# rm(list = ls())

# Setup
################################################################################

# Packages
library(TMB)

# FishLife
# devtools::install_github("james-thorson/FishLife")
# devtools::install_github("james-thorson/FishLife", ref="add-recruitment")
library(FishLife)

# Directories
# datadir <- "data/priors/fishlife2"

# FishLife 1.0 database
# Return0 <- database

# FishLife 2.0 database
# load(file.path(datadir, "Return.RData"))
load("data/priors/fishlife2/Return.RData")
Return_orig <- Return
rm(Return)
Return <- NULL
Return$beta_gv <- Return_orig$beta_gv
Return$Cov_gvv <- Return_orig$Cov_gvv
rm(Return_orig)

# Learn about database structure
################################################################################

# Inspect
# names(Return0$Y_ij)
# names(Return$Y_ij)
# colnames(Return$beta_gv)

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


# Function to calculate r prior
################################################################################

# For testing: r_prior(sciname)
# sciname <- "Allocyttus niger"

# R prior function
r_prior <- function(sciname, plot=F){
  
  # Match species to FishLife
  genus <- stringr::word(sciname, 1)
  nwords_in_spp <- length(strsplit(sciname, " ")[[1]])
  species <- stringr::word(sciname, start=2, end=nwords_in_spp)
  species <- ifelse(species=="spp", "predictive", species)
  suppressMessages( taxa_match <- FishLife::Search_species(Genus=genus, Species = species, add_ancestors=TRUE)$match_taxonomy)
  
  # Get predictions from FishLife (mean and covariance)
  # Use all predicted values (Loo to ln_MASPS, #1-11) and one derived value (intrinsic growth rate, #17)
  params <- colnames(Return$beta_gv)
  mus <- Return$beta_gv[rownames(Return$beta_gv)==taxa_match[[1]], c(1:11,17)]
  covs <- Return$Cov_gvv[rownames(Return$Cov_gvv)==taxa_match[[1]], c(1:11,17), c(1:11,17)]
  
  # R prior
  ln_r_mu <- mus["ln_r"]
  ln_r_sd <- sqrt(covs["ln_r", "ln_r"])
  r_mode <- exp(ln_r_mu - ln_r_sd^2)
  r_mean <- exp(ln_r_mu + ln_r_sd^2/2)
  r_median <- exp(ln_r_mu)
  r_lo <- qlnorm(0.025, meanlog=ln_r_mu, sdlog=ln_r_sd)
  r_hi <- qlnorm(0.975, meanlog=ln_r_mu, sdlog=ln_r_sd)
  
  # Only do the following if plotting
  if(plot==T){
  
    # Use mean and covariance to simulate R values
    # This frequently fails because the COV matrix isn't positive definite
    # I'm really only doing this to confirm that the univariate r prior lines up nicely with the multivariate simulations
    ndraws <- 10000
    draws <- try(data.frame(MASS::mvrnorm(n=ndraws, mu=mus, Sigma=covs)))
    
    # If no problems, plot similations and the prior
    if(!inherits(draws, "try-error")){
      
      # Plot simulations
      ln_rs <- draws$ln_r
      rs <- exp(ln_rs)
      xmin <- floor(min(rs)/0.2)*0.2
      xmax <- ceiling(max(rs)/0.2)*0.2
      hist(rs, col="grey80", border=F, 
           xlim=c(xmin, xmax), yaxt="n", xlab="Intrinsic growth rate, r", ylab="", main=sciname)
  
      # Add log-normal prior
      x <- seq(xmin, xmax, length.out=1000)
      hx <- dlnorm(x, meanlog=ln_r_mu, sdlog=ln_r_sd)
      par(new=T)
      plot(x=x, y=hx, type="l", xaxt="n", yaxt="n", frame=F, xlab="", ylab="")
      
      # Add lines
      abline(v=r_median, lty=2)
      abline(v=r_lo, lty=3)
      abline(v=r_hi, lty=3)
      
    # If there is a problem with the simulations, just plot the prior
    }else{
      
      # Plot log-normal prior
      xmin <- floor(qlnorm(0.02, meanlog=ln_r_mu, sdlog=ln_r_sd)/0.1)*0.1
      xmax <- ceiling(qlnorm(0.99, meanlog=ln_r_mu, sdlog=ln_r_sd)/0.1)*0.1
      x <- seq(xmin, xmax, length.out=1000)
      hx <- dlnorm(x, meanlog=ln_r_mu, sdlog=ln_r_sd)
      plot(x=x, y=hx, type="l", yaxt="n", frame=F, xlab="Intrinsic growth rate, r", ylab="", main=sciname)
      
    }
  
  }
  
  # Return log-normal prior parameters
  out <- data.frame(species=sciname, 
                    ln_r_mu=ln_r_mu, ln_r_sd=ln_r_sd, 
                    r_median=r_median, r_lo=r_lo, r_hi=r_hi)
  return(out)
  
}

