
# Life history priors from FishLife
################################################################################

#' Specify life history priors using FishLife
#'
#' Uses FishLife 2.0 to specify life history priors for >32,000 fish species.
#'
#' @param species Species to set priors for
#' @param traits Traits to get priors for
#' @return A data frame with priors for the listed traits.
#' @details
#' @references Thorson, J.T., Munch, S.B., Cope, J.M., Gao, J. (2017) 
#' Predicting life history parameters for all fishes worldwide. 
#' \emph{Ecological Applications} 27(8): 2262–2276.
#' 
#' Thorson, J.T. (in press) Predicting recruitment density dependence and intrinsic growth rate 
#' for all fishes worldwide using a data-integrated life-history model. \emph{Fish & Fisheries}
#' @examples
#' # Create vector of TOA scores and estimate status
# fl_priors(sciname="Gadus morhua", traits=c("r", "h", "amax", "m"))
#' @export
r_prior <- function(sciname){
  
  # Match species to FishLife
  genus <- stringr::word(sciname, 1)
  nwords_in_spp <- length(strsplit(sciname, " ")[[1]])
  species <- stringr::word(sciname, start=2, end=nwords_in_spp)
  species <- ifelse(species=="spp", "predictive", species)
  taxa_match <- FishLife::Search_species(Genus=genus, Species = species, add_ancestors=TRUE)$match_taxonomy
  
  # Get predictions from FishLife (mean and covariance)
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
  
  # Return log-normal prior parameters
  out <- data.frame(species=sciname, 
                    ln_r_mu=ln_r_mu, ln_r_sd=ln_r_sd, 
                    r_median=r_median, r_lo=r_lo, r_hi=r_hi)
  return(out)
  
}