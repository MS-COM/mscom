
# Function to plot MSCOM priors
plot_priors <- function(priors, true_values, species){
  
  # Unpack values
  r_priors <- priors[[1]]
  k_priors <- priors[[2]]
  d_priors <- priors[[3]]
  r_true <- true_values[[1]]
  k_true <- true_values[[2]]
  d_true <- true_values[[3]]
  
  # Plotting par
  nspp <- nrow(r_priors)
  par(mfcol=c(3,nspp))
  
  # Loop through parameters
  for(i in 1:nspp){
    
    # Plot r prior
    rvals <- rnorm(n=1000, mean=r_priors[i,1], sd=r_priors[i,2])
    hist(rvals, yaxt="n", xlab="r", ylab="", main=species[i], col="grey70", border=F)
    abline(v=r_true[i])
    
    # Plot K prior
    rvals <- rnorm(n=1000, mean=k_priors[i,1], sd=k_priors[i,2])
    hist(rvals, yaxt="n", xlab="K", ylab="", main=species[i], col="grey70", border=F)
    abline(v=k_true[i])
    
    # Plot depletion prior
    rvals <- rnorm(n=1000, mean=d_priors[i,1], sd=d_priors[i,2])
    hist(rvals, yaxt="n", xlab="Depletion", ylab="", main=species[i], col="grey70", border=F)
    abline(v=d_true[i])
    
  }
  

}