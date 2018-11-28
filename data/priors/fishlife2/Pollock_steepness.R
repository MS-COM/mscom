


# Load FishLife 2.0 package
devtools::install_github("james-thorson/FishLife", ref="add-recruitment")
library(FishLife)

# Load FishLife 2.0 package data
load( "data/priors/fishlife2/Return.RData" )


# Look at estimates
params = matrix(c("K", "M", "Winfinity", "Loo", "tmax", "tm", "Lm", "Temperature", "tm", "rho", "Winfinity", "h"), ncol=2, byrow=TRUE)
Cov_pred = Plot_taxa( Search_species(Genus="Theragra",Species="chalcogramma",add=FALSE)$match_taxonomy, params=params,
  Cov_gjj=Return$Cov_gvv, Mean_gj=Return$beta_gv, ParentChild_gz=Return$ParentChild_gz, Y_ij=Return$Y_ij, mfrow=c(3,2) )

# Mean steepness
Cov_pred[[1]]$Mean_pred['h']
sqrt( Cov_pred[[1]]$Cov_pred['h','h'] )

