
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(TMB)

# Directories
tmbdir <- "code/tmb"
datadir <- "data/simulated"
codedir <- "code/helper_functions"

## for windows
maindir <- "C:\\merrill\\MS-COM\\mscom"
tmbdir <- file.path(maindir, "code","tmb")
datadir <- file.path(maindir, "data","simulated")
codedir <- file.path(maindir, "code","helper_functions")

# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))

###############################
## Scenarios
###############################

## variability
var_vec <- c("deterministic")#, "variable")

## initial depletion
depl_vec <- c("unfished")#, "fished_same", "fished_diff")

## effort dynamics
effdyn_vec <- c("one-way")#, "two-way")

## magnitude of maximum fishing mortality
fmag_vec <- c("past_fmsy")#, "at_fmsy", 'below_fmsy')

## length of time series
time_vec <- 50

## strength of effort correlation
fdev_vec <- 0.1

## scenario combinations
scenarios <- expand.grid("PopDyn"=var_vec,
						"InitDepl"=depl_vec,
						"EffDyn"=effdyn_vec,
						"MaxF"=fmag_vec,
						"TimeSeries"=time_vec,
						"EffSD"=fdev_vec)

# Simulate data
################################################################################

# Fishery #1: Pelagic longline
# Tuna: Yellowfin tun (Thunnus albacares)
# Billfish = Striped marlin (Kajikia audax)
# Shark = Mako shark (Isurus oxyrinchus)

sims <- lapply(1:nrow(scenarios), function(x){
	# Parameters
	species <- c("tuna","billfish","shark")

	## initial depletion options
	if(scenarios[x,'InitDepl']=="unfished") init_depl <- rep(1,length(species))
	if(scenarios[x,'InitDepl']=="fished_same") init_depl <- rep(0.8,length(species))
	if(scenarios[x,'InitDepl']=="fished_diff") init_depl <- seq(from=0.4, to=0.8, length.out=length(species))

	## maximum F options
	if(scenarios[x,"MaxF"]=="past_fmsy") maxF <- rep(1.1,length(species))
	if(scenarios[x,"MaxF"]=="at_fmsy") maxF <- rep(1, length(species))
	if(scenarios[x,"MaxF"]=="below_fmsy") maxF <- rep(0.8, length(species))

	## recruitment variability
	if(scenarios[x,"PopDyn"]=="deterministic"){
		sigR <- 0
		rho <- 0
	}
	if(scenarios[x,"PopDyn"]=="variable"){
		sigR <- 0.7
		rho <- 0.4
	}

	input_df <- data.frame("SpeciesName"=species,
							"EffDyn"=scenarios[x,"EffDyn"],
							"InitDepl"=init_depl,
							"MSYscalar"=maxF,
							"SigmaR"=sigR, "rho"=rho,
							"SigmaF"=scenarios[x,"EffSD"],
							"r" = c(0.6, 0.67, 0.11),
							"K"=c(304, 227, 738),
							"p"=rep(0.2,length(species)))

	## simulate populations
	sim <- sim_pops(input_df=input_df, 
					nyears=scenarios[x,"TimeSeries"], 
					seed=1122, 
					model="biomass-dynamic")

	## plot simulation results
	ggplot(sim %>% dplyr::filter(variable %in% c("RelativeCatch","ExploitRate","RelativeEffort","Depletion"))) +
	geom_line(aes(x=Year, y=value, colour=Species), lwd=2) +
	facet_wrap(~variable, scales='free_y') +
	theme_lsd() +
	coord_cartesian(ylim=c(0,1.01))

	## plot simulation results
	ggplot(sim %>% dplyr::filter(variable %in% c("BBmsy","UUmsy"))) +
	geom_line(aes(x=Year, y=value, colour=Species), lwd=2) +
	facet_wrap(~variable, scales='free_y') +
	theme_lsd()

})


# Extract catch
catch_oneway <- sim_oneway %>% filter(variable=="Catch")
input_oneway <- sapply(1:length(species), function(x){
	sub <- catch_oneway %>% filter(Species==species[x])
	df <- sub$value
	return(df)
})
colnames(input_oneway) <- species

# Format data
true <- oneway_mat
data <- dcast(sim_oneway, Species + Year ~ variable, value.var="value")


# Export data
################################################################################

# Export data
save(data, true, file=file.path(datadir, "sim_pelagic_longline.Rdata"))







