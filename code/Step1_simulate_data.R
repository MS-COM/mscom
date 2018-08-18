
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
var_vec <- c("deterministic", "variable")

## initial depletion
depl_vec <- c("unfished", "fished")

## effort dynamics
effdyn_vec <- c("one-way", "two-way")

## magnitude of maximum fishing mortality
fmag_vec <- c("past_fmsy")#, "at_fmsy", 'below_fmsy')

## length of time series
nyears_vec <- 50

## strength of effort correlation
fdev_vec <- c(0,0.1)

## scenario combinations
scenarios <- expand.grid("PopDyn"=var_vec,
						"DataStart"=depl_vec,
						"EffDyn"=effdyn_vec,
						"MaxF"=fmag_vec,
						"Nyears"=nyears_vec,
						"EffSD"=fdev_vec)

# Simulate data
################################################################################

# Fishery #1: Pelagic longline
# Tuna: Yellowfin tun (Thunnus albacares)
# Billfish = Striped marlin (Kajikia audax)
# Shark = Mako shark (Isurus oxyrinchus)
itervec <- 1:10

## run multiple iterations
sim_byIter <- lapply(1:length(itervec), function(i){
  sim <- sim_scenarios(savedir=datadir, scen_df=scenarios, fishery="pelagic_longline", seed=itervec[i], iter=itervec[i])
  return(sim)
})

save(scenarios, file=file.path(datadir, "scenarios_pelagic_longline.Rdata"))

## plot example - first iteration, first scenario
	simdf <- melt(sim_byIter[[1]][[1]], id.var=c('Species','Year'))
	sim <- sim_byIter[[1]][[1]]

	# ## plot simulation results
	ggplot(simdf %>% dplyr::filter(variable %in% c("RelativeCatch","ExploitRate","RelativeEffort","Depletion"))) +
	geom_line(aes(x=Year, y=value, colour=Species), lwd=2) +
	facet_wrap(~variable, scales='free_y') +
	theme_lsd() +
	coord_cartesian(ylim=c(0,1.01))

	## plot simulation results
	ggplot(sim) +
	geom_line(aes(x=BBmsy, y=UUmsy, colour=Species), lwd=2) +
	theme_lsd() +
	geom_hline(yintercept=1) +
	geom_vline(xintercept=1)



