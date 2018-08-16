
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# Directories
tmbdir <- "code/tmb"
datadir <- "data/simulated"
codedir <- "code/helper_functions"

# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Simulate data
################################################################################

# Fishery #1: Pelagic longline
# Tuna: Yellowfin tun (Thunnus albacares)
# Billfish = Striped marlin (Kajikia audax)
# Shark = Mako shark (Isurus oxyrinchus)

# Parameters
species <- c("tuna", "billfish", "shark")
oneway_mat <- data.frame("SpeciesName"=species,
						"Fdynamics"="One-way",
						"InitialDepl"=1,
						"PercentFcrash"=c(0.4, 0.6, 0.8),
						"SigmaR"=0, "rho"=0,
						"SigmaF"=0,
						"r" = c(0.6, 0.67, 0.11),
						"K"=c(304, 227, 738))

# Run simulation
sim_oneway <- sim_pops(input_mat=oneway_mat,
				nyears=50,
				seed=123, 
				model="biomass-dynamic")

# Plot simulation results
p_oneway <- ggplot(sim_oneway %>% dplyr::filter(variable %in% c("RelativeCatch","ExploitRate","RelativeEffort","Depletion"))) +
	geom_line(aes(x=Year, y=value, colour=Species), lwd=2) +
	facet_wrap(~variable, scales='free_y') +
	theme_lsd() +
	coord_cartesian(ylim=c(0,1.01))
p_oneway

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







