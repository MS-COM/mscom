
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Turn off sci notation
options(scipen=999)

# Packages
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# Directories
datadir <- "data/simulated"
codedir <- "code/helper_functions"

# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Format data
################################################################################

# Read data
load(file.path(datadir, "sim_pelagic_longline1.Rdata"))

# Format data
data$bbmsy <- data$Biomass / data$Bmsy
data$uumsy <- data$ExploitRate / data$Umsy

# Format catch
catch <- dcast(data, Year ~ Species, value.var="Catch")
catch <- as.matrix(select(catch, -Year))

# Species info
# billfish, shark, tuna
species <- colnames(catch)
r_true <- c(0.67, 0.11, 0.60)
k_true <- c(227, 738, 304)
res <- c("Medium", "Low", "Medium")
yrs <- 1:nrow(catch)

# True
true <- NULL
true$r <- r_true
true$k <- k_true
true$b_ts <- dcast(data, Year ~ Species, value.var="Biomass")
true$er_ts <- dcast(data, Year ~ Species, value.var="ExploitRate")
true$bbmsy_ts <- dcast(data, Year ~ Species, value.var="bbmsy")


# Run functions
################################################################################

# Run and plot
out <- fit_mssra(catch=catch, years=yrs, stocks=species, res=res, id_fixed=F)
plot_mssra(out, true)



