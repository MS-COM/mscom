
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

# Load scenarios
load(file.path(datadir, "scenarios_pelagic_longline.Rdata"))

# Load simulatons
load(file.path(datadir, "sim1_pelagic_longline_byScenario.Rdata"))

# Which scenario?
s <- 4
sdata <- byScen[[s]]

# Format data for MS-cMSY
catch <- dcast(sdata, Year ~ Species, value.var="Catch")
yrs <- catch$Year
catch <- as.matrix(select(catch, -Year))
species <- colnames(catch)
res <- c("Medium", "Low", "Medium")

# True values
true <- NULL
true$r <- arrange(unique(sdata[,c("Species", "r")]), Species)$r
true$k <- arrange(unique(sdata[,c("Species", "K")]), Species)$K
true$b_ts <- dcast(sdata, Year ~ Species, value.var="Biomass")
true$er_ts <- dcast(sdata, Year ~ Species, value.var="ExploitRate")
true$bbmsy_ts <- dcast(sdata, Year ~ Species, value.var="BBmsy")

# Get true status for testing
bbmsy_end_true <- as.numeric(true$bbmsy_ts[nrow(true$bbmsy_ts),2:ncol(true$bbmsy_ts)])
status_end_true <- as.character(cut(bbmsy_end_true, breaks=c(0,0.5,1.5,999), labels=c("over", "fully", "under")))

# Fit MS-cMSY
out <- fit_mssra(catch=catch, years=yrs, stocks=species, res=res, id_fixed=T, npairs=20000)

# Plot MS-cMSY
plot_mssra(out, true)


