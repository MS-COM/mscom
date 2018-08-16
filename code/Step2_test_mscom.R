
# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(dplyr)
library(reshape2)
library(ggplot2)
library(TMB)
library(RColorBrewer)

# Directories
# tmbdir <- "code/tmb"
# datadir <- "data/simulated"
# ramdir <- "code/ram_stocks/data"
# codedir <- "code/helper_functions"
# plotdir <- "figures"

## for windows
maindir <- "C:\\merrill\\MS-COM\\mscom"
tmbdir <- file.path(maindir, "code","tmb")
datadir <- file.path(maindir, "data","simulated")
ramdir <- file.path(maindir, "code", "ram_stocks","data")
codedir <- file.path(maindir, "code", "helper_functions")
plotdir <- file.path(maindir,"figures")

setwd(tmbdir)
compile("mscom.cpp")


# Read helper functions
helpers <- list.files(codedir)
ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Format data
################################################################################

# Read data
load(file.path(datadir, "sim_pelagic_longline.Rdata"))

# Format data
data$bbmsy <- data$Biomass/data$Bmsy
data$uumsy <- data$ExploitRate/data$Umsy

# Format catch
catch <- dcast(data, Year ~ Species, value.var="Catch")
catch <- as.matrix(select(catch, -Year))
species <- colnames(catch)

# Format true values
true_orig <- true
B_ts <- as.matrix(select(dcast(data, Year ~ Species, value.var="Biomass"), -Year))
U_ts <- as.matrix(select(dcast(data, Year ~ Species, value.var="ExploitRate"), -Year))
BBMSY_ts <- as.matrix(select(dcast(data, Year ~ Species, value.var="bbmsy"), -Year))
UUMSY_ts <- as.matrix(select(dcast(data, Year ~ Species, value.var="uumsy"), -Year))
true <- NULL
true$B_ts <- B_ts
true$U_ts <- U_ts
true$Bmsy <- unique(data$Bmsy)
true$Umsy <- unique(data$Umsy)

# Fit model
################################################################################

# True parameters
r_true <- c(0.67, 0.11, 0.60)
k_true <- c(227, 738, 304)
id_true <- rep(1,length(k_true))

# Set priors
r_priors <- cbind(r_true, rep(10, ncol(catch)))
k_priors <- cbind(k_true, rep(1000, ncol(catch)))
id_priors <- cbind(id_true, rep(1, ncol(catch)))

# Plot priors
plot_priors(priors=list(r_priors, k_priors, id_priors),
            true_values=list(r_true, k_true, id_true),
            species=species)


# 1. MSCOM - ID fixed at 1.0 and defaults
fit1 <- fit_mscom(catch=catch, id_fixed=T)
plot_mscom(model=fit1, true=true)

# 2. MSCOM - ID estimated and defaults
fit2 <- fit_mscom(catch=catch, id_fixed=F)
plot_mscom(model=fit2, true=true)

# 3. MSCOM - ID fixed and r/k priors
fit3 <- fit_mscom(catch=catch, id_fixed=T, r_priors=r_priors, k_priors=k_priors)
plot_mscom(model=fit3, true=true)

# 4. MSCOM - ID fixed and r priors and K defaults
fit4 <- fit_mscom(catch=catch, id_fixed=T, r_priors=r_priors)
plot_mscom(model=fit4, true=true)


# Real data
################################################################################

# Read RAM Legacy Database
load(file.path(ramdir, "ram_multispecies_fisheries.Rdata"))

# Get DTS data
names(data)
dts <- data[["US West Coast DTS"]]
dts_catch <- as.matrix(select(dts$catch_use, -year))
dts_catch1 <- as.matrix(select(dts$catch_use, -year, -SSTHORNHPCOAST))

# Prep true values
true <- NULL
true$B_ts <- as.matrix(dts$biomass)
true$U_ts <- as.matrix(dts$er)
true$BBMSY_ts <- as.matrix(dts$bbmsy)
true$UUMSY_ts <- as.matrix(dts$ffmsy)


# 1. MSCOM - ID estimated and default
fit5 <- fit_mscom(catch=dts_catch, id_fixed=F)
plot_mscom1(model=fit5, years=dts$catch_use$year, stocks=colnames(dts_catch), true=true)

# 2. MSCOM - ID estimated and default (no shortspine thorny)
fit6 <- fit_mscom(catch=dts_catch1, id_fixed=F)
plot_mscom1(model=fit6, years=dts$catch_use$year, stocks=colnames(dts_catch1), true=true)





