
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Turn off sci notation
options(scipen=999)

# Packages
library(rio)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# Datalimited2
devtools::install_github("cfree14/datalimited2")
library(datalimited2)

# Directories
datadir <- "data/test_cases/data"
plotdir <- "figures"
codedir <- "code/helper_functions"

# Read data
data <- import(file.path(datadir, "gulf_of_thailand_trawl_catch_data.csv"))

# Read helper functions
# helpers <- list.files(codedir)
# ignore <- sapply(1:length(helpers), function(x) source(file.path(codedir, helpers[x])))


# Make predictions
################################################################################

# Reshape data
data_w <- na.omit(dcast(data, year ~ comm_name, value.var="catch"))

# Look up resilience
spp <- unique(data$sci_name)
res <- datalimited2::resilience(spp)

# Make predictions
out <- datalimited2::ms_cmsy(catch=data_w, 
               stocks=c("Bigeye scad", "Blue swimming crab", "Largehead hairtail"),
               res=c("High", "High", "Low"),
               id_fixed=F, npairs=5000)

# Plot data
datalimited2::plot_dlm(out)

# Export predictions
save(out, file=file.path(datadir, "got_ms_cmsy_predictions.Rdata"))



