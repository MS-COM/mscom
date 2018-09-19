
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Packages
library(plyr)
library(dplyr)

# Directories
datadir <- "data/priors"

# Read meta-analytics r values
r_vals <- read.csv(file.path(datadir, "neubauer_etal_2013_r_values.csv"), as.is=T)


# Format data
################################################################################

# Exponentiate (note that this is not correct for SDs though)
r_vals[,4:ncol(r_vals)] <- exp(r_vals[,4:ncol(r_vals)])

# SD of stock r's
r_sd_stock <- sd(r_vals$r_stock_mu)

# Summarize by family
r_vals_fam <- r_vals %>% 
  select(family, r_fam_mu) %>% 
  unique() %>% 
  rename(r_mu=r_fam_mu) %>% 
  mutate(r_sd=sd(r_mu)) %>% 
  arrange(desc(r_mu))

# Export data
write.csv(r_vals_fam, file.path(datadir, "r_priors_by_family.csv"), row.names=F)

# Plot data
################################################################################


# Setup figure 
figname <- "Fig1_r_meta_analytic_means.png"
png(paste(datadir, figname, sep="/"), width=6.5, height=4, units="in", res=600)
par(mfrow=c(1,1), mar=c(8, 4.5, 0.5, 0.5), mgp=c(3.3,1,0))

# Plot figure
plot(1:nrow(r_vals_fam), r_vals_fam$r_mu, bty="n", las=1, pch=16, 
     xaxt="n", xlab="", ylab="Intrinsic growth rate, r", ylim=c(0,1.2))
axis(side=1, at=1:nrow(r_vals_fam), labels=r_vals_fam$family, las=2, cex.axis=0.9)

# Off
dev.off()

