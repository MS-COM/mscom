
# Setup
################################################################################

# Clear workspace
rm(list = ls())

# Packages
library(freeR)
library(plyr)
library(dplyr)

# Directories
datadir <- "data/priors"

# Read meta-analytics r values
r_vals_orig <- read.csv(file.path(datadir, "neubauer_etal_2013_r_values.csv"), as.is=T)


# Format data
################################################################################

# Format r values
r_vals <- r_vals_orig %>% 
  # Remove SDs
  select(-c(r_fam_sd, r_spp_sd, r_stock_sd)) %>%
  # Rename columns
  rename(r_fam=r_fam_mu, r_spp=r_spp_mu, r_stock=r_stock_mu) %>% 
  # Fix scientific names and exponentiate means
  mutate(r_fam=exp(r_fam),
         r_spp=exp(r_spp),
         r_stock=exp(r_stock),
         r_orig=exp(r_orig),
         family=revalue(family, c("Scorpaenidae"="Sebastidae")),
         species=revalue(species, c("Chrysophrys auratus"="Pagrus auratus",
                                    "Clupea pallasii"="Clupea pallasii pallasii",
                                    "Neoplatycephalus richardsoni"="Platycephalus richardsoni",
                                    "Reinhardtius stomias"="Atheresthes stomias",
                                    "Solea vulgaris"="Solea solea")))

# Add resilience values
res <- datalimited2::resilience(r_vals$species)
r_vals1 <- r_vals %>% 
  left_join(res, by="species")

# Check names
# Species names correct?
freeR::suggest_names(r_vals$species)
# Family names correct?
r_taxa <- freeR::taxa(r_vals$species)
r_vals$family[!r_vals$family %in% r_taxa$family]

# SD of stock r's
r_sd_stock <- sd(r_vals$r_stock)

# Summarize by family
r_vals_fam <- r_vals %>% 
  select(family, r_fam) %>% 
  unique() %>% 
  rename(r_mu=r_fam) %>% 
  mutate(r_sd=sd(r_mu)) %>% 
  arrange(desc(r_mu))

# Export data
write.csv(r_vals_fam, file.path(datadir, "r_priors_by_family.csv"), row.names=F)

# Plot data
################################################################################

# Setup figure 
figname <- "Fig1_r_meta_analytic_means.png"
png(paste(datadir, figname, sep="/"), width=6.5, height=4.5, units="in", res=600)
par(mfrow=c(1,1), mar=c(8, 4.5, 0.5, 0.5), mgp=c(3.3,1,0))

# Plot figure
rmin <- pmax(0, r_vals_fam$r_mu - r_vals_fam$r_sd)
rmax <- r_vals_fam$r_mu + r_vals_fam$r_sd
plot(1:nrow(r_vals_fam), r_vals_fam$r_mu, type="n", bty="n", las=1, pch=16, 
     xaxt="n", xlab="", ylab="Intrinsic growth rate, r", ylim=c(0,1.4))
axis(side=1, at=1:nrow(r_vals_fam), labels=r_vals_fam$family, las=2, cex.axis=0.9)

# Add family SDs
sapply(1:length(rmin), function(i) lines(x=c(i,i), y=c(rmin[i], rmax[i]), col="grey60"))

# Add stock means
for(i in 1:nrow(r_vals_fam)){
  fam <- r_vals_fam$family[i]
  fdata <- filter(r_vals1, family==fam)
  points(x=rep(i, nrow(fdata)), y=fdata$r_stock, pch=4, col="grey30", cex=0.8)
}

# Legend
legend("topright", bty="n", cex=0.9,
       legend=c("Stock r", "Family r"), pch=c(4, 16), col=c("grey30", "black"))

# Add family means
points(x=1:nrow(r_vals_fam), y=r_vals_fam$r_mu, pch=16)

# Off
dev.off()

