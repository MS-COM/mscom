
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
datadir <- "data/simulated/data"
codedir <- "code/helper_functions"
outdir <- "data/output"

# Load scenarios
load(file.path(datadir, "simulated_multispecies_fisheries.Rdata"))

# Format "true" values
sims1 <- sims %>% 
  rename(biomass_true=biomass, catch_true=catch, er_true=er, q_true=q, effort_true=effort, sat_true=sat, uumsy_true=uumsy, bbmsy_true=bbmsy)


# Format data
################################################################################

# Output files
outfiles <- list.files(outdir, pattern=".Rdata")

# Loop through output files
i <- 1
for(i in 1:length(outfiles)){
  
  # Load file
  file_i <- outfiles[i]
  load(file.path(outdir, file_i))
  iter_i <- gsub(".Rdata", "", file_i)
  print(i)
  
  # Extract predictions
  preds <- out$preds
  
  # Extract truth for scenario
  obs <- filter(sims1, iter==iter_i)
  
  # Merge scenario details, true values, and predicted values
  combo <- obs %>% 
    left_join(preds, by=c("stock", "year")) %>% 
    select(-catch) %>% 
    left_join(select(scenarios, -c(fishery, scenario_id)), by=c("iter"="iter_id"))
  
  # Merge files
  if(i==1){data <- combo}else{data <- rbind(data, combo)}
  
}

# Format data
data1 <- data %>% 
  mutate(cmsy_diff=(bbmsy_cmsy-bbmsy_true)/bbmsy_true*100,
         mscmsy_diff=(bbmsy-bbmsy_true)/bbmsy_true*100)

# Export
################################################################################

# Export results
write.csv(data1, file.path(outdir, "simtest_results.csv"), row.names=F)



# Visualize performance
data2 <- filter(data1, year==2015)
boxplot(data2[,c("cmsy_diff", "mscmsy_diff")], frame=F, las=1, lty=1, ylim=c(-100, 400),
        names=c("cMSY", "MS-cMSY"), ylab="Percent error in estimate")
lines(x=c(0.5,3.5), y=c(0,0), lty=2, col="red")

# 
plot(bbmsy ~ bbmsy_true, data2, xlim=c(0,2), ylim=c(0,2))



