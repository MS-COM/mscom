
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
outdir <- "data/output"


# Format data
################################################################################

# Load MS-cMSY results
load(file.path(outdir, "pelagic_longline_simtest.Rdata"))

# Load simulation scenarios
load(file.path(datadir, "scenarios_pelagic_longline.Rdata"))
scenarios <- scenarios %>% 
  rename(dyn=PopDyn, id=DataStart, ed=EffDyn, maxf=MaxF, nyr=Nyears, effsd=EffSD)

# Simulation output: 1 file per iteration
# These were run through MS-cMSY in the same order
files <- list.files(datadir, pattern="byScenario")

# Loop through iterations: i <- 1; j <- 1
for(i in 1:length(files)){
  
  # MS-cMSY predictions
  preds_i <- output_all[[1]]
  
  # Load simulation output (to get true B/BMSY values)
  file_i <- files[i]
  load(file.path(datadir, file_i))
  
  # Loop through scenarios
  nscenarios <- length(byScen) # 1 element per scenario
  for(j in 1:nscenarios){
    
    # Observed status
    true <- byScen[[j]] %>% 
      filter(Year==max(Year)) %>% 
      select(Species, BBmsy) %>% 
      rename(stock=Species, true=BBmsy)
    
    # cMSY status
    cmsy <- sapply(preds_i[[j]]$bbmsy_v_median, cbind)
    colnames(cmsy) <- preds_i[[j]]$stocks
    cmsy_end <- cmsy[nrow(cmsy),]
    cmsy_end <- data.frame(stock=colnames(cmsy), cmsy=cmsy_end)
    
    # MS-cMSY status
    mscmsy <- sapply(preds_i[[j]]$bbmsy_vv_median, cbind)
    colnames(mscmsy) <- preds_i[[j]]$stocks
    mscmsy_end <- mscmsy[nrow(mscmsy),]
    mscmsy_end <- data.frame(stock=colnames(mscmsy), mscmsy=mscmsy_end)
    
    # Add scenario meta-data
    sdata <- true %>%
      left_join(cmsy_end, by="stock") %>% 
      left_join(mscmsy_end, by="stock")
    sdata1 <- cbind(scenarios[j,], sdata) # Raises WARNINGS that don't matter
    
    # Merge data
    if(j==1){data_i <- sdata1}else{data_i <- rbind(data_i, sdata1)}
    
  }
  
  # Format data 
  fdata <- data_i %>% 
    mutate(simfile=file_i, iter=i) %>% 
    select(simfile, dyn:effsd, iter, everything())
  
  # Merge
  if(i==1){data <- fdata}else{data <- rbind(data, fdata)}
  
}

# Format data
data <- data %>% 
  arrange(dyn, id, ed, maxf, effsd, iter) %>% 
  mutate(sim_id=paste(dyn, id, ed, maxf, effsd, iter, sep="-"),
         stock_id=paste(dyn, id, ed, maxf, effsd, iter, stock, sep="-"),
         cmsy_diff=(cmsy-true)/true*100,
         mscmsy_diff=(mscmsy-true)/true*100) %>% 
  select(sim_id, stock_id, dyn:mscmsy_diff)

# Export results
write.csv(data, file.path(outdir, "pelagic_longline_simtest_results.csv"), row.names=F)

# Visualize performance
boxplot(data[,c("cmsy_diff", "mscmsy_diff")], frame=F, las=1, lty=1, ylim=c(-300, 300),
        names=c("cMSY", "MS-cMSY"), ylab="Percent error in estimate")
lines(x=c(0.5,3.5), y=c(0,0), lty=2, col="red")






