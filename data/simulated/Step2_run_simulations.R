
# Clean
rm(list = ls())

# Setup
################################################################################

# Directories
datadir <- "data/simulated/data"


# Setup simulations
################################################################################

# Read fishery info
fdata <- read.csv(file.path(datadir, "demersal_trawl_fishery_10spp.csv"), as.is=T)
fisheries <- sort(unique(fdata$fishery))
  
# Factors
niter <- 1
id_vec <- c("unfished", "fished")
ed_vec <- c("one-way", "two-way")
sigmaP_vec <- c(0.05, 0.15)
sigmaF_vec <- c(0, 0.10)
nstocks_vec <- c(2,5,10)
nyrs_vec <- c(25, 50)

# Final year?
yr_end <- 2015

# Build scenario grid
scenarios_raw <- expand.grid(fishery=fisheries, id=id_vec, ed=ed_vec, sigmaP=sigmaP_vec, 
                             sigmaF=sigmaF_vec, nstocks=nstocks_vec, nyrs=nyrs_vec, iter=1:niter, stringsAsFactors=F)

# Format scenario grid
scenarios <- scenarios_raw %>% 
  arrange(fishery, id, ed, sigmaP, sigmaF, nstocks, nyrs, iter) %>% 
  mutate(scenario_id=paste(id, ed, sigmaP, sigmaF, nstocks, nyrs, sep="_"), 
         iter_id=paste(fishery, id, ed, sigmaP, sigmaF, nstocks, nyrs, iter, sep="_")) %>% 
  select(fishery, scenario_id, iter_id, iter, everything())

# Simulation loop
################################################################################

# Loop through scenarios
i <- 1
pb <- txtProgressBar(min = 0, max = nrow(scenarios), style = 3)
for(i in 1:nrow(scenarios)){
    
  # Scenario
  scenario <- scenarios[i,]
  iter_id <- scenario$iter_id
  scenario_id <- scenario$scenario_id
  setTxtProgressBar(pb, i)
    
  # Scenario details
  fshry <- scenario$fishery
  nstocks <- scenario$nstocks
  sdata <- fdata %>% filter(fishery==fshry) %>%  slice(1:nstocks)
  stocks <- sdata$comm_name
  rs <- sdata$r
  ks <- sdata$k
  nyr <- scenario$nyrs
  id <- scenario$id
  ed <- scenario$ed
  sigmaP <- scenario$sigmaP
  sigmaF <- scenario$sigmaF
  iter <- scenario$iter
  
  # Set start year based on fished state
  yr1 <- ifelse(id=="unfished", 1, 21)
  tyr <- ifelse(id=="unfished", nyr, nyr+yr1-1)
  
  # Set rho (autocorrelation) based on sigmaF
  rho <- 0.4
  
  # Loop through species
  pop <- lapply(1:nstocks, function(x){
    
    # Pop parameters
    stock <- stocks[x]
    r <- rs[x]
    k <- ks[x]
    p <- 0.2
    msy <- (r*k)/((p+1)^((p+1)/p))
    bmsy <- k * (1/(p+1))^(1/p)
    umsy <- r/p*(1-1/(p+1))
    
    # Effort dynamics
    ############################################
    
    # Derive effort
    umax <- umsy * 1.6 # how much does effort exceed MSY by?
    
    # Effort deviations
    dev <- rnorm(tyr, mean=-(sigmaF ^ 2)/2, sd=sigmaF)
    qdev <- unlist(sapply(1:length(dev), function(x){
      if(x==1) w <- dev[x]
      if(x>1) w <- dev[x-1] + dev[x]
      return(w)
    }))

    # One-way trip
    if(ed=="one-way"){
      et <- seq(0, 1, length=tyr)
      qt <- umax/max(et) * exp(qdev)
      ut <- qt*et
    }
    
    # Two-way trip
    if(ed=="two-way"){
      et1 <- seq(0,to=1,length=ceiling(tyr/2))
      et2 <- seq(max(et1),to=0.5,length=floor(tyr/2))
      et <- c(et1,et2)
      qt <- umax/max(et) * exp(qdev)
      ut <- qt*et
    }
    
    # Biomass dynamics
    ############################################
    
    # Process deviations
    dev <- rnorm(tyr, mean= -(sigmaP ^ 2)/2, sd=sigmaP)
    pdev <- rep(NA, tyr)
    pdev[1] <- dev[1]
    for(t in 2:length(dev)){
      pdev[t] <- pdev[t-1] * rho + sqrt(1 - rho ^ 2) * dev[t]
    }

    # Setup containers
    bt <- ct <- pt <- rep(NA, tyr)
    
    # Insert first year
    bt[1] <- k
    ct[1] <- ut[1] * bt[1]
    pt[1] <- (r/p)*bt[1]*(1-(bt[1]/k)^p)
    
    # Simulate following years
    for(t in 2:tyr){
      ct[t] <- ut[t-1] * bt[t-1]
      pt[t] <- (r/p)*bt[t-1]*(1-(bt[t-1]/k)^p)
      # bt[t] <- (bt[t-1] + pt[t] - ct[t-1]) * exp(pdev[t]) # Merrill's version
      bt[t] <- bt[t-1] + pt[t] * exp(pdev[t]) - ct[t-1] # Chris' version
    }
    
    # Trim to first year (if fished)
    bt_out <- bt[yr1:length(bt)]
    ct_out <- ct[yr1:length(bt)]
    ut_out <- ut[yr1:length(bt)]
    et_out <- et[yr1:length(bt)]
    qt_out <- qt[yr1:length(bt)]
    
    # Collate results
    out <- data.frame(fishery=fshry,
                      scenario=scenario_id,
                      iter=iter_id,
                      stock=stock,
                      year=(yr_end-length(bt_out)+1):yr_end,
                      biomass=bt_out, 
                      catch=ct_out,
                      er=ut_out,
                      q=qt_out, 
                      effort=et_out,
                      sat=bt_out/k,
                      uumsy=ut_out/umsy,
                      bbmsy=bt_out/bmsy)
    
    # Return
    return(out)
    
  })
  
  # Collapse list of stock histories into single dataframe
  pop_out <- do.call(rbind, pop)
  
  # Merge simulation results
  if(i==1){sims <- pop_out}else{sims <- rbind(sims, pop_out)}
  
}
close(pb)

# Remove iterations w/ negative biomass
prob_iter <- as.character(sort(unique(sims$iter[!is.na(sims$biomass) & sims$biomass <0])))
sims_final <- filter(sims, !iter %in% prob_iter)

# Inspect distribution of final B/BMSY values
stats <- sims_final %>% 
  group_by(iter, stock) %>% 
  summarize(bbmsy=bbmsy[year==yr_end])
hist(stats$bbmsy, xlim=c(0,2), breaks=seq(0,2,0.1))

# Export simulations
sims <- sims_final
save(scenarios, sims, fdata, file=file.path(datadir, "simulated_multispecies_fisheries.Rdata"))
  