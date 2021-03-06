
# Mode
mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Calculate r prior
calc_r_priors <- function(res){
  for(i in 1:length(res)){
    if(res[i]=="High"){r_prior <- c(0.6,1.5)}
    if(res[i]=="Medium"){r_prior <- c(0.2,0.8)}
    if(res[i]=="Low"){r_prior <- c(0.05,0.5)}
    if(res[i]=="Very low"){r_prior <- c(0.015,0.1)}
    if(i==1){r_priors <- r_prior}else{r_priors <- rbind(r_priors, r_prior)}
  }
  colnames(r_priors) <- c("r_lo", "r_hi")
  rownames(r_priors) <- NULL
  return(r_priors)
}

# Calculate k priors
calc_k_priors <- function(C_mat){
  cmax <- apply(C_mat, 2, max)
  k_lo <- cmax
  k_hi <- cmax * 50
  k_priors <- cbind(k_lo=k_lo, k_hi=k_hi)
  return(k_priors)
}

# Calculate initial year depletion prior
calc_sat1_priors <- function(C_mat){
  nyr <- nrow(C_mat)
  nstocks <- ncol(C_mat)
  if(nyr>50){
    s1_priors <- cbind(s1_lo=rep(0.7, nstocks), s1_hi=rep(1.0, nstocks))
  }else{
    s1_priors <- cbind(s1_lo=rep(0.4, nstocks), s1_hi=rep(1.0, nstocks))
  }
}

# Calculate final year depletion prior
# >0.8 - 0.4-0.8
# 0.5-0.8 - 0.2-0.6
# <0.5 - 0.01-0.04
calc_sat2_priors <- function(C_mat){
  c_ends <- C_mat[nrow(C_mat),]
  c_maxs <- apply(C_mat, 2, max)
  c_ratios <- c_ends / c_maxs
  b_catgs <- cut(c_ratios, breaks=c(0,0.05, 0.15, 0.35, 0.5, 0.8, 1),
                 labels=c("very very very low", "very very low", "very low", "low", "medium", "high"))
  s_prior_key <- matrix(data=c(0.01, 0.1, # very very very low
                                0.01, 0.2, # very very low
                                0.01, 0.3, # very low
                                0.01, 0.4, # low
                                0.2, 0.6, # medium
                                0.4, 0.8), # high
                         ncol=2, byrow=T, dimnames=list(NULL, c("s_lo", "s_hi")))
  s2_priors <- s_prior_key[b_catgs,]
  return(s2_priors)
}


# Multi-species catch-only model
# For testing: catch<-catch; years<-yrs; stocks<-species; res<-res; id_fixed <- F; npairs <- 1000
fit_mscmsy <- function(catch, years, stocks, res, id_fixed, npairs=10000){
  
  # Time series info
  nyrs <- length(years)
  nstocks <- length(stocks)
  
  # Calculate r and k priors
  # R prior based on resilience; K prior based on max catch and r prior
  r_priors <- calc_r_priors(res)
  k_priors <- calc_k_priors(catch)
  r_priors_ln <- log(r_priors)
  k_priors_ln <- log(k_priors)
  
  # Calculate saturation priors
  s1_priors <- calc_sat1_priors(catch)
  s2_priors <- calc_sat2_priors(catch)

  # Randomly sample r-k pairs in log-space
  npairs <- npairs
  ri <- sapply(1:nstocks, function(x) exp(runif(npairs, r_priors_ln[x,1], r_priors_ln[x,2])))
  ki <- sapply(1:nstocks, function(x) exp(runif(npairs, k_priors_ln[x,1], k_priors_ln[x,2])))
  
  # Lists to hold viable r/k pairs and associated biomass/exploitation trajectories
  id_rk_try <- list()
  id_rk_v <- list()
  b_mats_v <- list()
  bbmsy_mats_v <- list()
  er_mats_v <- list()
  
  # Loop through stocks to identify viable r/k pairs and trajectoris
  for(i in 1:nstocks){
    
    # Get info
    stock <- stocks[i]
    c_vec <- catch[,i]
    # print(stock)
    
    # Initial depletions to evaluate
    if(id_fixed==T){
      ids <- 1
    }else{
      s1_prior <- s1_priors[i,]
      ids <- seq(s1_prior[1], s1_prior[2], 0.1)
    }
    
    # Get r/k pairs to evaluate
    rs <- ri[,i]
    ks <- ki[,i]
    rk_pairs <- cbind(r=rs, k=ks, viable=rep(NA, npairs))
    
    # Build ID/r/k combos to evaluate
    id_rk_combos <- as.data.frame(do.call("rbind",
                       lapply(ids, function(x) cbind(id=rep(x, nrow(rk_pairs)), rk_pairs))))
    
    # Loop through r/k pairs to see if viable
    p <- 0.2
    sigmaP <- 0.1
    b_mat <- matrix(data=NA, nrow=nyrs, ncol=nrow(id_rk_combos))
    for(j in 1:nrow(id_rk_combos)){
      id <- id_rk_combos$id[j]
      r <- id_rk_combos$r[j]
      k <- id_rk_combos$k[j]
      b_mat[1,j] <- k * id
      for(yr in 2:nyrs){
        b_mat[yr,j] <- b_mat[yr-1,j] +  r*b_mat[yr-1,j]/p*(1-(b_mat[yr-1,j]/k)^p) - c_vec[yr-1]
        # b_mat[yr,j] <- b_mat[yr-1,j] +  r*b_mat[yr-1,j]/p*(1-(b_mat[yr-1,j]/k)^p)*exp(rnorm(1,0,sigmaP)) - c_vec[yr-1]
      }
    }
    
    # Create saturation matrix
    s_mat <- t(t(b_mat)/ks)
    
    # Reduce to viable r/k pairs and trajectories
    check0 <- apply(b_mat, 2, function(x) sum(x<0, na.rm=T)==0) # check B doesn't go below 0
    s2_lo <- s2_priors[i,1]
    s2_hi <- s2_priors[i,2]
    s2_vec <- s_mat[nrow(s_mat),]
    checkS <- s2_vec > s2_lo & s2_vec < s2_hi & !is.na(s2_vec) # check final yr saturation inside prior
    viable <- check0 & checkS # merge positive biomass and fina saturation checks
    id_rk_combos[,"viable"] <- viable
    nviable <- sum(viable)
    b_mat_viable <- b_mat[,viable]
    id_rk_viable <- id_rk_combos[viable,]
    
    # Derive B/BMSY
    bmsy <- id_rk_viable[,"k"] * (1 / (p+1))^(1/p)
    bbmsy_mat_viable <- t(t(b_mat_viable) / bmsy)
    
    # # Plot viable r/k pairs
    # par(mfrow=c(1,1))
    # plot(k ~ r, rk_pairs, log="xy", bty="n", las=1, pch=15,
    #      xlim=r_priors[i,], ylim=k_priors[i,], xlab="r", ylab="k", col="gray80")
    # points(rk_viable[,1], rk_viable[,2], pch=15, col="black")
    # 
    # # Plot viable biomass trajectories
    # plot(b_mat_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, max(b_mat_viable, na.rm=T)), xlab="Year", ylab="Biomass")
    # for(k in 1:ncol(b_mat_viable)){lines(x=yrs, y=b_mat_viable[,k], col="grey80")}
    #
    # # Plot viable B/BMSY trajectories
    # plot(bbmsy_mat_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, max(bbmsy_mat_viable, na.rm=T)), xlab="Year", ylab="B/BMSY")
    # for(k in 1:ncol(bbmsy_mat_viable)){lines(x=yrs, y=bbmsy_mat_viable[,k], col="grey80")}

    # Calculate exploitation rate
    # I name the rows and columns so that I can validate covariance matrix below
    er_mat_viable <- c_vec / b_mat_viable
    rownames(er_mat_viable) <- yrs
    colnames(er_mat_viable) <- paste0(LETTERS[i], 1:ncol(er_mat_viable))
    
    # Plot viable exploitation trajectories
    # plot(er_mat_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, 1), xlab="Year", ylab="Exploitation rate")
    # for(k in 1:ncol(er_mat_viable)){lines(x=yrs, y=er_mat_viable[,k], col="grey80")}
    
    # Record results
    id_rk_v[[i]] <- id_rk_viable
    b_mats_v[[i]] <- b_mat_viable
    bbmsy_mats_v[[i]] <- bbmsy_mat_viable
    er_mats_v[[i]] <- er_mat_viable
    
  }
  
  # Measure correlation
  ###################################################
  
  # Scale exploitation rates to mean:   apply(er_mats_v_scaled[[]], 2, mean)
  er_mats_v_scaled <- sapply(er_mats_v, function(x) t(t(x) / apply(x, 2, mean)))
  
  
  # Calculate variance ratio
  vr <- function(mat){var(apply(mat, 1, sum)) / sum(apply(mat, 2, var))}
  

  # Build grid of parameter combos to evaluate
  id_rk_v_n_all <- sapply(id_rk_v, function(x) 1:nrow(x)) # number of viable r/k pairs for each stock
  combos <- expand.grid(id_rk_v_n_all)
  colnames(combos) <- paste0("index", 1:ncol(combos))
  combos_sample <- sample_n(combos, nrow(combos)*0.1)
  
  # Loop through combos
  vrs <- NA
  print("Calculating variance ratios:")
  pb <- txtProgressBar(min=0, max=nrow(combos_sample), initial = 0) 
  for(i in 1:nrow(combos_sample)){
    
    # Build data
    setTxtProgressBar(pb,i)
    inds <- as.numeric(combos_sample[i,])
    erdata <- do.call("cbind", lapply(1:length(inds), function(x) er_mats_v_scaled[[x]][,inds[x]]))
    
    # Calculate and record variance ratio
    vrs[i] <- vr(erdata)
    
  }
  
  # Add variance ratio to r/K/ID combo key
  combos_sample$vr <- vrs  # hist(combos_sample$vr)
  
  # Identify top X%
  top_p <- 0.01
  top_n <- ceiling(nrow(combos_sample) * top_p)
  top_n <- 5
  top_combos <- combos_sample %>% 
    arrange(desc(vr)) %>%
    slice(1:top_n)
  
  # Get biomass trajectories of top 10%
  # (also sneak in calculation of cMSY prediction)
  bbmsy_v_meds <- NULL
  er_mats_vv <- NULL
  bbmsy_mats_vv <- NULL
  bbmsy_vv_meds <- NULL
  for(i in 1:length(b_mats_v)){
    vv_index <- unlist(top_combos[,paste0("index", i)])
    er_mat_v <- er_mats_v[[i]]
    er_mat_vv <- er_mat_v[,vv_index]
    er_mats_vv[[i]] <- er_mat_vv
    bbmsy_mat_v <- bbmsy_mats_v[[i]]
    bbmsy_mat_vv <- bbmsy_mat_v[,vv_index]
    bbmsy_mats_vv[[i]] <- bbmsy_mat_vv
    bbmsy_v_med <- apply(bbmsy_mat_v, 1, median)
    bbmsy_vv_med <- apply(bbmsy_mat_vv, 1, median)
    bbmsy_v_meds[[i]] <- bbmsy_v_med
    bbmsy_vv_meds[[i]] <- bbmsy_vv_med
  }

  # Things to return
  out <- list(stocks=stocks,
              yrs=yrs,
              r_priors=r_priors,
              k_priors=k_priors,
              s1_priors=s1_priors,
              s2_priors=s2_priors,
              id_try=ids,
              r_try=ri, 
              k_try=ki, 
              id_rk_v=id_rk_v, 
              b_v=b_mats_v, 
              bbmsy_v=bbmsy_mats_v,
              bbmsy_v_median=bbmsy_v_meds,
              er_v=er_mats_v,
              er_vv=er_mats_vv,
              top_corr=top_corr,
              bbmsy_vv=bbmsy_mats_vv,
              bbmsy_vv_median=bbmsy_vv_meds)
  return(out)
  
}

# Plot MS-SRA results
plot_mssra <- function(out, true){
  
  # Loop through species
  spp <- species
  nspp <- length(species)
  par(mfcol=c(3,nspp), mar=c(3,4,2,1), mgp=c(2.5,0.7,0), xpd=NA)
  for(i in 1:length(species)){
    
    # Subset data
    yrs <- out$yrs
    ids <- out$id_try
    rs <- out$r_try[,i]
    ks <- out$k_try[,i]
    id_rk_viable <- out$id_rk_v[[i]]
    b_viable <- out$b_v[[i]]
    bbmsy_viable <- out$bbmsy_v[[i]]
    er_viable <- out$er_v[[i]]
    s1_priors <- out$s1_priors
    s2_priors <-out$s2_priors
    r_priors <- out$r_priors
    k_priors <- out$k_priors
    er_vv <- out$er_vv[[i]]
    bbmsy_vv <- out$bbmsy_vv[[i]]
    bbmsy_vv_median <- out$bbmsy_vv_median[[i]]
    top_corr <- out$top_corr
    
    # Plot r/k pairs
    #########################################
    
    # Plot r/k pairs
    plot(ks ~ rs, log="xy", type="n", bty="n", las=1, pch=15, col="gray80", xpd=NA,
         xlim=r_priors[i,], ylim=k_priors[i,], xlab="r", ylab="k", main=spp[i])
    
    # Add viable r/k pairs
    # Potentially reduce this to unique r/k pairs
    # There could be redundancy when evaluating multiple IDs
    points(id_rk_viable$r, id_rk_viable$k, pch=15, col="grey70")
    
    # Add most highly correlated pairs
    id_rk_v_ind <- unlist(top_corr[,paste0("index", i)])
    rk_corr <- subset(id_rk_viable, index %in% id_rk_v_ind)
    points(x=rk_corr$r, y=rk_corr$k, pch=15, col=freeR::tcolor("darkorange", 0.6))
    
    # # Add most common highly correlated pair
    # rk_mode <- mode(rk_v_ind)
    # rk_corr <- subset(rk_viable, index==rk_mode)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="green")
    
    # # Add repeatedly highly correlated r/k pairs
    # rk_corr <- subset(rk_viable, ncorr>=5)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="darkorange")
  
    # Add true value
    points(x=true$r[i], y=true$k[i], pch=15, col="red", cex=1.3)
    
    # Add legend
    if(i==1){
      legend("bottomright", bty="n", pch=15, pt.cex=1.3, cex=0.9, 
             col=c("grey70", "darkorange", "red"),
             legend=c("Viable", "Effort correlated", "True value"))
    }
    
    # Plot BBMSY trajectories
    #########################################
    
    # Plot BBMSY trajectories
    ymax <- freeR::ceiling1(max(bbmsy_viable, true$bbmsy_ts[,i+1], na.rm=T), 0.5)
    plot(bbmsy_viable[,1] ~ yrs, type="n", bty="n", las=1,
         ylim=c(0, ymax), xlab="Year", ylab=expression("B/B"["MSY"]))
    for(k in 1:ncol(bbmsy_viable)){lines(x=yrs, y=bbmsy_viable[,k], col="grey70")}
    for(k in 1:ncol(bbmsy_vv)){lines(x=yrs, y=bbmsy_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
    lines(x=yrs, y=bbmsy_vv_median, lwd=1.5, col="black")
    lines(x=yrs, y=true$bbmsy_ts[,i+1], lwd=1.5, col="red", lty=3) # true B/BMSY
    lines(x=c(0, max(yrs)), y=c(0.5, 0.5), lty=3)
    lines(x=yrs, y=apply(bbmsy_viable,1,median), lwd=1.5, lty=3, col="green")
    
    # Add legend
    if(i==1){
      legend("bottomright", bty="n", lty=c(1,1,3), lwd=1.5, cex=0.9, 
             col=c("grey70", "darkorange", "red"),
             legend=c("Viable", "Effort correlated", "True value"))
    }
    
    # Plot biomass trajectories
    #########################################
    
    # # Plot biomass trajectories
    # plot(b_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, max(b_viable)), xlab="Year", ylab="Biomass")
    # for(k in 1:ncol(b_viable)){lines(x=yrs, y=b_viable[,k], col="grey80")}
    # lines(x=yrs, y=true$b_ts[,i+1], lwd=1.5, col="red")
    
    # Plot exploitation trajectories
    #########################################
    
    # Plot exploitation trajectories
    plot(er_viable[,1] ~ yrs, type="n", bty="n", las=1,
         ylim=c(0, 1), xlab="Year", ylab="Exploitation rate")
    for(k in 1:ncol(er_viable)){lines(x=yrs, y=er_viable[,k], col="grey80")}
    for(k in 1:ncol(er_vv)){lines(x=yrs, y=er_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
    lines(x=yrs, y=true$er_ts[,i+1], lwd=1.5, col="red", lty=3)
    
  }
  
  
}
