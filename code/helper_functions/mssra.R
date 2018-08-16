
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

# Calculate k prior
calc_k_priors <- function(C_mat, r_priors){
  for(i in 1:ncol(C_mat)){
    catch <- C_mat[,i]
    r_prior <- r_priors[i,]
    r_lo <- r_prior[1]
    r_hi <- r_prior[2]
    c_max <- max(catch)
    k_lo <- c_max / r_hi
    k_hi <- 12 * c_max / r_lo
    k_prior <- c(k_lo, k_hi)
    if(i==1){k_priors <- k_prior}else{k_priors <- rbind(k_priors, k_prior)}
  }
  colnames(k_priors) <- c("k_lo", "k_hi")
  rownames(k_priors) <- NULL
  return(k_priors)
}


# Calculate initial year depletion prior
calc_sat1_priors <- function(C_mat){
  nyr <- nrow(C_mat)
  
  
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
  s_priors <- s_prior_key[b_catgs,]
  return(s_priors)
}


# Multi-species catch-only model
# For testing: catch<-catch; years<-yrs; stocks<-species; res<-res
fit_mssra <- function(catch, years, stocks, res){
  
  # Time series info
  nyrs <- length(years)
  nstocks <- length(stocks)
  
  # Calculate r and k priors
  # R prior based on resilience; K prior based on max catch and r prior
  r_priors <- calc_r_priors(res)
  k_priors <- calc_k_priors(catch, r_priors)
  r_priors_ln <- log(r_priors)
  k_priors_ln <- log(k_priors)
  
  # Calculate final saturation priors
  s2_priors <- calc_sat2_priors(catch)
  
  # Randomly sample r-k pairs in log-space
  npairs <- 10000
  ri <- sapply(1:nstocks, function(x) exp(runif(npairs, r_priors_ln[x,1], r_priors_ln[x,2])))
  ki <- sapply(1:nstocks, function(x) exp(runif(npairs, k_priors_ln[x,1], k_priors_ln[x,2])))
  
  # Lists to hold viable r/k pairs and associated biomass/exploitation trajectories
  rk_try <- list()
  rk_v <- list()
  b_mats_v <- list()
  bbmsy_mats_v <- list()
  er_mats_v <- list()
  
  # Loop through stocks to identify viable r/k pairs and trajectoris
  for(i in 1:nstocks){
    
    # Get info
    stock <- stocks[i]
    c_vec <- catch[,i]
    print(stock)
    
    # Get r/k pairs to evaluate
    ri1 <- ri[,i]
    ki1 <- ki[,i]
    rk_pairs <- cbind(r=ri1, k=ki1, viable=rep(NA, npairs))
    
    # Loop through r/k pairs to see if viable
    # Currently, population begins at carrying capacity
    b_mat <- matrix(data=NA, nrow=nyrs, ncol=npairs)
    for(j in 1:npairs){
      r <- ri1[j]
      k <- ki1[j]
      b_mat[1,j] <- k
      for(yr in 2:nyrs){
        b_mat[yr,j] <- b_mat[yr-1,j] +  r*b_mat[yr-1,j] * (1-b_mat[yr-1,j]/k) - c_vec[yr-1]
      }
    }
    
    # Create saturation matrix
    s_mat <- t(t(b_mat)/ki1)
    
    # Reduce to viable r/k pairs and trajectories
    check0 <- apply(b_mat, 2, function(x) sum(x<0)==0) # check B doesn't go below 0
    s2_lo <- s2_priors[i,1]
    s2_hi <- s2_priors[i,2]
    s2_vec <- s_mat[nrow(s_mat),]
    checkS <- s2_vec > s2_lo & s2_vec < s2_hi # check final yr saturation inside prior
    viable <- check0 & checkS
    rk_pairs[,"viable"] <- viable
    nviable <- sum(viable)
    b_mat_viable <- b_mat[,viable]
    rk_viable <- rk_pairs[viable,]
    
    # Derive BBMSY
    bmsy <- rk_viable[,"k"] / 2
    bbmsy_mat_viable <- t(t(b_mat_viable) / bmsy)
    
    # # Plot viable r/k pairs
    # par(mfrow=c(1,1))
    # plot(k ~ r, rk_pairs, log="xy", bty="n", las=1, pch=15,
    #      xlim=r_priors[i,], ylim=k_priors[i,], xlab="r", ylab="k", col="gray80")
    # points(rk_viable[,1], rk_viable[,2], pch=15, col="black")
    # 
    # # Plot viable biomass trajectories
    # plot(b_mat_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, max(b_mat_viable)), xlab="Year", ylab="Biomass")
    # for(k in 1:ncol(b_mat_viable)){lines(x=yrs, y=b_mat_viable[,k], col="grey80")}
    
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
    rk_v[[i]] <- rk_viable
    b_mats_v[[i]] <- b_mat_viable
    bbmsy_mats_v[[i]] <- bbmsy_mat_viable
    er_mats_v[[i]] <- er_mat_viable
    
  }
  
  # Measure correlation
  ###################################################
  
  # Build index (start and end points) for each stock
  # This is used to index the correlation matrix below
  rk_v_n_all <- sapply(rk_v, nrow) # number of viable r/k pairs for each stock
  finishes <- cumsum(rk_v_n_all)
  starts <- c(0,finishes[1:(length(finishes)-1)])+1
  indices <- cbind(starts, finishes)
  
  # Merge viable effort time series for all stocks
  er_v_all <- t(do.call("cbind", er_mats_v)) # transposed so that G=F*Ft works
  
  # Calculate covariance matrix then correlation matrix
  cov_mat <- er_v_all %*% t(er_v_all)
  corr_mat <- cov2cor(cov_mat)
  
  # # Overwrite self-comparisons
  # dim(corr_mat)
  # corr_mat1 <- corr_mat
  # for(j in 1:nrow(indices)){
  #   pt1 <- indices[j,1]
  #   pt2 <- indices[j,2]
  #   corr_mat1[pt1:pt2, pt1:pt2] <- NA
  # }
  # # Plot correlation matrix
  # corr_mat2 <- t(apply(corr_mat1, 2, rev))
  # image(x=1:ncol(corr_mat2), y=1:ncol(corr_mat2), z=corr_mat2, 
  #       xaxt="n", yaxt="n", xlab="", ylab="")
  # abline(v=finishes, lwd=1.5)
  # abline(h=ncol(corr_mat2)-finishes, lwd=1.5)
  
  
  # Setup containter to hold highest correlation coefficients per row
  # Columns: row, index1, index2, index3, corr12, corr13, corr23, corr_sum
  spp <- 1:nstocks # each species gets a number
  spp_vec <- rep(1:nstocks, rk_v_n_all) # indexes species identity (1, 2, or 3 etc)
  spp_ind <- unlist(sapply(rk_v_n_all, function(x) 1:x)) # indexes index of r/k pair
  hi_corr_mat <- matrix(nrow=nrow(corr_mat), ncol=1+nstocks+nstocks, data=NA, 
                        dimnames=list(NULL, c("row", 
                                              paste0("index", spp), 
                                              paste0("corr", apply(combn(1:3, 2),2, function(x) paste(x, collapse=""))))))
  hi_corr_mat[,"row"] <- 1:nrow(hi_corr_mat)
  
  # Loop through rows of correlation matix and identify highest correlation in each row-block
  for(j in 1:nrow(corr_mat)){
    
    # Which species am I working on?
    spp_curr <- spp_vec[j]
    spp_others <- spp[!(spp %in% spp_curr)]
    
    # Record index of current species
    hi_corr_mat[j, paste0("index", spp_curr)] <- spp_ind[j]
    
    # Loop through other species blocks
    for(k in 1:length(spp_others)){
      spp_other <- spp_others[k]
      spp_other_ind <- indices[spp_other, 1]:indices[spp_other,2]
      corrs <- corr_mat[j,spp_other_ind]
      hi_corr_mat[j,paste0("index", spp_other)] <- which.max(corrs)
      corr_col <- paste0("corr", paste(sort(c(spp_curr, spp_other)), collapse=""))
      hi_corr_mat[j,corr_col] <- max(corrs)
    }
    
  }
  
  # Sum correlations
  corr_cols <- colnames(hi_corr_mat)[grepl("corr", colnames(hi_corr_mat))]
  corr_sums <- apply(hi_corr_mat[,corr_cols ], 1, sum, na.rm=T) 
  hi_corr_mat <- cbind(hi_corr_mat, corr_sum=corr_sums)
  
  # Mark viable r/k pairs that show high correlation
  n_each_index <- melt(hi_corr_mat[,c("index1", "index2", "index3")], 
                       variable.name="stock", value.name="index") %>% 
    select(-Var1) %>% 
    rename(stock=Var2) %>% 
    mutate(stock=gsub("index", "", stock)) %>% 
    group_by(stock, index) %>% 
    summarize(ncorr=n()) %>% 
    ungroup()

  # Add correlation counts to viable r/k pair table
  for(i in 1:length(rk_v)){
    rk_v_do <- as.data.frame(rk_v[[i]])
    rk_v_do1 <- rk_v_do %>% 
      mutate(index=1:nrow(rk_v_do)) %>% 
      left_join(select(filter(n_each_index, stock==i), index, ncorr), by="index") %>% 
      select(index, r, k, ncorr)
    rk_v[[i]] <- rk_v_do1
  }
  
  # Identify top 10% most highly correlated effort time series
  top_p <- 0.10
  top_n <- ceiling(nrow(hi_corr_mat) * top_p)
  top_corr <- as.data.frame(hi_corr_mat) %>% 
    arrange(desc(corr_sum)) %>% 
    slice(1:top_n)

  # Things to return
  out <- list(stocks=stocks,
              yrs=yrs,
              r_priors=r_priors,
              k_priors=k_priors,
              s2_priors=s2_priors,
              r_try=ri, 
              k_try=ki, 
              rk_viable=rk_v, 
              b_viable=b_mats_v, 
              bbmsy_viable=bbmsy_mats_v,
              er_viable=er_mats_v,
              top_corr)
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
    rs <- out$r_try[,i]
    ks <- out$k_try[,i]
    rk_viable <- out$rk_viable[[i]]
    b_viable <- out$b_viable[[i]]
    bbmsy_viable <- out$bbmsy_viable[[i]]
    er_viable <- out$er_viable[[i]]
    r_priors <- out$r_priors
    k_priors <- out$k_priors
    
    # Plot r/k pairs
    #########################################
    
    # Plot r/k pairs
    plot(ks ~ rs, log="xy", bty="n", las=1, pch=15, col="gray80", xpd=NA,
         xlim=r_priors[i,], ylim=k_priors[i,], xlab="r", ylab="k", main=spp[i])
    
    # Add viable r/k pairs
    points(rk_viable$r, rk_viable$k, pch=15, col="grey30")
    
    # # Add most highly correlated pairs
    # rk_v_ind <- unlist(top_corr[,paste0("index", i)])
    # rk_corr <- subset(rk_viable, index %in% rk_v_ind)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="darkorange")
    
    # # Add most common highly correlated pair
    # rk_mode <- mode(rk_v_ind)
    # rk_corr <- subset(rk_viable, index==rk_mode)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="green")
    
    # Add repeatedly highly correlated r/k pairs
    rk_corr <- subset(rk_viable, ncorr>=5)
    points(x=rk_corr$r, y=rk_corr$k, pch=15, col="darkorange")
  
    # Add true value
    points(x=true$r[i], y=true$k[i], pch=15, col="red", cex=1.3)
    
    # Plot BBMSY trajectories
    #########################################
    
    # Plot BBMSY trajectories
    plot(bbmsy_viable[,1] ~ yrs, type="n", bty="n", las=1,
         ylim=c(0, 2), xlab="Year", ylab=expression("B/B"["MSY"]))
    for(k in 1:ncol(bbmsy_viable)){lines(x=yrs, y=bbmsy_viable[,k], col="grey80")}
    lines(x=yrs, y=true$bbmsy_ts[,i+1], lwd=1.5, col="red")
    lines(x=c(0, max(yrs)), y=c(0.5, 0.5), lty=3)
    
    # Plot biomass trajectories
    #########################################
    
    # Plot biomass trajectories
    plot(b_viable[,1] ~ yrs, type="n", bty="n", las=1,
         ylim=c(0, max(b_viable)), xlab="Year", ylab="Biomass")
    for(k in 1:ncol(b_viable)){lines(x=yrs, y=b_viable[,k], col="grey80")}
    lines(x=yrs, y=true$b_ts[,i+1], lwd=1.5, col="red")
    
    # Plot exploitation trajectories
    #########################################
    
    # # Plot exploitation trajectories
    # plot(er_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, 1), xlab="Year", ylab="Exploitation rate")
    # for(k in 1:ncol(er_viable)){lines(x=yrs, y=er_viable[,k], col="grey80")}
    # lines(x=yrs, y=true$er_ts[,i+1], lwd=1.5, col="red")
    
  }
  
  
}
