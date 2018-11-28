
################################################################################
# MS-cMSY priors
################################################################################

# Read r priors by family
r_vals_fam <- read.csv("data/priors/data/r_priors_by_family.csv", as.is=T)

# Source FishLife r prior code
source("data/priors/fishlife2/fishlife_r_priors.R")

# R priors
calc_r_priors <- function(key){
  
  # Get FishLife r priors
  fl_r_list <- lapply(key$species, function(x) r_prior(x))
  fl_r_df <- do.call("rbind", fl_r_list)
  fl_r_df_format <- fl_r_df %>% 
    select(species, ln_r_mu, ln_r_sd) %>% 
    mutate(species=as.character(species))
  
  # Build r prior table
  r_priors <- key %>% 
    select(stock, species, family, resilience) %>% 
    # Add FishLife r priors
    left_join(fl_r_df_format, by="species") %>% 
    # Add Neubauer family-level r priors
    left_join(r_vals_fam, by="family") %>% 
    # Set resilience priors and record final prior source
    mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")),
           r_lo=c(0.01, 0.01, 0.01, 0.50)[resilience],
           r_hi=c(0.15, 0.50, 1.00, 1.25)[resilience],
           use=ifelse(!is.na(ln_r_mu), "FishLife", ifelse(!is.na(r_mu), "family", "resilience")))
  
  # Print message if resilience priors are used
  if(any(r_priors$use=="resilience")){
    stocks_with_res_r <- r_priors$stock[r_priors$use=="resilience"]
    print(paste("r priors for the following stocks are based on resilience instead of family: ", paste( stocks_with_res_r, collapse=", ")))
  }
  return(r_priors)
}

# K priors
calc_k_priors <- function(data){
  k_priors <- data %>% 
    group_by(stock) %>% 
    summarize(cmax=max(catch),
              k_lo=cmax*2,
              k_hi=cmax*25)
  return(k_priors)
}

# Initial saturation priors
# yrs <- 1900:2015
# s1 <- calc_sat1_priors(data=data.frame(year=yrs, stock=yrs))
# plot(x=yrs, y=s1$s1_lo, type="l", ylim=c(0,2)); lines(x=yrs, y=s1$s1_hi)
calc_sat1_priors <- function(data){
  s1_priors <- data %>% 
    group_by(stock) %>% 
    summarize(yr1=min(year),
              s1_lo=ifelse(yr1<=1945, 0.75, NA),
              s1_lo=ifelse(yr1>1945 & yr1<1980, 0.75+(0.2-0.75)/(1980-1945)*(yr1-1945), s1_lo),
              s1_lo=ifelse(yr1>=1980, 0.2, s1_lo),
              s1_hi=1)
  return(s1_priors)
}

# Final saturation priors
calc_sat2_priors <- function(data){
  s2_priors <- data %>% 
    group_by(stock) %>% 
    summarize(cfinal=catch[year==max(year)],
              cmax=max(catch),
              cratio=cfinal/cmax,
              s2_lo=0 + 0.2*cratio,
              s2_hi=0.5 + 0.4*cratio)
  return(s2_priors)
}


################################################################################
# Diagnostic plots
################################################################################

# Plot viable r/k pairs
plot_rk <- function(id_rk_viable){
  par(mfrow=c(1,1))
  rmin <- pmax(0.01, floor(min(id_rk_viable$r)/0.1)*0.1)
  rmax <- ceiling(max(id_rk_viable$r)/0.1)*0.1
  plot(k/1000 ~ r, id_rk_viable, log="xy", bty="n", las=1, pch=15,
       xlim=c(rmin, rmax), ylim=unlist(k_prior[, c("k_lo", "k_hi")])/1000,
       xlab="Intrinsic growth rate, r", ylab="Carrying capacity, K", col="gray80")
}

# Plot viable biomass trajectories
plot_b <- function(b_mat_viable, yrs){
  par(mfrow=c(1,1))
  yr1 <- floor(min(yrs) / 10) * 10
  yr2 <- ceiling(max(yrs) / 10) * 10
  plot(b_mat_viable[,1]/1000 ~ yrs, type="n", bty="n", las=2, xlim=c(yr1, yr2),
       ylim=c(0, max(b_mat_viable/1000, na.rm=T)), xlab="Year", ylab="Biomass (1000s)")
  for(k in 1:ncol(b_mat_viable)){lines(x=yrs, y=b_mat_viable[,k]/1000, col="grey80")}
}

# Plot viable B/BMSY trajectories
plot_bbmsy <- function(bbmsy_mat_viable, yrs){
  par(mfrow=c(1,1))
  yr1 <- floor(min(yrs) / 10) * 10
  yr2 <- ceiling(max(yrs) / 10) * 10
  plot(bbmsy_mat_viable[,1] ~ yrs, type="n", bty="n", las=2, xlim=c(yr1, yr2),
       ylim=c(0, max(bbmsy_mat_viable, na.rm=T)), xlab="Year", ylab=expression("B/B"["MSY"]))
  for(k in 1:ncol(bbmsy_mat_viable)){lines(x=yrs, y=bbmsy_mat_viable[,k], col="grey80")}
  lines(x=c(yr1, yr2), y=c(0.5, 0.5), lty=3)
}

# Plot viable exploitation trajectories
plot_er <- function(er_mat_viable, yrs){
  par(mfrow=c(1,1))
  yr1 <- floor(min(yrs) / 10) * 10
  yr2 <- ceiling(max(yrs) / 10) * 10
  plot(er_mat_viable[,1] ~ yrs, type="n", bty="n", las=2, xlim=c(yr1, yr2),
       ylim=c(0, max(er_mat_viable)), xlab="Year", ylab="Exploitation rate")
  for(k in 1:ncol(er_mat_viable)){lines(x=yrs, y=er_mat_viable[,k], col="grey80")}
}

# Plot viable U/UMSY trajectories
plot_uumsy <- function(uumsy_mat_viable, yrs){
  par(mfrow=c(1,1))
  yr1 <- floor(min(yrs) / 10) * 10
  yr2 <- ceiling(max(yrs) / 10) * 10
  plot(uumsy_mat_viable[,1] ~ yrs, type="n", bty="n", las=2, xlim=c(yr1, yr2),
       ylim=c(0, max(uumsy_mat_viable, na.rm=T)), xlab="Year", ylab=expression("U/U"["MSY"]))
  for(k in 1:ncol(uumsy_mat_viable)){lines(x=yrs, y=uumsy_mat_viable[,k], col="grey80")}
  lines(x=c(yr1, yr2), y=c(1, 1), lty=3)
}

################################################################################
# Helper functions
################################################################################

# Look up families
get_family <- function(species){
  
  # FB/SLB taxa key
  taxa_key_fb <- rfishbase::load_taxa(server="https://fishbase.ropensci.org") %>% mutate(type="fish") %>% select(type, everything())
  taxa_key_slb <- rfishbase::sealifebase %>% mutate(type="invert") %>% select(type, everything())
  taxa_key <-  taxa_key_fb %>%
    bind_rows(taxa_key_slb) %>%
    setNames(tolower(names(.))) %>%
    mutate(sciname=paste(genus, species)) %>%
    select(type, class, order, family, genus, species, sciname) %>%
    unique()
  
  # Build family key
  spp <- unique(species)
  family_key <- taxa_key %>%
    filter(sciname %in% spp) %>% 
    select(sciname, family) %>% 
    rename(species=sciname)
  return(family_key)
  
}


################################################################################
# Fit MS-cMSY
################################################################################

# Multi-species catch-only model
# data <- fdata; key <- fkey; npairs <- 1000
ms_cmsy <- function(data, key, npairs=5000){
  
  # Alphabetize
  key <- arrange(key, stock)
  data <- arrange(data, stock, year)
  
  # Get info
  stocks <- key$stock
  nstocks <- length(stocks)
  
  # Check data inputs
  ###################################################
  
  # Data must contain the following columns: stock, year, catch
  if( any(!c("stock", "year", "catch")%in%colnames(data))){
    stop("The 'data' object must include the following column names: stock, year, catch")}
  
  # Key must contain the following columns: stock, family, resilience, id_fixed
  if( any(!c("stock", "family", "resilience", "id_fixed")%in%colnames(key))){
    stop("The 'key' object must include the following column names: stock, family, resilience, id_fixed")}
  
  # Check resilience values
  if(any(!key$resilience%in%c("Very low", "Low", "Medium", "High", NA))){
    stop("Resilience values must be one of the following: Very low, Low, Medium, High, NA")}
  
  # Catch time series must be continuous (i.e., no NAs between 1st and last year)
  cstats <- data %>% 
    group_by(stock) %>% 
    summarize(yr1=min(year), 
              yr2=max(year),
              nyr=n(),
              nyr_check=length(yr1:yr2),
              pass=nyr==nyr_check)
  fail <- cstats$stock[cstats$pass==F]
  if(length(fail)>0){
    stop(paste("The following stocks have non-continuous catch time series (i.e, NAs in catch data):", paste(fail, collapse=", ")))
  }
  
  # Must be catch data for stocks in the key
  key_stocks <- key$stock
  data_stocks <- sort(unique(data$stock))
  if(any(!key_stocks%in%data_stocks)){stop("There are stocks in the 'key' without catch data.")}
  if(any(!data_stocks%in%key_stocks)){stop("There are stocks with catch data but missing from the 'key'.")}
  
  # Calculate priors
  ###################################################
  
  # Calculate r and k priors
  r_priors <- calc_r_priors(key)
  k_priors <- calc_k_priors(data)
  
  # Calculate saturation priors
  s1_priors <- calc_sat1_priors(data)
  s2_priors <- calc_sat2_priors(data)
  
  # Simulate biomass trajectories
  ###################################################
  
  # Lists to hold viable combos and associated biomass/exploitation trajectories
  id_rk_v <- list()
  b_mats_v <- list()
  bbmsy_mats_v <- list()
  er_mats_v <- list()
  uumsy_mats_v <- list()
  
  # Loop through stocks to identify viable r/k pairs and trajectoris
  for(i in 1:nstocks){
    
    # Get info
    stock_i <- stocks[i]
    sdata <- filter(data, stock==stock_i)
    sinfo <- filter(key, stock==stock_i)
    yrs <- sdata$year
    nyrs <- length(yrs)
    c_vec <- sdata$catch
    
    # Get priors
    r_prior <- filter(r_priors, stock==stock_i)
    k_prior <- filter(k_priors, stock==stock_i)
    s1_prior <- filter(s1_priors, stock==stock_i)
    s2_prior <- filter(s2_priors, stock==stock_i)
    
    # Get initial depletions (IDs) to evaluate
    id_fixed <- sinfo$id_fixed
    if(id_fixed==T){
      ids <- 1
    }else{
      ids <- seq(s1_prior$s1_lo, s1_prior$s1_hi, length.out=5)
    }
    
    # Get r/k pairs to evaluate
    r_prior_method <- r_prior$use
    # FishLife-based r priors
    if(r_prior_method=="FishLife"){
      rs <- rlnorm(npairs, meanlog=r_prior$ln_r_mu, sdlog=r_prior$ln_r_sd)
    }
    # Family-based r priors
    if(r_prior_method=="family"){
      rs <- rlnorm(npairs, meanlog=log(r_prior$r_mu), sdlog=r_prior$r_sd)
    }
    # Resilience-based r priors
    if(r_prior_method=="resilience"){
      rs <- exp(runif(npairs, log(r_prior$r_lo), log(r_prior$r_hi)))
    }
    ks <- exp(runif(npairs, log(k_prior$k_lo), log(k_prior$k_hi)))
    rk_pairs <- cbind(r=rs, k=ks, viable=rep(NA, npairs))
    
    # Build ID/r/k combos to evaluate
    id_rk_combos <- as.data.frame(do.call("rbind",
                       lapply(ids, function(x) cbind(id=rep(x, nrow(rk_pairs)), rk_pairs))))
    
    # Loop through ID/r/k combos to see if they produce viable biomass trajectories
    p <- 0.2 # Pella-Tomlinson shape parameter (max growth @ 40% of carrying capacity)
    sigmaP <- 0.1 # process error on productivity term
    b_mat <- matrix(data=NA, nrow=nyrs, ncol=nrow(id_rk_combos), dimnames=list(yrs, 1:nrow(id_rk_combos)))
    for(j in 1:nrow(id_rk_combos)){
      id <- id_rk_combos$id[j]
      r <- id_rk_combos$r[j]
      k <- id_rk_combos$k[j]
      b_mat[1,j] <- k * id
      for(yr in 2:nyrs){
        # b_mat[yr,j] <- b_mat[yr-1,j] +  r*b_mat[yr-1,j]/p*(1-(b_mat[yr-1,j]/k)^p) - c_vec[yr-1] # without process error
        b_mat[yr,j] <- b_mat[yr-1,j] +  r*b_mat[yr-1,j]/p*(1-(b_mat[yr-1,j]/k)^p)*exp(rnorm(1,0,sigmaP)) - c_vec[yr-1] # with process error
      }
    }
    
    # Create saturation matrix
    s_mat <- t(t(b_mat)/ks)
    
    # Reduce to viable r/k pairs and trajectories
    check0 <- apply(b_mat, 2, function(x) sum(x<0, na.rm=T)==0) # check B doesn't go below 0
    s2_lo <- s2_prior$s2_lo
    s2_hi <- s2_prior$s2_hi
    s2_vec <- s_mat[nrow(s_mat),]
    checkS <- s2_vec > s2_lo & s2_vec < s2_hi & !is.na(s2_vec) # check final yr saturation inside prior
    viable <- check0 & checkS # merge positive biomass and final saturation checks
    id_rk_combos[,"viable"] <- viable
    nviable <- sum(viable)
    b_mat_viable <- b_mat[,viable]
    id_rk_viable <- id_rk_combos[viable,]
    if(nviable==0){print(paste("No viable trajectories found for stock:", stock_i))}
    
    # Derive B/BMSY
    bmsy <- id_rk_viable[,"k"] * (1 / (p+1))^(1/p)
    bbmsy_mat_viable <- t(t(b_mat_viable) / bmsy)
    
    # Calculate exploitation rate
    # I name the rows and columns so that I can validate the covariance matrix below
    er_mat_viable <- c_vec / b_mat_viable
    rownames(er_mat_viable) <- yrs
    colnames(er_mat_viable) <- paste0(LETTERS[i], 1:ncol(er_mat_viable))
    
    # Derive U/UMSY
    umsy <- id_rk_viable[,"r"] / p * (1-1/(p+1))
    uumsy_mat_viable <- t(t(er_mat_viable) / umsy)
    
    # Diagnostic plots
    # plot_rk(id_rk_viable)
    # plot_b(b_mat_viable, yrs)
    # plot_bbmsy(bbmsy_mat_viable, yrs)
    # plot_er(er_mat_viable, yrs)
    # plot_uumsy(uumsy_mat_viable, yrs)
    
    # Record results
    id_rk_v[[i]] <- id_rk_viable
    b_mats_v[[i]] <- b_mat_viable
    bbmsy_mats_v[[i]] <- bbmsy_mat_viable
    er_mats_v[[i]] <- er_mat_viable
    uumsy_mats_v[[i]] <- uumsy_mat_viable
    
  }
  
  # Measure correlation
  ###################################################
  
  # Build index (start and end points) for each stock
  # This is used to index the correlation matrix below
  id_rk_v_n_all <- sapply(id_rk_v, nrow) # number of viable r/k pairs for each stock
  finishes <- cumsum(id_rk_v_n_all)
  starts <- c(0,finishes[1:(length(finishes)-1)])+1
  indices <- cbind(starts, finishes)
  
  # Reduce viable effort time series to period of time overlap
  yrs_overlap <- na.omit(dcast(data, year ~ stock, value.var="catch"))$year
  er_mats_v_overlap <- sapply(er_mats_v, function(x) x[rownames(x)%in%yrs_overlap,])
  
  # Merge viable effort time series for all stocks
  er_v_all <- t(do.call("cbind", er_mats_v_overlap)) # transposed so that G=F*Ft works
  
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
  spp_vec <- rep(1:nstocks, id_rk_v_n_all) # indexes species identity (1, 2, or 3 etc)
  spp_ind <- unlist(sapply(id_rk_v_n_all, function(x) 1:x)) # indexes index of r/k pair
  # I create the column names first to use as a method for setting up the number of columns
  hi_corr_mat_cols <- c("row", 
                        paste0("index", spp), 
                        paste0("corr", apply(combn(1:nstocks, 2),2, function(x) paste(x, collapse=""))))
  hi_corr_mat <- matrix(nrow=nrow(corr_mat), ncol=length(hi_corr_mat_cols), data=NA, 
                        dimnames=list(NULL, hi_corr_mat_cols))
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
  
  # Average correlations if >2 stocks
  # If there are 2 stocks, there is only 1 correlation
  corr_cols <- colnames(hi_corr_mat)[grepl("corr", colnames(hi_corr_mat))]
  if(nstocks==2){
    corr_avgs <- hi_corr_mat[,corr_cols ]
  }else{
    corr_avgs <- apply(hi_corr_mat[,corr_cols ], 1, mean, na.rm=T) 
  }
  hi_corr_mat <- cbind(hi_corr_mat, corr_avg=corr_avgs)
  
  # Mark viable r/k pairs that show high correlation
  # n_each_index <- melt(hi_corr_mat[,paste0("index", 1:nstocks)], 
  #                      variable.name="stock", value.name="index") %>% 
  #   select(-Var1) %>% 
  #   rename(stock=Var2) %>% 
  #   mutate(stock=gsub("index", "", stock)) %>% 
  #   group_by(stock, index) %>% 
  #   summarize(ncorr=n()) %>% 
  #   ungroup()

  # Add correlation counts to viable r/k pair table
  for(i in 1:length(id_rk_v)){
    id_rk_v_do <- as.data.frame(id_rk_v[[i]])
    id_rk_v_do1 <- id_rk_v_do %>%
      mutate(index=1:nrow(id_rk_v_do)) %>%
      select(index, id, r, k)
    id_rk_v[[i]] <- id_rk_v_do1
  }
  
  # Identify top-XX% most highly correlated effort time series
  top_p <- 0.10
  top_n <- ceiling(nrow(hi_corr_mat) * top_p)
  top_corr <- as.data.frame(hi_corr_mat) %>%
    arrange(desc(corr_avg)) %>%
    slice(1:top_n)
  
  # Identify combos producing > 0.4 mean correlation
  # corr_thresh <- 0.8
  # top_corr <- as.data.frame(hi_corr_mat) %>% 
  #   arrange(desc(corr_avg)) %>% 
  #   filter(corr_avg>=corr_thresh)
  # if(nrow(top_corr)==0){print("No highly correlated effort time series found.")}
  
  # Get effort-constrained trajectories
  # and calculate/assemble final predictions
  er_mats_vv <- list()
  bbmsy_mats_vv <- list()
  uumsy_mats_vv <- list()
  for(i in 1:nstocks){
    # Index of effort-constrained combos
    vv_index <- unlist(top_corr[,paste0("index", i)])
    # ER time series of effort-constrained combos
    er_mat_v <- er_mats_v[[i]]
    er_mat_vv <- er_mat_v[,vv_index]
    er_mats_vv[[i]] <- er_mat_vv
    # B/BMSY time series of effort-constrained combos
    bbmsy_mat_v <- bbmsy_mats_v[[i]]
    bbmsy_mat_vv <- bbmsy_mat_v[,vv_index]
    bbmsy_mats_vv[[i]] <- bbmsy_mat_vv
    # U/UMSY time series of effort constrained combos
    uumsy_mat_v <- uumsy_mats_v[[i]]
    uumsy_mat_vv <- uumsy_mat_v[,vv_index]
    uumsy_mats_vv[[i]] <- uumsy_mat_vv
    # Calculate final predictions
    bbmsy_cmsy <- t(apply(bbmsy_mat_v, 1, function(x) quantile(x, probs=c(0.5, 0.025, 0.975))))
    bbmsy_mscmsy <- t(apply(bbmsy_mat_vv, 1, function(x) quantile(x, probs=c(0.5, 0.025, 0.975))))
    uumsy_mscmsy <- t(apply(uumsy_mat_vv, 1, function(x) quantile(x, probs=c(0.5, 0.025, 0.975))))
    er_mscmsy <- t(apply(er_mat_vv, 1, function(x) quantile(x, probs=c(0.5, 0.025, 0.975))))
    # Assign column names to final predictions
    colnames(bbmsy_cmsy) <- c("bbmsy_cmsy", "bbmsy_lo_cmsy", "bbmsy_hi_cmsy")
    colnames(bbmsy_mscmsy) <- c("bbmsy", "bbmsy_lo", "bbmsy_hi")
    colnames(uumsy_mscmsy) <- c("uumsy", "uumsy_lo", "uumsy_hi")
    colnames(er_mscmsy) <- c("er", "er_lo", "er_hi")
    # Combine final predictions
    f_preds <- data.frame(stock=stocks[i], year=as.numeric(rownames(bbmsy_cmsy)),
                          bbmsy_cmsy, bbmsy_mscmsy, uumsy_mscmsy, er_mscmsy, stringsAsFactors=F)
    if(i==1){preds <- f_preds}else{preds <- rbind(preds, f_preds)}
  }
  
  # Add catch to final predictions
  preds_c <- preds %>% 
    left_join(select(data, stock, year, catch), by=c("stock", "year")) %>%
    select(stock, year, catch, everything())
    
  # Things to return
  out <- list(key=key,
              data=data,
              preds=preds_c,
              # Priors
              r_priors=r_priors, 
              k_priors=k_priors, 
              s1_priors=s1_priors,
              s2_priors=s2_priors,
              # Viable trajectories
              id_rk_v=id_rk_v,
              b_v=b_mats_v, 
              er_v=er_mats_v,
              bbmsy_v=bbmsy_mats_v,
              uumsy_v=uumsy_mats_v,
              # Effort-constrained trajectories
              top_corr=top_corr,
              er_vv=er_mats_vv,
              bbmsy_vv=bbmsy_mats_vv,
              uumsy_vv=uumsy_mats_vv)
  
  # Return
  return(out)
  
}


################################################################################
# Plot MS-cMSY
################################################################################

# Plot MS-cMSY results
plot_ms_cmsy <- function(output, true){
  
  # Extract info
  key <- output$key
  preds <- output$preds
  stocks <- key$stock
  nstocks <- length(stocks)
  
  # Loop through species
  par(mfcol=c(4,nstocks), mar=c(3,2,2,1), oma=c(0,2,0,0), mgp=c(2.5,0.7,0), xpd=NA)
  for(i in 1:nstocks){
    
    # Stock
    stock_i <- stocks[i]
    
    # Year info
    #########################################
    
    # Year info
    sdata <- filter(preds, stock==stock_i)
    yrs <- sdata$year
    yr1 <- floor(min(yrs) / 10) * 10
    yr2 <- ceiling(max(yrs) / 10) * 10
    
    # Plot catch
    #########################################
    
    # Plot catch
    ymax <- max(sdata$catch)
    ylabel <- ifelse(i==1, "Catch", "")
    plot(catch ~ year, sdata, type="l", bty="n", las=2, xaxt="n",
         xlim=c(yr1, yr2), ylim=c(0, ymax), xlab="", ylab=ylabel, main=stock_i)
    axis(side=1, at=seq(yr1, yr2, 10), las=2)
    
    # Plot r/k pairs
    #########################################
    
    # Extract r/k info
    id_rk_v <- output$id_rk_v[[i]]
    id_rk_v$k <-   id_rk_v$k
    top_corr <- output$top_corr
    r_prior <- unlist(select(filter(output$r_priors, stock==stock_i), r_lo, r_hi))
    k_prior <- unlist(select(filter(output$k_priors, stock==stock_i), k_lo, k_hi))
    
    # Plot viable r/k pairs
    # Potentially reduce this to unique r/k pairs
    # There could be redundancy when evaluating multiple IDs
    rmin <- floor(min(id_rk_v$r) / 0.1) * 0.1
    rmax <- ceiling(max(id_rk_v$r) / 0.1) * 0.1
    kmin <- min(id_rk_v$k) 
    kmax <- max(id_rk_v$k)
    ylabel <- ifelse(i==1, "Carrying capacity, K", "")
    plot(k ~ r, id_rk_v, bty="n", las=1, pch=15, col="gray70",
         xlim=c(rmin, rmax), ylim=c(kmin, kmax), 
         xlab="Intrinsic growth rate, r", ylab=ylabel)

    # Add most highly correlated pairs
    id_rk_v_ind <- unlist(top_corr[,paste0("index", i)])
    rk_corr <- subset(id_rk_v, index %in% id_rk_v_ind)
    points(x=rk_corr$r, y=rk_corr$k, pch=15, col=freeR::tcolor("darkorange", 0.6))

    # # Add most common highly correlated pair
    # rk_mode <- mode(rk_v_ind)
    # rk_corr <- subset(rk_viable, index==rk_mode)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="green")

    # # Add repeatedly highly correlated r/k pairs
    # rk_corr <- subset(rk_viable, ncorr>=5)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="darkorange")
    
    # Add truth if available
    if(!missing(true)){
      rk <- filter(true$rk, comm_name==stock_i)
      points(x=rk$r, y=rk$k, col="red", pch=15, cex=1.2)
    }

    # Add legend
    if(i==1){
      legend("topright", bty="n", pch=15, pt.cex=1.3, cex=0.9,
             col=c("grey70", "darkorange"),
             legend=c("Viable", "Effort constrained"))
    }
    
    # Plot B/BMSY trajectories
    #########################################
    
    # Extract B/BMSY trajectories
    bbmsy_v <- output$bbmsy_v[[i]]
    bbmsy_vv <- output$bbmsy_vv[[i]]
    
    # Plot B/BMSY trajectories
    if(!missing(true)){
      ymax <- ceiling(max(bbmsy_v, true$ts$bbmsy, na.rm=T)/0.5) * 0.5
    }else{
      ymax <- ceiling(max(bbmsy_v, na.rm=T)/0.5) * 0.5
    }
    ylabel <- ifelse(i==1, expression("B/B"["MSY"]), "")
    plot(bbmsy_v[,1] ~ yrs, type="n", bty="n", las=2, xaxt="n",
         xlim=c(yr1, yr2), ylim=c(0, ymax), xlab="", ylab=ylabel)
    axis(side=1, at=seq(yr1, yr2, 10), las=2)
    # Viable trajectories
    for(k in 1:ncol(bbmsy_v)){lines(x=yrs, y=bbmsy_v[,k], col="grey70")}
    # Effort constrained trajectories
    for(k in 1:ncol(bbmsy_vv)){lines(x=yrs, y=bbmsy_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
    lines(x=sdata$year, y=sdata$bbmsy, lwd=1.5, col="brown")
    # cMSY trajectory
    lines(x=sdata$year, y=sdata$bbmsy_cmsy, lwd=1.2, col="black")
    # Overfishing line
    lines(x=c(yr1, yr2), y=c(0.5, 0.5), lty=3)
    lines(x=c(yr1, yr2), y=c(1, 1), lty=2)
    
    # Add truth if available
    if(!missing(true)){
      ts <- filter(true$ts, stock==stock_i)
      lines(x=ts$year, y=ts$bbmsy, col="red", lwd=1.3)
    }

    # Plot exploitation trajectories
    #########################################
    
    # Extract ER trajectories
    er_v <- output$er_v[[i]]
    er_vv <- output$er_vv[[i]]

    # Plot exploitation trajectories
    if(!missing(true)){
      ymax <- ceiling(max(er_v, true$ts$er, na.rm=T)/0.5) * 0.5
    }else{
      ymax <- ceiling(max(er_v, na.rm=T)/0.5) * 0.5
    }
    ylabel <- ifelse(i==1, "Exploitation rate", "")
    plot(er_v[,1] ~ yrs, type="n", bty="n", las=2, xaxt="n",
         xlim=c(yr1, yr2), ylim=c(0, ymax), xlab="", ylab=ylabel)
    axis(side=1, at=seq(yr1, yr2, 10), las=2)
    for(k in 1:ncol(er_v)){lines(x=yrs, y=er_v[,k], col="grey80")}
    for(k in 1:ncol(er_vv)){lines(x=yrs, y=er_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
    lines(x=sdata$year, y=sdata$er, lwd=1.5, col="brown")
    
    # Add truth if available
    if(!missing(true)){
      ts <- filter(true$ts, stock==stock_i)
      lines(x=ts$year, y=ts$er, col="red", lwd=1.3)
    }

  }
  
}
