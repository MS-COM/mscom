
# Function to fit the MSCOM
# To compile model: compile("mscom.cpp")
# catch <- catch; id_fixed <- T; r_priors=NA; k_priors=NA; id_priors=NA
fit_mscom <- function(catch, id_fixed=NA,
                      r_priors=NA, k_priors=NA, id_priors=NA){
  
  # Load model
  wd <- getwd()
  setwd(tmbdir)
  dyn.load(dynlib("mscom"))
  
  # Useful info
  nyr <- nrow(catch)
  nspp <- ncol(catch)
  cmaxs <- apply(catch, 2, max)
  
  # Unpack priors
  r_prior_yn <- ifelse(length(r_priors)>1, 1, 0)
  k_prior_yn <- ifelse(length(k_priors)>1, 1, 0)
  id_prior_yn <- ifelse(length(id_priors)>1 & id_fixed==F, 1, 0)
  # R priors
  if(r_prior_yn==1){
    r_means <- r_priors[,1]
    r_sds <- r_priors[,2]
  }else{
    r_means <- rep(0.3, nspp) # default start values
    r_sds <- rep(NA, nspp)
  }
  # K priors
  if(k_prior_yn==1){
    k_means <- k_priors[,1]
    k_sds <- k_priors[,2]
  }else{
    k_means <- cmaxs *10 # default start values
    k_sds <- rep(NA, nspp)
  }
  # Depletion priors
  if(id_prior_yn==1){
    id_means <- id_priors[,1] 
    id_sds <- id_priors[,2]
  }else{
    id_means <- rep(1, nspp) # default start values
    id_sds <- rep(NA, nspp)
  }

  # Print statement
  id_text <- ifelse(id_fixed, "Initial depletion fixed at zero.", "Initial depletion estimated.")
  prior_text <- paste0("r=", ifelse(r_prior_yn==1, "prior", "default"), 
                       "; K=", ifelse(k_prior_yn==1, "prior", "default"), 
                       "; depletion=", ifelse(id_fixed==T, "fixed", ifelse(id_prior_yn==1, "prior", "default")))
  cat(id_text, prior_text, sep="\n")
  
  # Package data
  data <- list("n_t"=nrow(catch), # number of years
               "n_s"=ncol(catch), # number of stocks
               "C_ts"=catch, # catch time series [year, stock]
               "r_prior"=r_prior_yn, # use r priors? 0=no, 1=yes
               "r_means"=r_means, # r prior means
               "r_sds"=r_sds, # r prior sds
               "K_prior"=k_prior_yn, # use K priors? 0=no, 1=yes
               "K_means"=k_means, # K prior means
               "K_sds"=k_sds, # K prior sds
               "delta_prior"=id_prior_yn, # use initial depeletion prior? 0=no, 1=yes
               "delta_means"=id_means, # initial depletion prior means
               "delta_sds"=id_sds) # initial depletion prior sds
  
  # Package starting values
  params <- list("logr"=log(r_means), # r, intrinsic growth rate (1 per stock)
                 "logK"=log(k_means), # K, carrying capacity (1 per stock)
                 "delta_s"=rep(1, ncol(catch)), # initial depletion (1 per stock)
                 "logq"=rep(log(1e-2), ncol(catch)), # q, catchability (1 per stock)
                 "logsigmaC"=log(10), # sd of exploitation likelihood (1 overall)
                 "lE_t"=rep(log(1), nrow(catch))) # effort time series (shared btw stocks)
  
  # Fix some parameters
  # Fix logsigmaC to 10 to aid model convergence
  map <- list()
  map[["logsigmaC"]] <- NA
  map[["logsigmaC"]] <- factor(map[["logsigmaC"]])
  # Fix initial depletion to 1 if user says so
  if(id_fixed==T){
    map[["delta_s"]] <- rep(NA, ncol(catch))
    map[["delta_s"]] <- factor(map[["delta_s"]])
  }
  
  # Setup model
  obj <- MakeADFun(data=data, parameters=params, map=map, DLL="mscom", silent=TRUE)
  # check <- TMBhelper::Check_Identifiable(obj)
  
  # If no priors, add upper and lower limits on parameters
  lwr <- rep(-Inf, length(obj$par))
  upr <- rep(Inf, length(obj$par))
  if(r_prior_yn==0){
    lwr[which(names(obj$par)=="logr")] <- log(0.001)
    upr[which(names(obj$par)=="logr")] <- log(0.99)
  }
  if(k_prior_yn==0){
    lwr[which(names(obj$par)=="logK")] <- log(cmaxs)
    upr[which(names(obj$par)=="logK")] <- log(cmaxs*50)
  }
  if(id_prior_yn==0){
    lwr[which(names(obj$par)=="delta_s")] <- 0.25
    upr[which(names(obj$par)=="delta_s")] <- 1.01
  }
  
  # Fit model
  opt <- TMBhelper::Optimize(obj=obj, loopnum=3, lower=lwr, upper=upr)
  
  # Generate reports
  report <- obj$report()
  sdreport <- sdreport(obj)
  
  # Reset working directory
  setwd(wd)
  
  # Return model and reports
  out <- list(opt, report, sdreport, obj)
  return(out)

}
