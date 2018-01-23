rm(list=ls())

library(TMB)
library(RColorBrewer)

##-------------------------
## Directories
##--------------------------

main_dir <- "C:\\merrill\\MS-COM\\mscom\\code\\tmb"
src_dir <- file.path(main_dir, "src")

funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))


##--------------------------
## load data file
##--------------------------

simtest <- readRDS(file.path(main_dir, "simtest.rds"))
true <- readRDS(file.path(main_dir, "simtest_true.rds"))

##-------------------------
## template file
##--------------------------
setwd(src_dir)
compile("mscom.cpp")

dyn.load( dynlib("mscom") )

## ----------------------
## starting values
## ---------------------

r_true <- true$r
K_true <- true$K

r_random <- runif(3, 0.05,0.5)
K_random <- runif(3, 50, 1500)

r_error <- sapply(1:length(true$r), function(x) rlnorm(1, log(true$r[x]), 1))
K_error <- sapply(1:length(true$K), function(x) rlnorm(1, log(true$K[x]), 1))


### -------(1) priors slightly off, starting values at truth, normal likelihood minimizing squared deviation of Ets and Uts
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=length(simtest$years), "n_s"=ncol(simtest$catch), "C_ts"=simtest$catch, "choose_ref"=1, "like_type"=1, "rk_prior"=1, "r_means"=r_error, "r_sds"=rep(10,3), "K_means"=K_error, "K_sds"=rep(1000,3))
Params <- list("logr"=log(r_true), "logK"=log(K_true), "logsigma"=log(10))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
Upr <- rep(Inf, length(Obj$par))
Upr[which(names(Obj$par)=="logr")] <- log(0.99)

# Opt <- TMBhelper::Optimize( obj=Obj, newtonsteps=3 )
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr)
Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)

Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

bts_sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="B_ts"),]
bts_sd[,2][which(is.na(bts_sd[,2]))] <- 0
sd <- read_sdreport(bts_sd, log=FALSE)

basic_diagnostics(years=simtest$years, true=true, Report=Report)



### ------- (2) priors slightly off, starting values slightly off, normal likelihood minimizing squared deviation of Ets and Uts
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=length(simtest$years), "n_s"=ncol(simtest$catch), "C_ts"=simtest$catch, "choose_ref"=1, "like_type"=1, "rk_prior"=1, "r_means"=r_error, "r_sds"=rep(10,3), "K_means"=K_error, "K_sds"=rep(1000,3))
Params <- list("logr"=log(r_error), "logK"=log(K_error), "logsigma"=log(10))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
Upr <- rep(Inf, length(Obj$par))
Upr[which(names(Obj$par)=="logr")] <- log(0.99)

# Opt <- TMBhelper::Optimize( obj=Obj, newtonsteps=3 )
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr)
Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)

Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

bts_sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="B_ts"),]
bts_sd[,2][which(is.na(bts_sd[,2]))] <- 0
sd <- read_sdreport(bts_sd, log=FALSE)

basic_diagnostics(years=simtest$years, true=true, Report=Report)


### ------- (3) priors slightly off, uninformative starting values, normal likelihood minimizing squared deviation of Ets and Uts
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=length(simtest$years), "n_s"=ncol(simtest$catch), "C_ts"=simtest$catch, "choose_ref"=1, "like_type"=1, "rk_prior"=1, "r_means"=r_error, "r_sds"=rep(10,3), "K_means"=K_error, "K_sds"=rep(1000,3))
Params <- list("logr"=log(r_random)), "logK"=log(K_random), "logsigma"=log(10))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
Upr <- rep(Inf, length(Obj$par))
Upr[which(names(Obj$par)=="logr")] <- log(0.99)

# Opt <- TMBhelper::Optimize( obj=Obj, newtonsteps=3 )
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr)
Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)

Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

bts_sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="B_ts"),]
bts_sd[,2][which(is.na(bts_sd[,2]))] <- 0
sd <- read_sdreport(bts_sd, log=FALSE)

basic_diagnostics(years=simtest$years, true=true, Report=Report)



### ------- (4) estimate sigma, priors slightly off, starting values slightly off, normal likelihood minimizing squared deviation of Ets and Uts
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=length(simtest$years), "n_s"=ncol(simtest$catch), "C_ts"=simtest$catch, "choose_ref"=1, "like_type"=1, "rk_prior"=1, "r_means"=r_error, "r_sds"=rep(10,3), "K_means"=K_error, "K_sds"=rep(1000,3))
Params <- list("logr"=log(r_error), "logK"=log(K_error), "logsigma"=log(5))
Map <- list()
    # Map[["logsigma"]] <- NA
    # Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
Upr <- rep(Inf, length(Obj$par))
Upr[which(names(Obj$par)=="logr")] <- log(1)

# Opt <- TMBhelper::Optimize( obj=Obj, newtonsteps=3 )
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr)
Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)

Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

bts_sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="B_ts"),]
bts_sd[,2][which(is.na(bts_sd[,2]))] <- 0
sd <- read_sdreport(bts_sd, log=FALSE)

basic_diagnostics(years=simtest$years, true=true, Report=Report)



### NOTES
## ask Jim about:
### -- lognormal prior TMB
### -- Suggestions for likelihood
## use r meta-analysis for stock status priors project
## test likelihood comparing effort to Eref
## add option to pool many species in simulation
## how correlated does effort have to be?
## how to determine if effort is correlated?
## multispecies with one assessed