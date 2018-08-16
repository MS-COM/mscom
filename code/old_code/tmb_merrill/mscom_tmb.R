rm(list=ls())

library(TMB)
library(RColorBrewer)

##-------------------------
## Directories
##--------------------------

main_dir <- "code/old_code/tmb_merrill"
src_dir <- file.path(main_dir, "src")

fig_dir <- file.path(main_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))


##--------------------------
## load data file
##--------------------------

simtest <- readRDS(file.path(main_dir, "simtest.rds"))
true <- readRDS(file.path(main_dir, "simtest_true.rds"))
catch <- simtest$catch

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

r_unif <- runif(ncol(catch), 0.05,0.5)
K_unif <- runif(ncol(catch), 50, 1500)

set.seed(444)
r_error <- sapply(1:length(true$r), function(x) min(rlnorm(1, log(true$r[x]), 0.6), 0.5))
K_error <- sapply(1:length(true$K), function(x) rlnorm(1, log(true$K[x]), 0.7))


### -------(1) priors at truth, starting values at truth
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=nrow(catch), 
			"n_s"=ncol(catch), 
			"C_ts"=catch, 
			# "choose_ref"=1,
			"rk_prior"=1, 
			"r_means"=r_true, "r_sds"=rep(10,ncol(catch)), 
			"K_means"=K_true, "K_sds"=rep(1000,ncol(catch)))
Params <- list("logr"=log(r_true), 
				"logK"=log(K_true),
				"logq"=rep(log(1e-2), ncol(catch)), 
				"logsigma"=log(10),
				"lE_t"=rep(log(1), nrow(catch)))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
# Upr <- rep(Inf, length(Obj$par))
# Upr[which(names(Obj$par)=="logr")] <- log(0.99)
# Upr[which(names(Obj$par)=="logK")] <- log(1000)

# Lwr <- rep(-Inf, length(Obj$par))
# Lwr[which(names(Obj$par)=="logr")] <- log(0.001)

Opt <- TMBhelper::Optimize( obj=Obj, loopnum=3 )
# Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr)
Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)
which(abs(df) > 0.001)
 
Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

png(file.path(fig_dir, "CorrectPrior_CorrectStart.png"), height=8, width=10, res=200, units="in")
basic_diagnostics(years=simtest$years, true=true, Report=Report)
dev.off()


### ------- (2) priors slightly off, starting values at truth
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=nrow(catch), 
			"n_s"=ncol(catch), 
			"C_ts"=catch, 
			# "choose_ref"=1,
			"rk_prior"=1, 
			"r_means"=r_error, "r_sds"=rep(10,ncol(catch)), 
			"K_means"=K_error, "K_sds"=rep(1000,ncol(catch)))
Params <- list("logr"=log(r_true), 
				"logK"=log(K_true),
				"logq"=rep(log(1e-2), ncol(catch)), 
				"logsigma"=log(10),
				"lE_t"=rep(log(1), nrow(catch)))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
# Upr <- rep(Inf, length(Obj$par))
# Upr[which(names(Obj$par)=="logr")] <- log(0.99)
# Upr[which(names(Obj$par)=="logK")] <- log(1000)

# Lwr <- rep(-Inf, length(Obj$par))
# Lwr[which(names(Obj$par)=="logr")] <- log(0.001)

Opt <- TMBhelper::Optimize( obj=Obj, loopnum=3 )
# Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr)
Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)
which(abs(df) > 0.001)
 
Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

png(file.path(fig_dir, "IncorrectPrior_CorrectStart.png"), height=8, width=10, res=200, units="in")
basic_diagnostics(years=simtest$years, true=true, Report=Report)
dev.off()

### ------- (3) priors at truth, starting values slightly off
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=nrow(catch), 
			"n_s"=ncol(catch), 
			"C_ts"=catch, 
			# "choose_ref"=1,
			"rk_prior"=1, 
			"r_means"=r_true, "r_sds"=rep(10,ncol(catch)), 
			"K_means"=K_true, "K_sds"=rep(1000,ncol(catch)))
Params <- list("logr"=log(r_error), 
				"logK"=log(K_error),
				"logq"=rep(log(1e-2), ncol(catch)), 
				"logsigma"=log(10),
				"lE_t"=rep(log(1), nrow(catch)))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
# Upr <- rep(Inf, length(Obj$par))
# Upr[which(names(Obj$par)=="logr")] <- log(0.99)
# Upr[which(names(Obj$par)=="logK")] <- log(1000)

# Lwr <- rep(-Inf, length(Obj$par))
# Lwr[which(names(Obj$par)=="logr")] <- log(0.001)

Opt <- TMBhelper::Optimize( obj=Obj, loopnum=3 )
# Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr)
Opt[["final_gradient"]] = Obj$gr( Opt$opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)
which(abs(df) > 0.001)
 
Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

png(file.path(fig_dir, "CorrectPrior_IncorrectStart.png"), height=8, width=10, res=200, units="in")
basic_diagnostics(years=simtest$years, true=true, Report=Report)
dev.off()

### ------- (4) priors slightly off, starting values slightly off
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=nrow(catch), 
			"n_s"=ncol(catch), 
			"C_ts"=catch, 
			# "choose_ref"=1,
			"rk_prior"=1, 
			"r_means"=r_error, "r_sds"=rep(10,ncol(catch)), 
			"K_means"=K_error, "K_sds"=rep(1000,ncol(catch)))
Params <- list("logr"=log(r_error), 
				"logK"=log(K_error),
				"logq"=rep(log(1e-2), ncol(catch)), 
				"logsigma"=log(10),
				"lE_t"=rep(log(1), nrow(catch)))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
# Upr <- rep(Inf, length(Obj$par))
# Upr[which(names(Obj$par)=="logr")] <- log(0.99)
# Upr[which(names(Obj$par)=="logK")] <- log(1000)

# Lwr <- rep(-Inf, length(Obj$par))
# Lwr[which(names(Obj$par)=="logr")] <- log(0.001)

Opt <- TMBhelper::Optimize( obj=Obj, loopnum=3 )
# Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr)
Opt[["final_gradient"]] = Obj$gr( Opt$opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)
which(abs(df) > 0.001)
 
Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

png(file.path(fig_dir, "IncorrectPrior_IncorrectStart.png"), height=8, width=10, res=200, units="in")
basic_diagnostics(years=simtest$years, true=true, Report=Report)
dev.off()

### ------- (4) priors and starting values bad
##--------------------------
## build inputs and object
##--------------------------

Data <- list("n_t"=nrow(catch), 
			"n_s"=ncol(catch), 
			"C_ts"=catch, 
			# "choose_ref"=1,
			"rk_prior"=1, 
			"r_means"=r_unif, "r_sds"=rep(10,ncol(catch)), 
			"K_means"=K_unif, "K_sds"=rep(1000,ncol(catch)))
Params <- list("logr"=log(r_unif), 
				"logK"=log(K_unif),
				"logq"=rep(log(1e-2), ncol(catch)), 
				"logsigma"=log(10),
				"lE_t"=rep(log(1), nrow(catch)))
Map <- list()
    Map[["logsigma"]] <- NA
    Map[["logsigma"]] <- factor(Map[["logsigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
# Upr <- rep(Inf, length(Obj$par))
# Upr[which(names(Obj$par)=="logr")] <- log(0.99)
# Upr[which(names(Obj$par)=="logK")] <- log(1000)

# Lwr <- rep(-Inf, length(Obj$par))
# Lwr[which(names(Obj$par)=="logr")] <- log(0.001)

Opt <- TMBhelper::Optimize( obj=Obj, loopnum=3 )
# Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr)
Opt[["final_gradient"]] = Obj$gr( Opt$opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)
which(abs(df) > 0.001)
 
Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

png(file.path(fig_dir, "RandomPrior_RandomStart.png"), height=8, width=10, res=200, units="in")
basic_diagnostics(years=simtest$years, true=true, Report=Report)
dev.off()

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