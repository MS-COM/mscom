rm(list=ls())

library(dplyr)
library(reshape2)
library(ggplot2)
library(TMB)
library(RColorBrewer)

main_dir <- "code/old_code/tmb_merrill"
src_dir <- file.path(main_dir, "src")

fig_dir <- file.path(main_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

funs <- list.files(file.path(main_dir, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(main_dir,"R", funs[x])))

##--------------------------
## Simulate data
##--------------------------
## tuna = yellowfin Thunnus albacares
## billfish = striped marlin Kajikia audax
## shark = mako Isurus oxyrinchus
species <- c("tuna", "billfish", "shark")
oneway_mat <- data.frame("SpeciesName"=species,
						"Fdynamics"="One-way",
						"InitialDepl"=1,
						"PercentFcrash"=c(0.4, 0.6, 0.8),
						"SigmaR"=0, "rho"=0,
						"SigmaF"=0,
						"r" = c(0.6, 0.67, 0.11),
						"K"=c(304, 227, 738))

sim_oneway <- sim_pops(input_mat=oneway_mat,
				nyears=50,
				seed=123, 
				model="biomass-dynamic")

p_oneway <- ggplot(sim_oneway %>% dplyr::filter(variable %in% c("RelativeCatch","ExploitRate","RelativeEffort","Depletion"))) +
	geom_line(aes(x=Year, y=value, colour=Species), lwd=2) +
	facet_wrap(~variable, scales='free_y') +
	theme_lsd() +
	coord_cartesian(ylim=c(0,1.01))
	
## get catch data
catch_oneway <- sim_oneway %>% filter(variable=="Catch")
input_oneway <- sapply(1:length(species), function(x){
	sub <- catch_oneway %>% filter(Species==species[x])
	df <- sub$value
	return(df)
})
colnames(input_oneway) <- species

catch <- input_oneway

##-------------------------
## template file
##--------------------------
setwd(src_dir)
compile("mscom.cpp")

dyn.load( dynlib("mscom") )

##--------------------------
## build inputs and object
##--------------------------
### Round 1 testing: Initial depletion fixed to 1, no r or K prior, likelihood based on catch
Data <- list("n_t"=nrow(catch), 
			"n_s"=ncol(catch), 
			"C_ts"=catch, 
			"rk_prior"=1, 
			"r_means"=rep(0.5,ncol(catch)), "r_sds"=rep(10,ncol(catch)), 
			"K_means"=rep(500,ncol(catch)), "K_sds"=rep(1000,ncol(catch)),
			"delta_prior"=1,
			"delta_means"=rep(1,ncol(catch)), "delta_sds"=rep(1,ncol(catch)))
Params <- list("logr"=rep(log(0.3), ncol(catch)), 
				"logK"=rep(log(500), ncol(catch)),
				"delta_s"=rep(1, ncol(catch)),
				"logq"=rep(log(1e-2), ncol(catch)), 
				"logsigmaC"=log(10),
				"lE_t"=rep(log(1), nrow(catch)))
Map <- list()
    Map[["logsigmaC"]] <- NA
    Map[["logsigmaC"]] <- factor(Map[["logsigmaC"]])

    Map[["delta_s"]] <- rep(NA, ncol(catch))
    Map[["delta_s"]] <- factor(Map[["delta_s"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")
##--------------------------
## optimize
##--------------------------
Upr <- rep(Inf, length(Obj$par))
Upr[which(names(Obj$par)=="logr")] <- log(0.99)
Upr[which(names(Obj$par)=="logK")] <- log(1000)
# Upr[which(names(Obj$par)=="delta_s")] <- 1.01

Lwr <- rep(-Inf, length(Obj$par))
Lwr[which(names(Obj$par)=="logr")] <- log(0.001)
Lwr[which(names(Obj$par)=="logK")] <- log(50)
# Lwr[which(names(Obj$par)=="delta_s")] <- 0.25


Opt <- TMBhelper::Optimize( obj=Obj, loopnum=3, upper=Upr, lower=Lwr )
# Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr)
if("opt" %in% names(Opt)) Opt[["final_gradient"]] = Obj$gr( Opt$opt$par ) 
if("par" %in% names(Opt)) Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)
which(abs(df) > 0.001)
 
Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty

png(file.path(fig_dir, "Round1testing_OnewaySim.png"), height=8, width=10, res=200, units="in")
basic_diagnostics(years=1:nrow(catch), true=sim_oneway, Report=Report, species=species)
dev.off()
