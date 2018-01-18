rm(list=ls())

library(TMB)
library(RColorBrewer)

##-------------------------
## Directories
##--------------------------

main_dir <- "C:\\merrill\\MS-COM\\mscom\\code\\tmb"
src_dir <- file.path(main_dir, "src")

##--------------------------
## load data file
##--------------------------

simtest <- readRDS(file.path(main_dir, "simtest.rds"))
true <- readRDS(file.path(main_dir, "simtest_true.rds"))

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


##-------------------------
## template file
##--------------------------
setwd(src_dir)
compile("mscom.cpp")

##--------------------------
## build inputs and object
##--------------------------
dyn.load( dynlib("mscom") )

Data <- list("n_t"=length(simtest$years), "n_s"=ncol(simtest$catch), "C_ts"=simtest$catch, "choose_ref"=1, "rk_prior"=1)
Params <- list("logr"=log(rep(0.3,3)), "logK"=log(rep(300,3)), "e_sigma"=10)
Map <- list()
    Map[["e_sigma"]] <- NA
    Map[["e_sigma"]] <- factor(Map[["e_sigma"]])

Obj <- MakeADFun( data=Data, parameters=Params, map=Map, DLL="mscom")


##--------------------------
## optimize
##--------------------------
Upr = rep(Inf, length(Obj$par))
# Upr[which(names(Obj$par)=="logr")] = log(1)
# Upr[which(names(Obj$par)=="logK")] = log(1000)

Lwr = rep(-Inf, length(Obj$par))
# Lwr[which(names(Obj$par)=="logr")] = log(0.0001)
# Lwr[which(names(Obj$par)=="logK")] = log(50)
# Lwr[which(names(Obj$par)=="e_sigma")] = 0.001

# Opt <- TMBhelper::Optimize( obj=Obj, newtonsteps=3 )
Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lwr, upper=Upr )
Opt[["final_gradient"]] = Obj$gr( Opt$par ) 
df <- Opt[["final_gradient"]]
colnames(df) <- names(Obj$par)


Report <- Obj$report()
Sdreport <- sdreport(Obj)

## priors vs. true values
## parameter estimates vs uncertainty
read_sdreport <- function(InputMat, log=TRUE){
          index <- which(is.na(InputMat[,2])==FALSE)
          if(log==TRUE) return(c( exp(InputMat[index,1]-1.96*InputMat[index,2]), rev(exp(InputMat[index,1]+1.96*InputMat[index,2]))))
          if(log==FALSE) return(c( InputMat[index,1]-1.96*InputMat[index,2], rev(InputMat[index,1]+1.96*InputMat[index,2])))
} 
bts_sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="B_ts"),]
bts_sd[,2][which(is.na(bts_sd[,2]))] <- 0
sd <- read_sdreport(bts_sd, log=FALSE)


par(mfrow=c(2,2))
years <- simtest$years
cols <- brewer.pal(Report$n_s, "Set1")
## biomass
bts <- Report$B_ts
tbts <- true$B_ts
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="Biomass")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}

## exploit rate
bts <- Report$U_ts
tbts <- true$U_ts
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="Exploitation Rate")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}

## BBmsy
bts <- Report$BBmsy_ts
tbts <- sapply(1:ncol(true$B_ts), function(x){
	true$B_ts[,x]/true$Bmsy[x]
})
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="B/Bmsy")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}

## UUmsy
bts <- Report$UUmsy_ts
tbts <- sapply(1:ncol(true$U_ts), function(x){
	true$U_ts[,x]/true$Umsy[x]
})
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="U/Umsy")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}


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