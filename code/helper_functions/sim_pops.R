#' Operating model
#'
#' \code{sim_pops} Simulate population dynamics for multiple species in the same fishery
#'
#' @author M.B. Rudd, K. Kleisner, et al.
#' @param input_df data frame of species-specific information; columns "SpeciesName", "EffDyn" (one-way or two-way), "InitDepl" (value <= 1), "MSYscalar" multiplier on maximum catch compared to MSY, "SigmaR" recruitment standard deviation, "rho" recruitment autocorrelation, "SigmaF" standard deviation for fishing mortality deviates, "r" intrinsic rate of growth, "K" carrying capacity, 'p' = pella-tomlinson shape parameter
#' @param nyears number of years of time series
#' @param seed seed for random deviates
#' @param model currently only coded for 'biomass-dynamic' but could include 'age-structured' in the future
#' @importFrom FishLife
#' @return named list of attributes of true population/data
#' @export

sim_pops <- function(input_df,
					nyears,
					seed,
					model='biomass-dynamic'){

	tyears <- nyears

#######################
## life history types
#######################
	## **** calc r and K from life history types required for age-structure

	species <- unique(input_df$SpeciesName)
	nspecies <- length(input_df$SpeciesName)

#######################
## Random variables
#######################

	set.seed(seed)

	## Recruitment deviations
	RecDev <- lapply(1:nspecies, function(x){
		sp <- input_df %>% dplyr::filter(species==species[x])
		raw <- with(sp, rnorm(tyears, mean= -(SigmaR ^ 2)/2, sd=SigmaR))

		out <- rep(NA, tyears)
		out[1] <- raw[1]
		for(t in 2:length(raw)){
			out[t] <- with(sp, out[t-1] * rho + sqrt(1 - rho ^ 2) * raw[t])
		}
		return(out)
	})

	FishDev <- lapply(1:nspecies, function(x){
		sp <- input_df %>% dplyr::filter(species==species[x])
		dev <- with(sp, rnorm(tyears, mean= -(SigmaF ^ 2)/2, sd=SigmaF))
		walk <- sapply(1:length(dev), function(x){
			if(x==1) w <- dev[x]
			if(x>1) w <- dev[x-1] + dev[x] 
			return(w)
		})
		return(walk)
	})


#######################################
## initiate fishing mortality dynamics
#######################################

	if(model == "biomass-dynamic"){

		pop <- lapply(1:nspecies, function(x){
			sp <- input_df %>% dplyr::filter(species==species[x])

			r <- sp$r
			K <- sp$K
			p <- sp$p
			msy <- (r*K)/((p+1)^((p+1)/p))
			cmax <- msy * sp$MSYscalar

			if(tolower(sp$EffDyn)=="constant"){
				et <- rep(1, nyears) * exp(FishDev[[x]])
				ct <- et * cmax 
				# q <- umax/max(et)
				# ut <- q*et * exp(FishDev[[x]])
			}
			if(tolower(sp$EffDyn)=="one-way"){
				et <- seq(0.0001,1,length=nyears) * exp(FishDev[[x]])
				ct <- et * cmax 
				# q <- umax/max(et)
				# ut <- q*et * exp(FishDev[[x]])
			}
			if(tolower(sp$EffDyn)=="two-way"){
				et1 <- seq(0.0001,to=1,length=ceiling(nyears/2))
				et2 <- seq(max(et1),to=0.5,length=floor(nyears/2))
				et <- c(et1,et2) * exp(FishDev[[x]])
				ct <- et * cmax 
			}

			bt <- ut <- pt <- rep(NA, tyears)
			bt[1] <- sp$K*sp$InitDepl
			ut[1] <- ct[1]/bt[1]
			pt[1] <- (r/p)*bt[1]*(1-(bt[1]/K)^p)

			for(t in 2:tyears){
				if(ct[t-1] < bt[t-1]) ut[t] <- ct[t-1] / bt[t-1]
				if(ct[t-1] >= bt[t-1]) stop(paste0("population crash for ", species[x],": change MSYscalar"))
				pt[t] <- (r/p)*bt[t-1]*(1-(bt[t-1]/K)^p)
				bt[t] <- bt[t-1] + (pt[t] * exp(RecDev[[x]][t])) - ct[t-1]
			}

			q <- ut[which(et==max(et))]/max(et)

			## find fmsy and bmsy
			bvec <- seq(0.01,1,by=0.01)
			surprod <- sapply(1:length(bvec),function(x) (r/p)*(bvec[x]*K)*(1-(bvec[x])^p))
			bmsy <- bvec[which(surprod==max(surprod))]*K
			umsy <- msy/bmsy


			out <- data.frame("Species"=species[x], "Year"=1:tyears, "Biomass"=bt, "Catch"=ct, "ExploitRate"=ut, "q"=q, "Effort"=et, "Depletion"=bt/sp$K,"Bmsy"=bmsy, "Umsy"=umsy,"MSY"=msy, "UUmsy"=ut/umsy,"BBmsy"=bt/bmsy)
			out2 <- melt(out, id.vars=c("Species","Year"))

			return(out2)
		})

		pop_out <- do.call(rbind, pop)

		cpop <- pop_out %>% dplyr::filter(variable=="Catch")
		cpop2 <- cpop
		cpop2$variable <- "RelativeCatch"
		cpop2$value <- cpop$value/max(cpop$value)

		epop <- pop_out %>% dplyr::filter(variable=="Effort")
		epop2 <- epop
		epop2$variable <- "RelativeEffort"
		epop2$value <- epop$value/max(epop$value)

		pop_out2 <- rbind(pop_out, epop2, cpop2)

		return(pop_out2)

	}
}