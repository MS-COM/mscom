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

#######################
## life history types
#######################
	## **** calc r and K from life history types required for age-structure

	species <- unique(input_df$SpeciesName)
	nspecies <- length(input_df$SpeciesName)

#######################
## Data collection start
#######################

	Tyears <- lapply(1:nspecies, function(x){
		start <- input_df[x,"YearStart"]
		if(start==1) tyears <- nyears
		if(start>1) tyears <- nyears + start
		return(tyears)
	})

#######################
## Random variables
#######################

	set.seed(seed)

	## Process error
	ProcDev <- lapply(1:nspecies, function(x){
		sp <- input_df %>% dplyr::filter(species==species[x])
		raw <- with(sp, rnorm(Tyears[[x]], mean= -(SigmaR ^ 2)/2, sd=SigmaR))

		out <- rep(NA, Tyears[[x]])
		out[1] <- raw[1]
		for(t in 2:length(raw)){
			out[t] <- with(sp, out[t-1] * rho + sqrt(1 - rho ^ 2) * raw[t])
		}
		return(out)
	})

	qDev <- lapply(1:nspecies, function(x){
		sp <- input_df %>% dplyr::filter(species==species[x])
		dev <- with(sp, rnorm(Tyears[[x]], mean= -(SigmaF ^ 2)/2, sd=SigmaF))
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

			## find fmsy and bmsy
			bvec <- seq(0.01,1,by=0.01)
			surprod <- sapply(1:length(bvec),function(x) (r/p)*(bvec[x]*K)*(1-(bvec[x])^p))
			bmsy <- bvec[which(surprod==max(surprod))]*K
			umsy <- msy/bmsy
			umax <- umsy * sp$MSYscalar

			if(tolower(sp$EffDyn)=="constant"){
				et <- rep(1, Tyears[[x]])
				qt <- umax/max(et) * exp(qDev[[x]])
				ut <- qt*et 
			}
			if(tolower(sp$EffDyn)=="one-way"){
				et <- seq(0,1,length=Tyears[[x]])
				qt <- umax/max(et) * exp(qDev[[x]])
				ut <- qt*et
			}
			if(tolower(sp$EffDyn)=="two-way"){
				et1 <- seq(0,to=1,length=ceiling(Tyears[[x]]/2))
				et2 <- seq(max(et1),to=0.5,length=floor(Tyears[[x]]/2))
				et <- c(et1,et2)
				qt <- umax/max(et) * exp(qDev[[x]])
				ut <- qt*et
			}

			bt <- ct <- pt <- rep(NA, Tyears[[x]])
			bt[1] <- sp$K
			ct[1] <- ut[1] * bt[1]
			pt[1] <- (r/p)*bt[1]*(1-(bt[1]/K)^p)

			for(t in 2:Tyears[[x]]){
				ct[t] <- ut[t-1] * bt[t-1]
				pt[t] <- (r/p)*bt[t-1]*(1-(bt[t-1]/K)^p)
				bt[t] <- (bt[t-1] + pt[t] - ct[t-1]) * exp(ProcDev[[x]][t])
			}

			## cut time series to start in fished condition
			btout <- bt[input_df$YearStart[x]:length(bt)]
			ctout <- ct[input_df$YearStart[x]:length(ct)]
			utout <- ut[input_df$YearStart[x]:length(ut)]
			etout <- et[input_df$YearStart[x]:length(et)]
			qtout <- qt[input_df$YearStart[x]:length(qt)]


			out <- data.frame("Species"=species[x], "Year"=seq_along(input_df$YearStart[x]:Tyears[[x]]), "Biomass"=btout, "Catch"=ctout, "ExploitRate"=utout, "qt"=qt, "Effort"=etout, "Depletion"=btout/sp$K,"Bmsy"=bmsy, "Umsy"=umsy,"MSY"=msy, "UUmsy"=utout/umsy,"BBmsy"=btout/bmsy)

			return(out)
		})

		pop_out <- do.call(rbind, pop)

		pop_out <- pop_out %>%
					mutate("RelativeCatch"=Catch/max(Catch)) %>%
					mutate("RelativeEffort"=Effort/max(Effort))

		return(pop_out)

	}
}