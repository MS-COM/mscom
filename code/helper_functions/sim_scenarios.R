## Function to run simulation scenarios
## requires data frame with specific scenario settings, name of fishery (hard-coded parameter values), and seed for iteration

sim_scenarios <- function(savedir, scen_df, fishery, seed, iter){


	byScen <- lapply(1:nrow(scen_df), function(x){
	# for(x in 1:nrow(scen_df)){

		if(fishery=="pelagic_longline"){

			species <- c("tuna","billfish","shark")		
		
			## maximum F options
			if(scen_df[x,"MaxF"]=="past_fmsy") maxF <- rep(1.2,length(species))
			if(scen_df[x,"MaxF"]=="at_fmsy") maxF <- rep(1, length(species))
			if(scen_df[x,"MaxF"]=="below_fmsy") maxF <- rep(0.8, length(species))		
			
			## process error in population biomass
			if(scen_df[x,"PopDyn"]=="deterministic"){
				sigR <- 0
				rho <- 0
			}
			if(scen_df[x,"PopDyn"]=="variable"){
				sigR <- 0.15
				rho <- 0.4
			}	

			## start data collection when population is unfished or fished (initial depletion)
			if(scen_df[x,'DataStart']=="unfished") ystart <- rep(1,length(species))
			if(scen_df[x,"DataStart"]=="fished") ystart <- rep(21,length(species))
		
			input_df <- data.frame("SpeciesName"=species,
								"EffDyn"=scen_df[x,"EffDyn"],
								"YearStart"=ystart,
								"MSYscalar"=maxF,
								"SigmaR"=sigR, "rho"=rho,
								"SigmaF"=scen_df[x,"EffSD"],
								"r" = c(0.6, 0.67, 0.11),
								"K"=c(304, 227, 738),
								"p"=rep(0.2,length(species)))	

		}

		## simulate populations
		sim <- sim_pops(input_df=input_df, 
						nyears=scen_df[x,"Nyears"], 
						seed=seed, 
						model="biomass-dynamic")		
		return(sim)

	})
	save(byScen, file=file.path(savedir, paste0("sim", iter, "_", fishery, "_byScenario.Rdata")))
	return(byScen)
}