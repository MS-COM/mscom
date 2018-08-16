
# Get FishLife life history
################################################################################

# Packages
library(stringr)
library(FishLife)

# Gather FishLife life history
fishlife <- function(species){
  
  # Setup container
  spp <- sort(unique(species))
  fl <- data.frame(species=spp, linf_cm=NA, k=NA, winf_g=NA, tmax_yr=NA, tmat_yr=NA, 
                   m=NA, lmat_cm=NA, temp_c=NA, stringsAsFactors=F)
  
  # Loop through species
  for(i in 1:nrow(fl)){
    
    # Get spp info
    sciname <- fl$species[i]
    genus <- word(sciname, 1)
    nwords_in_spp <- length(strsplit(sciname, " ")[[1]])
    species <- word(sciname, start=2, end=nwords_in_spp)
    species <- ifelse(species=="spp", "predictive", species)
    
    # Get and plot life history info
    spp_info <- Plot_taxa(Search_species(Genus=genus, Species=species)$match_taxonomy)
    
    # Values are in log-scale except temperature
    spp_lh_vals_log <- spp_info[[1]]$Mean_pred
    spp_lh_vals <- c(exp(spp_lh_vals_log[1:7]), spp_lh_vals_log[8])
    fl[i,2:ncol(fl)] <- spp_lh_vals

  }
  
  # Return
  return(fl)
  
}

# Get FishBase resilience
################################################################################

# Packages
library(rfishbase)

# Get FishBase resilience
# Use "mode" resilience and if multiple categories are tied, use the lowest
resilience <- function(species){
  
  # FB/SLB taxa key
  taxa_key_fb <- load_taxa(server="https://fishbase.ropensci.org") %>% mutate(type="fish") %>% select(type, everything())
  taxa_key_slb <- sealifebase %>% mutate(type="invert") %>% select(type, everything())
  taxa_key <-  taxa_key_fb %>%
    bind_rows(taxa_key_slb) %>% 
    setNames(tolower(names(.))) %>% 
    mutate(sciname=paste(genus, species)) %>% 
    select(type, class, order, family, genus, species, sciname) %>% 
    unique()
  
  # Check whether species are in FishBase
  spp <- species
  spp_do <- spp[spp %in% taxa_key$sciname]
  spp_not_in_fb <- spp[!spp %in% taxa_key$sciname]
  print(paste0("The following species are not in FishBase: ", paste(spp_not_in_fb, collapse=", ")))
  
  # Divide into fish and inverts
  taxa_do <- taxa_key %>% 
    filter(sciname %in% spp_do)
  spp_fin <- taxa_do$sciname[taxa_do$type=="fish"]
  spp_inv <- taxa_do$sciname[taxa_do$type=="invert"]
  
  # Finfish
  #######################################
  
  # If there are finfish
  if(length(spp_fin)>0){
  
    # Get resilience info from Fishbase
    options(FISHBASE_API = "https://fishbase.ropensci.org")
    fin_orig <- stocks(spp_fin)
    lh_fin <- fin_orig %>% 
      select(sciname, Resilience) %>% 
      filter(!is.na(Resilience)) %>% 
      rename(species=sciname, resilience=Resilience) %>% 
      mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")))
    
    # Calculate mode (default is to resolve ties with lower value)
    lh_fin_sum <- as.data.frame.matrix(table(lh_fin$species, lh_fin$resilience))
    lh_fin_sum$resilience <- apply(lh_fin_sum, 1, function(x) colnames(lh_fin_sum)[which.max(x)])
    
    # Format for export
    lh_fin1 <- lh_fin_sum %>% 
      mutate(species=row.names(lh_fin_sum)) %>% 
      select(species, resilience)
    
    # Invertebrates
    #######################################
    
    # If there are inverts
    if(length(spp_fin)>0){
    
      # Get resilience info from Fishbase
      options(FISHBASE_API = "https://fishbase.ropensci.org/sealifebase")
      inv_orig <- stocks(spp_inv)
      lh_inv <- inv_orig %>% 
        select(sciname, Resilience) %>% 
        filter(!is.na(Resilience)) %>% 
        rename(species=sciname, resilience=Resilience) %>% 
        mutate(resilience=factor(resilience, levels=c("Very low", "Low", "Medium", "High")))
      
      # Calculate mode (default is to resolve ties with lower value)
      lh_inv_sum <- as.data.frame.matrix(table(lh_inv$species, lh_inv$resilience))
      lh_inv_sum$resilience <- apply(lh_inv_sum, 1, function(x) colnames(lh_inv_sum)[which.max(x)])
      
      # Format for export
      lh_inv1 <- lh_inv_sum %>% 
        mutate(species=row.names(lh_inv_sum)) %>% 
        select(species, resilience)
      
    }
    
  }
  
  # Merge results
  #######################################
  
  # Merge finfish and inverts
  if(length(spp_fin)>0 & length(spp_inv)>0){ # fish and inverts
    out <- rbind(lh_fin1, lh_inv1) %>% arrange(species)
  }
  if(length(spp_fin)>0 & length(spp_inv)==0){ # only fish
    out <- lh_fin1
  }
  if(length(spp_fin)==0 & length(spp_inv)>0){ # only inverts
    out <- lh_inv1
  }
  
  # Make resilience a factor
  out$resilience <- factor(out$resilience, levels=c("Very low", "Low", "Medium", "High"))
  
  # Return
  return(out)
  
}



# Get invertebrate life history info
options(FISHBASE_API = "https://fishbase.ropensci.org/sealifebase")
inverts <- spp_lh_key$sci_name_fb[spp_lh_key$type=="invertebrate"]
inverts_gen <- spp_lh_key$genus[spp_lh_key$type=="invertebrate"]
inverts_get <- subset(taxa_key, genus %in% inverts_gen)$sciname
species_inv <- species(inverts_get)
popgrowth_inv <- popgrowth(inverts_get)
stocks_inv <- stocks(inverts_get)


