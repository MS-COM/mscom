
# Clean
rm(list = ls())

# Setup simulations
################################################################################

# Directories
datadir <- "data/simulated/data"
plotdir <- "data/simulated/figures"

# Load simulations
load(file.path(datadir, "simulated_multispecies_fisheries.Rdata"))

# Setup simulations
################################################################################

# Scenarios to plot
toplot <- scenarios %>% 
  filter(iter==1 & nstocks==10)

# Setup figure
figname <- "simulation_plots.pdf"
pdf(file.path(plotdir, figname), width=8.5, height=11)
par(mfrow=c(10,3), mar=c(2,4,0.5,0.5), oma=c(2,2,4,1))

# Loop through scenarios to plot
i <- 1
# for(i in 1:2){
for(i in 1:nrow(toplot)){
  
  # Subset data
  scen <- toplot$scenario_id[i]
  sdata <- filter(sims, scenario==scen)
  
  # Build scenario label
  scen_details <- scenarios[scenarios$scenario_id==scen,]
  scen_label <- paste(scen_details$id, ", ", scen_details$ed, ", ", scen_details$sigmaP, " sigmaP, ",
                      scen_details$sigmaF, " sigmaF, ", scen_details$nyrs)
  
  # Max year
  yr2 <- max(sdata$year)
  
  # Loop through species
  species <- sort(unique(sdata$stock))
  for(j in 1:length(species)){
    
    # Subset data
    spp <- species[j]
    spdata <- filter(sdata, stock==spp)
    iters <- sort(unique(spdata$iter))
    
    # Plot effort
    ymax <- freeR::ceiling1(max(spdata$er), 0.2)
    plot(er ~ year, spdata, type="n", bty="n", las=1,
         xlab="", ylab="", xlim=c(0, yr2), ylim=c(0,ymax))
    for(k in 1:length(iters)){
      idata <- filter(spdata, iter==iters[k])
      lines(idata$year, idata$er)
    }
    text(x=0, y=ymax, pos=4, xpd=NA, label=spp)
    
    # Plot catch
    ymax <- freeR::ceiling1(max(spdata$catch), 10)
    plot(catch ~ year, spdata, type="n", bty="n", las=1, 
         xlab="", ylab="", xlim=c(0, yr2), ylim=c(0,ymax))
    for(k in 1:length(iters)){
      idata <- filter(spdata, iter==iters[k])
      lines(idata$year, idata$catch)
    }
    
    # Plot status
    ymax <- freeR::ceiling1(max(spdata$bbmsy), 0.5)
    plot(bbmsy ~ year, spdata, type="n", bty="n", las=1, 
         xlab="", ylab="", xlim=c(0, yr2), ylim=c(0,ymax))
    for(k in 1:length(iters)){
      idata <- filter(spdata, iter==iters[k])
      lines(idata$year, idata$bbmsy)
    }
    lines(x=c(0,yr2), y=c(0.5, 0.5), lty=3)
    
  }
  
  # Add axis labels
  mtext("Exploitation rate", outer=T, side=2, line=-1, adj=0.5, cex=0.8)
  mtext("Catch (1000s mt)", outer=T, side=2, line=-22, adj=0.5, cex=0.8)
  mtext(expression("B/B"["MSY"]), outer=T, side=2, line=-42, adj=0.5, cex=0.8)
  mtext(scen_label, outer=T, side=3, adj=0, line=1, font=2, cex=0.9)
  
  
}


# Off
dev.off()




