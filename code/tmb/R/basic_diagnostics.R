basic_diagnostics <- function(years, true, Report, species){
	par(mfrow=c(2,2), mar=c(4,4,2,2))
# years <- simtest$years
cols <- brewer.pal(Report$n_s, "Set1")
## biomass
bts <- Report$B_ts
true_find <- true %>% filter(variable=="Biomass")
tbts <- sapply(1:length(species), function(x){
	sub <- true_find %>% filter(Species==species[x])
	df <- sub$value
	return(df)
})
colnames(bts) <- colnames(tbts) <- species
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="Biomass")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}

## exploit rate
bts <- Report$Upred_ts
true_find <- true %>% filter(variable=="ExploitRate")
tbts <- sapply(1:length(species), function(x){
	sub <- true_find %>% filter(Species==species[x])
	df <- sub$value
	return(df)
})
colnames(bts) <- colnames(tbts) <- species
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="Exploitation Rate")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}

## BBmsy
bts <- Report$BBmsy_ts
true_find <- true %>% filter(variable %in% c("Biomass","Bmsy"))
tbts <- sapply(1:length(species), function(x){
	sub <- true_find %>% filter(Species==species[x])
	sub1 <- sub %>% filter(variable=="Biomass")
	sub2 <- sub %>% filter(variable=="Bmsy")
	df <- sub1$value/sub2$value
	return(df)
})
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="B/Bmsy")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}

## UUmsy
bts <- Report$UUmsy_qE_ts
true_find <- true %>% filter(variable %in% c("ExploitRate","Umsy"))
tbts <- sapply(1:length(species), function(x){
	sub <- true_find %>% filter(Species==species[x])
	sub1 <- sub %>% filter(variable=="ExploitRate")
	sub2 <- sub %>% filter(variable=="Umsy")
	df <- sub1$value/sub2$value
	return(df)
})
plot(x=1, y=1, type="n", xlim=c(min(years), max(years)), ylim=c(0, max(c(max(bts),max(tbts)))*1.1), xlab="Year", ylab="U/Umsy")
for(i in 1:ncol(bts)){
	lines(x=years, y=bts[,i], col=cols[i], lwd=2)
	lines(x=years, y=tbts[,i], col=paste0(cols[i],80), lwd=5, lty=3)
}
legend("topleft", legend=c(paste0(species, " estimate"), paste0(species, " true")), col=c(cols[1:Report$n_s], paste0(cols[1:Report$n_s], 80)), lwd=c(rep(3,Report$n_s), rep(5,Report$n_s)), lty=c(rep(1,Report$n_s), rep(3,Report$n_s)))

}