basic_diagnostics <- function(years, true, Report){
	par(mfrow=c(2,2))
# years <- simtest$years
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
legend("topleft", legend=c(paste0("species ", 1:Report$n_s, " estimate"), paste0("species ", 1:Report$n_s, " true")), col=c(cols, paste0(cols, 80)), lwd=c(rep(3,Report$n_s), rep(5,Report$n_s)), lty=c(rep(1,Report$n_s), rep(3,Report$n_s)))

}