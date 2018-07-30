lh_list <- function(linf, vbk, t0, lwa, lwb, M, M50, S1, S2, h, R0, selex_type){

	Amax <- ceiling(-log(0.01)/M)
	ages <- 0:Amax

	L_a <- linf*(1 - exp(-vbk*(ages - t0)))
	W_a <- lwa * L_a ^ lwb

	Mat_a <- c(0.001, 1 / (1 + exp(M50 - ages[-1])))

	if(selex_type=="logistic") S_a <- (1 / (1 + exp(-log(19) * (ages - S1) / (S2 - S1))))
	if(selex_type=="dome"){
		S_a <- dnorm(ages, mean=8, sd=2)/max(dnorm(ages, mean=8, sd=2))
		S_a[which(S_a > 0.6)] <- 1
		plot(S_a)
		abline(v=S1)
		abline(v=S2)
	}

	out <- NULL
	out$Amax <- Amax
	out$ages <- ages
	out$L_a <- L_a
	out$W_a <- W_a
	out$Mat_a <- Mat_a
	out$M <- M
	out$h <- h
	out$S_a <- S_a
	out$R0 <- R0
	return(out)

}