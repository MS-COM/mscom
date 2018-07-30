#' Per recruit calculation and spawning potential ratio
#'
#' \code{calc_ref} Calculates derived spawning potential ratio: lifetime total egg production in fished:unfished states
#'
#' @author M.B. Rudd
#' @return List, a tagged list of potentially useful benchmarks
#' @details Use this function with uniroot to find the value of F that results in SPR matching the specified reference value (e.g. 0.30 to find F30)
#' @export
calc_ref <- function(K, r, u, nburn, ref=FALSE){
    bt <- ct <- rep(NA, nburn)
    bt[1] <- K
    for(t in 2:nburn){
        ct[t-1] <- bt[t-1] * u
        bt[t] <- bt[t-1] + bt[t-1] * r * (1 - bt[t-1]/K) - ct[t-1]
    }
    depl <- bt[length(bt)]/K

    if(ref==FALSE) return(depl)
            
    ## can use uniroot function on call to calc_ref to calculate the fishing mortality rate that results in a specified SPR, then compare current fishing mortality with this reference point
    if(ref!=FALSE){
        diff <- ref - depl
        return(diff)
    }
}