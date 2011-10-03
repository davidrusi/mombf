priorp2g <- function(priorp,q,nu=1,prior="iMom") {
  if (!(prior %in% c('normalMom','tMom','iMom'))) stop("Currently only prior=='normalMom', 'tMom' or 'iMom' are implemented")
  if (prior == "normalMom") {
    e <- function(logg) { return((1-2*pmom(-abs(q),tau=exp(logg))-priorp[i])^2) }
    ans <- double(length(priorp))
    for (i in 1:length(priorp)) {
      ans[i] <- exp(nlminb(start=0,objective=e)$par)
    }
  } else if (prior == 'tMom') {
    stop("prior=='tMom' is not currently implemented in priorp2g")
  } else if (prior == "iMom") {
    p <- (1-priorp)/2
    ans <- qgamma(2*p,nu/2,1)*q^2
  }
  return(ans)
}
