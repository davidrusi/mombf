#######################################################################
##
## ASYMMETRIC LAPLACE ROUTINES
##
#######################################################################

dalapl <- function(x, th=0, scale=1, alpha=0, logscale=FALSE) {
  ans <- -abs(th-x)/(sqrt(scale)*ifelse(x<th,1+alpha,1-alpha)) - 0.5*log(scale) - log(2)
  if (!logscale) ans <- exp(ans)
  return(ans)
}

palapl <- function(x, th=0, scale=1, alpha=0) {
    ans <- ifelse(x<=th, exp(-abs(x-th)/(sqrt(scale)*(1+alpha))) * (1+alpha)/2, (1+alpha)/2 + (1 - exp(-abs(x-th)/(sqrt(scale)*(1-alpha)))) * (1-alpha)/2)
    return(ans)
}

ralapl <- function(n, th=0, scale=1, alpha=0) {
    s <- runif(n) < palapl(0,scale=scale,alpha=alpha)
    ans <- double(n)
    ans[s] <- th - rexp(sum(s)) * sqrt(scale)*(1+alpha)
    ans[!s] <- th + rexp(sum(!s)) * sqrt(scale)*(1-alpha)
    return(ans)
}
