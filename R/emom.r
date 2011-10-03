setGeneric("demom",function(x, tau=1, phi=1, heavyTail= FALSE, logscale=FALSE) standardGeneric("demom"))

setMethod("demom",signature(x='vector'),function(x, tau=1, phi=1, heavyTail= FALSE, logscale=FALSE) {
  V1 <- 1
  pen <- -tau*phi/x^2
  if (!heavyTail) {
    normct <- sqrt(2)
    ans <- pen + dnorm(x,mean=0,sd=sqrt(tau*phi*V1),log=TRUE) + normct
  } else {
    bt <- .5*(x^2/(tau*phi)+1) * tau*phi/x^2
    num <- sqrt(2) + log(2) + .5*log(bt) + log(besselK(sqrt(4*bt),nu=1,expon.scaled=TRUE)) - sqrt(4*bt)
    den <- lgamma(.5) + .5*(log(pi)+log(tau)+log(phi)) + log(1+x^2/(tau*phi))
    ans <- num - den
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
)

setMethod("demom",signature(x='data.frame'),function(x, tau=1, phi=1, heavyTail= FALSE, logscale=FALSE) {
  demom(as.matrix(x),tau=tau,phi=phi,heavyTail=heavyTail,logscale=logscale)
}
)

setMethod("demom",signature(x='matrix'),function(x, tau=1, phi=1, heavyTail= FALSE, logscale=FALSE) {
  require(mvtnorm)
  p <- ncol(x)
  V1 <- diag(p)
  pen <- -tau*phi*rowSums(1/x^2)
  if (!heavyTail) {
    normct <- p*sqrt(2)
    ans <- pen + dmvnorm(x,mean=rep(0,p),sigma=tau*phi*V1,log=TRUE) + normct
  } else {
    stop('heavyTail option not yet implemented for multivariate emom')
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
)


pemom <- function(q, tau = 1, heavyTail = FALSE) integrate(demom,-Inf,q,tau=tau,heavyTail=heavyTail)$value
