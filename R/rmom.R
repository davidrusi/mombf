##################################################################################
## Routines to simulate from MOM prior and posterior
##################################################################################


rmvmomPost <- function(y, x, niter=10^3, burnin=round(niter/10), thinning=1, priorCoef, priorVar) {
  p <- ncol(x); n <- length(y)
  if (nrow(x) != n) stop('Dimensions of y and x do not match')
  tau <- as.double(priorCoef@priorPars['tau'])
  r <- as.integer(ifelse(priorCoef@priorDistr=='pMOM', priorCoef@priorPars['r'], 0))
  if (priorVar@priorDistr=='invgamma') {
    a_phi <- as.double(priorVar@priorPars['alpha'])
    b_phi <- as.double(priorVar@priorPars['lambda'])
  } else stop("Only invgamma prior for residual variance is currently implemented")
  ans <- .Call("rmvmomPostCI",as.integer(niter),as.integer(burnin),as.integer(thinning),as.double(y),as.double(x),as.integer(p),as.integer(r),tau,a_phi,b_phi)
  ans <- matrix(ans,ncol=p+1)
  if (is.null(colnames(x))) colnames(ans) <- c(paste('beta',1:ncol(x),sep=''),'phi') else colnames(ans) <- c(colnames(x),'phi')
  return(ans)
}



##################################################################################
## Routines to simulate from truncated Normal
##################################################################################

rtnorm <- function(n, m, sd, lower, upper) {
  if (length(lower) != length(upper)) stop('length of lower and upper must match')
  if (length(lower)>0) {
    lower[1] <- max(-1.0e10,lower[1])
    upper[length(upper)] <- min(1.0e10,upper[length(upper)])
    ans <- .Call("rnorm_truncMultCI",as.integer(n),as.double(lower),as.double(upper),as.double(m),as.double(sd))
  } else {
    ans <- rnorm(n, m, sd)
  }
  return(ans)
}

rtmvnorm <- function(n, m, Sigma, SigmaInv, lower, upper, within, method='Gibbs', burnin=round(.1*n)) {
# Multivariate normal samples under rectangular constraint
# Input
# - n: number of draws
# - m: multivariate normal mean
# - Sigma: multivariate normal covariance
# - lower: vector with lower truncation points
# - upper: vector with upper truncation points
# - within: if TRUE, each variable is truncated to be >=lower and <= upper. If FALSE, it's truncated to be <lower or >upper
# - method: set method=='Gibbs' for Gibbs sampling, and method=='MH' for independent proposal MH
# - burnin: number of burn-in iterations
# Output: n draws obtained via Gibbs sampling after orthogonalization
  if (length(lower)==1) lower <- rep(1,length(m))
  if (length(upper)==1) upper <- rep(1,length(m))
  if (length(lower)!=length(m)) stop('Length of lower and m do not match')
  if (length(upper)!=length(m)) stop('Length of upper and m do not match')
  if (nrow(Sigma)!=length(m) | ncol(Sigma)!=length(m)) stop('Dimensions of m and Sigma do no match')
  if (!(method %in% c('Gibbs','MH'))) stop('Method should be Gibbs or MH')
  method <- as.integer(ifelse(method=='Gibbs',1,2))
  ans <- .Call("rtmvnormCI",as.integer(n), as.double(m), as.double(Sigma), as.double(lower), as.double(upper), as.integer(within), method)
  matrix(ans,ncol=length(m))
}

rtmvnormProd <- function(n, m, Sigma, k=1, lower=0, upper=Inf, burnin=round(.1*n)) {
# Multivariate normal samples under product contraint  lower <= prod(x^k) <= upper
# Input
# - m: multivariate normal mean
# - Sigma: multivariate normal covariance
# - k: power of product
# - lower: lower truncation point
# - upper: upper truncation point
# - burnin: number of burn-in iterations
# Output: n draws obtained via Gibbs sampling
  lowtrunc <- ifelse(lower==0,as.integer(0),as.integer(1))
  uptrunc <- ifelse(upper==Inf,as.integer(0),as.integer(1))
  if (upper==Inf) { uptrunc <- as.integer(0); upper <- as.double(0) } else { uptrunc <- as.integer(1); upper <- as.double(upper) }
  ans= .Call("rtmvnormProdCI",as.integer(n), as.double(m), as.double(Sigma), as.integer(k), as.double(lower), upper, lowtrunc, uptrunc, as.integer(burnin));
  matrix(ans,nrow=n)
}

#rmvnormTrunc0 <- function(n, m, S, Sinv, t0) {
#  if (length(m) != nrow(S)) stop('Dimensions of m and S must match')
#  if (length(t0)==1) t0 <- as.double(rep(t0,length(m)))
#  if (length(t0)!= length(m)) stop('Dimensions of m and t0 must match')
#  if (missing(Sinv)) Sinv <- solve(S)
#  ans <- .Call("rmvnorm_trunc0R", as.integer(n), as.double(m), as.double(Sinv), as.double(t0))
#  matrix(ans,nrow=n)
#}
