####### MARGINAL LIKELIHOOD FOR NORMAL PRIOR

normalidMarginalK <- function(sel, y, x, phi, tau, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under normal prior
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel));
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  phi <- as.double(phi); tau <- as.double(tau)
  logscale <- as.integer(logscale)
  ngroups= p; nvaringroup= as.integer(rep(1,p))
  ans <- .Call("normalidMarginalKI", sel, nsel, n, p, y, sumy2, XtX, ytX, phi, tau, logscale, ngroups, nvaringroup)
  return(ans)
}


normalidMarginalKR <- function(y, x, phi, tau, logscale=TRUE) {
  #Marginal likelihood for normal prior
  n <- length(y); p <- ncol(x)
  if (p==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    S <- t(x) %*% x + diag(p) / tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    ans <- -.5*(sum(y^2) - t(m) %*% S %*% m)/phi - .5*n*log(2*pi*phi) - .5*p*log(tau) - log(sqrt(det(S)))
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}


normalidMarginalU <- function(sel, y, x, alpha=0.001, lambda=0.001, tau=1, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under normal prior (unknown variance)
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) { x <- matrix(x,ncol=1) } else { x <- as.matrix(x) }
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel));
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  tau <- as.double(tau)
  logscale <- as.integer(logscale)
  alpha <- as.double(alpha); lambda <- as.double(lambda)
  ngroups= p; nvaringroup= as.integer(rep(1,p))
  ans <- .Call("normalidMarginalUI",sel,nsel,n,p,y,sumy2,x,XtX,ytX,tau,logscale,alpha,lambda,ngroups,nvaringroup)
  return(ans);
}

normalidMarginalUR <- function(y, x, alpha=0.001, lambda=0.001, tau, logscale=TRUE) {
  #Marginal likelihood for normal prior
  n <- length(y); p <- ncol(x)
  if (ncol(x)==0) {
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    S <- t(x) %*% x + diag(p) / tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    nu <- n + alpha
    ss <- as.numeric(lambda + sum(y^2) - t(m) %*% S %*% m)
    #
    num <- lgamma(nu/2) + .5*alpha*log(lambda/2) + .5*nu*log(2) - .5*nu*log(ss)
    den <- .5*n*log(2*pi) + .5*log(det(S)) + (.5*p)*log(tau) + lgamma(alpha/2)
    ans <- num - den
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
