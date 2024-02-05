##############################################################################################
##
## ROUTINES TO COMPUTE INTEGRATED LIKELIHOODS
##
##############################################################################################

nlpMarginal <- function(
  sel, y, x, data, smoothterms, nknots=9, groups=1:ncol(x), family="normal",
  priorCoef, priorGroup, priorVar=igprior(alpha=0.01,lambda=0.01),
  priorSkew=momprior(tau=0.348), phi, method='auto', adj.overdisp='intercept', hess='asymp',
  optimMethod, optim_maxit, initpar='none', B=10^5, logscale=TRUE, XtX, ytX
) {
  #Check input
  # format input data
  tmp <- formatInputdata(y=y,x=x,data=data,smoothterms=smoothterms,nknots=nknots,family=family)
  x <- tmp$x; y <- tmp$y; is_formula <- tmp$is_formula
  splineDegree <- tmp$splineDegree
  if (!is.null(tmp$groups)) groups <- tmp$groups
  hasgroups <- tmp$hasgroups
  if (!is.null(tmp$constraints)) constraints <- tmp$constraints
  outcometype <- tmp$outcometype; uncens <- tmp$uncens; ordery <- tmp$ordery
  typeofvar <- tmp$typeofvar
  p= ncol(x); n= length(y)
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  sumy2 <- as.double(sum(y^2)); sumy <- as.double(sum(y))
  colsumsx <- as.double(colSums(x))
  #
  familyint= formatFamily(family, issurvival= length(uncens)>0)$familyint
  if (familyint == 22) { sumlogyfact= as.double(sum(lgamma(y+1))) } else { sumlogyfact= as.double(0) } #Poisson regression
  # check prior and set defaults if necessary
  if (missing(priorCoef)) {
      defaultprior= defaultmom(outcometype=outcometype,family=family)
      priorCoef= defaultprior$priorCoef; priorVar= defaultprior$priorVar
  }
  if (missing(priorGroup)) { if (length(groups)==length(unique(groups))) { priorGroup= priorCoef } else { priorGroup= groupzellnerprior(tau=n) } }

  if (missing(phi)) { knownphi <- as.integer(0); phi <- double(0) } else { knownphi <- as.integer(1); phi <- as.double(phi) }

  # format arguments for .Call
  thinit= getthinit(y=y, x=x, family=family, initpar=initpar, enumerate=TRUE)
  usethinit= thinit$usethinit; thinit= thinit$thinit

  method <- formatmsMethod(method=method, usethinit=usethinit, optimMethod=optimMethod, optim_maxit=optim_maxit, priorCoef=priorCoef, priorGroup=priorGroup, knownphi=0, outcometype=outcometype, family=family, hasgroups=hasgroups, adj.overdisp=adj.overdisp, hess=hess)
  optimMethod <- method$optimMethod; optim_maxit <- method$optim_maxit; adj.overdisp <- method$adj.overdisp; hesstype <- method$hesstype; method <- method$method
  #hesstype <- as.integer(ifelse(hess=='asympDiagAdj',2,1)); optimMethod <- as.integer(ifelse(optimMethod=='CDA',2,1))
    
  B <- as.integer(B)
  tmp= codeGroupsAndConstraints(p=p,groups=groups)
  ngroups= tmp$ngroups; constraints= tmp$constraints; invconstraints= tmp$invconstraints; nvaringroup=tmp$nvaringroup; groups=tmp$groups
  tmp= formatmsPriorsMarg(priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, priorSkew=priorSkew, n=n)
  r= tmp$r; prior= tmp$prior; priorgr= tmp$priorgr; tau=tmp$tau; taugroup=tmp$taugroup; alpha=tmp$alpha; lambda=tmp$lambda; taualpha=tmp$taualpha; fixatanhalpha=tmp$fixatanhalpha

  if (!is_formula) {
    sel <- check_sel_groups(sel, groups)
    sel <- as.integer(sel-1); nsel <- as.integer(length(sel))
  } else {
    if (!missing(sel)) warning("y is of type formula: ignoring sel argument")
    sel <- as.integer(seq(p)-1)
    nsel <- length(sel)
  }

  ans <- nlpMarginalCI(knownphi, sel, nsel, familyint, prior, priorgr, n, p, y, uncens, sumy2, sumy, sumlogyfact, x, colsumsx, XtX, ytX, method, adj.overdisp, hesstype, optimMethod, optim_maxit, thinit, usethinit, B, alpha, lambda, tau, taugroup, taualpha, fixatanhalpha, r, groups, ngroups, nvaringroup, constraints, invconstraints, logscale)
  #ans <- .Call("nlpMarginalCI", knownphi, sel, nsel, familyint, prior, priorgr, n, p, y, uncens, sumy2, sumy, sumlogyfact, x, colsumsx, XtX, ytX, method, adj.overdisp, hesstype, optimMethod, optim_maxit, thinit, usethinit, B, alpha, lambda, tau, taugroup, taualpha, fixatanhalpha, r, groups, ngroups, nvaringroup, constraints, invconstraints, logscale)
  return(ans)
}

check_sel_groups <- function(sel, groups) {
  p <- length(groups); seqp <- seq(p)
  if (is.logical(sel)) sel <- seqp[sel]
  if (any(sel > p)) stop("found index in sel larger than ncol(x). Please make sure all indexes refer to existing variables")
  invsel <- seqp[!(seqp %in% sel)]
  if (any(groups[sel] %in% groups[invsel])) stop("selected indexes incompatible with defined groups. Make sure each group is selected or discarded at once")
  return(sel)
}

##############################################################################################
## INTEGRATED LIKELIHOOD FOR LINEAR MODELS UNDER NORMAL RESIDUALS
##############################################################################################


###
### margpmom.R
###

eprod <- function(m, S, power=1, dof= -1) {
  #Mean of prod (x_i)^(2*power) when x_i ~ T_dof(mu,sigma). Set dof=-1 for N(mu,sigma). Written by John Cook
  ans <- eprod_I(as.double(m),as.double(S), as.integer(length(m)), as.integer(power), as.double(dof))
  ans
}

fmomNeg <- function(th, m, S, phi, tau, r, logscale=TRUE) .5*mahalanobis(th, center=m, cov=S, inverted=TRUE)/phi - r*sum(log(th^2))
fpmomNeg <- function(th, m, S, phi, tau, r) S %*% matrix(th-m, ncol=1)/phi - 2*r/th
fppmomNeg <- function(th, m, S, phi, tau, r) S/phi + 2*r*diag(1/th^2,nrow=length(th))

pmomIntegralApproxR <- function(m, S, phi, tau, r, logscale=TRUE) {
  #Laplace approx to integral N(th; m, phi*solve(S)) prod (th/(phi*tau))^2r wrt th
  opt <- nlminb(m, objective=fmomNeg, gradient=fpmomNeg, m=m, S=S, phi=phi, tau=tau, r=r)$par
  fopt <- -fmomNeg(opt,m=m,S=S,phi=phi,tau=tau,r=r)
  hess <- fppmomNeg(opt,m=m,S=S,phi=phi,tau=tau,r=r)
  ans <- fopt + .5*log(det(S)) - .5*log(det(hess)) - .5*length(m)*log(phi)
  if (!logscale) ans <- exp(ans)
  return(ans)
}

pmomMarginalK <- function(sel, y, x, phi, tau, r=1, method='auto', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product mom prior (known variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - phi: residual variance
# - tau: prior dispersion parameter
# - r: prior power parameter is 2*r
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' for Monte Carlo. 'Plug-in' for plug-in estimate. 'auto' for exact calculation if p<=10, else Laplace approx
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
# - XtX, ytX: optionally, X'X and y'X can be specified to speed up computations
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel));
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  phi <- as.double(phi); tau <- as.double(tau); r <- as.integer(r)
  if (method=='auto') method=-1 else if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else stop("Invalid 'method'")
  method <- as.integer(method)
  B <- as.integer(B); logscale <- as.integer(logscale)
  ngroups= p; nvaringroup= as.integer(rep(1,p))
  ans <- pmomMarginalKI(sel, nsel, n, p, y, sumy2, XtX, ytX, phi, tau, r, method, B, logscale, ngroups, nvaringroup)
  return(ans)
}


pmomMarginalKR <- function(y, x, phi, tau, r=1, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product moment prior (variance phi known)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior proportional to N(th; 0, tau*phi*I) * prod(th^2/(phi*tau))^r
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm)
  n <- length(y); p <- ncol(x)
  if (p==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    ans <- -.5*(sum(y^2) - t(m) %*% S %*% m)/phi - p*(sum(log(seq(1,2*r-1,by=2)))) - .5*n*log(2*pi*phi) - .5*p*log(tau) - log(sqrt(det(S))) - r*p*log(tau*phi)
    if (method=='Laplace') {
      I <- pmomIntegralApproxR(m=m, S=S, phi=phi, tau=tau, r=r, logscale=TRUE)
    } else if (method=='1storder') {
      I <- r*sum(log(m^2))
    } else if (method=='MC') {
      thsim <- rmvnorm(B,m,phi*solve(S))
      eprod <- exp(rowSums(log(thsim^(2*r))))
      I <- log(mean(eprod))
    }
    ans <- ans + I
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}


pmomMarginalU <- function(sel, y, x, alpha=0.001, lambda=0.001, tau=1, r=1, method='auto', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (unknown variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - alpha, lambda: prior for phi is IGamma(alpha/2,lambda/2)
# - tau: prior dispersion parameter
# - r: prior power parameter is 2*r
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' for Monte Carlo. 'Plug-in' for plug-in estimate.
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) { x <- matrix(x,ncol=1) } else { x <- as.matrix(x) }
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel));
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  tau <- as.double(tau); r <- as.integer(r)
  if (method=='auto') method=-1 else if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else stop("Invalid 'method'")
  method <- as.integer(method)
  B <- as.integer(B); logscale <- as.integer(logscale)
  alpha <- as.double(alpha); lambda <- as.double(lambda)
  ngroups= p; nvaringroup= as.integer(rep(1,p))
  ans <- pmomMarginalUI(sel,nsel,n,p,y,sumy2,x,XtX,ytX,tau,r,method,B,logscale,alpha,lambda,ngroups,nvaringroup)
  return(ans);
}

pmomMarginalUR <- function(y, x, r, alpha=0.001, lambda=0.001, tau, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product moment prior (variance phi unknown)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior for th proportional to N(th; 0, tau*phi*I) * prod(th^2/(phi*tau))^r
  # - Prior for phi: IGamma(alpha/2,lambda/2)
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm)
  n <- length(y); p <- ncol(x)
  if (ncol(x)==0) {
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    nu <- 2*r*p + n + alpha
    ss <- as.numeric(lambda + sum(y^2) - t(m) %*% S %*% m)
    V <- S*nu/ss
    #
    if (method=='Laplace') {
      I <- pmomIntegralApproxR(m=m, S=S, phi=nu/(nu-2), tau=tau, r=r, logscale=TRUE)
    } else if (method=='1storder') {
      I <- r*sum(log(m^2))
    } else if (method=='MC') {
      cholV <- t(chol(solve(V)))
      z <- rmvnorm(B,rep(0,p),diag(p))
      thsim <- as.vector(m) + (cholV %*% t(z)) * sqrt(nu/rchisq(B,df=nu))
      eprod <- exp(colSums(log(thsim^(2*r))))
      seprod <- sd(eprod)/sqrt(length(eprod))
      I <- log(mean(eprod))
    } else {
      stop("Only 'Laplace', '1storder' and 'MC' methods are implemented")
    }
    #
    num <- lgamma(nu/2) + .5*alpha*log(lambda/2) + .5*nu*log(2) - .5*nu*log(ss)
    den <- p*(sum(log(seq(1,2*r-1,by=2)))) + .5*n*log(2*pi) + .5*log(det(S)) + (.5*p+r*p)*log(tau) + lgamma(alpha/2)
    ans <- I + num - den
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}


###
### margpimom.R
###

#fimomNeg works both for th of class vector and class matrix
#fpimomNeg, fppimomNeg work only for vector th
fimomNeg <- function(th, XtX, ytX, phi, tau) {
  p <- ncol(XtX)
  if (p==1) th <- matrix(th,ncol=1)
  if (is.vector(th)) {
    ans <- as.vector(.5*(mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE) - 2*ytX %*% matrix(th,ncol=1))/phi + tau*phi*sum(1/th^2) + sum(log(th^2)))
  } else {
    ans <- as.vector(.5*(mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE) - 2*ytX %*% t(th))/phi + tau*phi*rowSums(1/th^2) + rowSums(log(th^2)))
  }
  return(ans)
}

fpimomNeg <- function(th, XtX, ytX, phi, tau) (XtX %*% matrix(th,ncol=1) - t(ytX))/phi - 2*tau*phi/th^3 + 2/th
fppimomNeg <- function(th, XtX, ytX, phi, tau) XtX/phi + diag(6*tau*phi/th^4 - 2/th^2, ncol=length(th))

imomModeK <- function(thini, XtX, ytX, phi, tau) {
  #piMOM mode when phi is known using gradient algorithm
  th <- thini
  err <- 1; niter <- 0
  while ((err > 0.001) & (niter<50)) {
    err <- 0; niter <- niter+1
    for (i in 1:length(th)) {
      a <- c(2*tau*phi, 0, -2, (ytX[i]-sum(XtX[i,-i]*th[-i]))/phi, -XtX[i,i]/phi)
      thnew <- polyroot(a)
      thnew <- Re(thnew[abs(Im(thnew))< 1e-7])
      thnew <- thnew[sign(thnew)==sign(th[i])]
      err <- err+abs(th[i]-thnew)
      th[i] <- thnew
    }
  }
  return(th)
}

imomIntegralApprox <- function(XtX, ytX, phi, tau, logscale=TRUE) {
#Laplace approx to product imom marginal (uses gradient search to find mode)
  m <- as.vector(solve(XtX + tau*diag(nrow(XtX))) %*% t(ytX))
  m <- imomModeK(m, XtX=XtX, ytX=ytX, phi=phi, tau=tau)
  V <- fppimomNeg(m, XtX=XtX, ytX=ytX, phi=phi, tau=tau)
  fopt <- fimomNeg(m,XtX=XtX,ytX=ytX,phi=phi,tau=tau)
  ans <- -fopt - .5*as.numeric(determinant(V,logarithm=TRUE)$modulus)
  if (!logscale) ans <- exp(ans)
  return(list(ans=ans,thopt=m,Vopt=V,objective=fopt))
}


pimomMarginalK <- function(sel, y, x, phi, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (known variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - phi: residual variance
# - tau: prior dispersion parameter
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' to correct Laplace approx via Importance Sampling based on multivariate Cauchy
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
# - XtX, ytX: optionally, X'X and y'X can be specified to speed up computations
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
  if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else stop("Invalid argument 'method'")
  method <- as.integer(method); B <- as.integer(B); logscale <- as.integer(logscale)
  ngroups= p; nvaringroup= as.integer(rep(1,p))
  ans <- pimomMarginalKI(sel, nsel, n, p, y, sumy2, XtX, ytX, phi, tau, method, B, logscale, ngroups, nvaringroup)
  return(ans)
}

pimomMarginalKR <- function(y, x, phi, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (known variance)
# - y: response variable
# - x: design matrix
# - phi: residual variance
# - tau: prior dispersion parameter
# - method: method to approximate the integral. 'Laplace' for Laplace approx. 'MC' to correct Laplace approx via Importance Sampling based on multivariate Cauchy
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
# - XtX, ytX: optionally, X'X and y'X can be specified to speed up computations
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  p1 <- ncol(x); n <- nrow(x)
  if (p1==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    if (n != length(y)) stop("Dimensions of x and y don't match")
    if (missing(XtX)) { XtX <- t(x) %*% x }
    if (missing(ytX)) { ytX <- t(y) %*% x }
    ILaplace <- imomIntegralApprox(XtX=XtX,ytX=ytX,phi=phi,tau=tau,logscale=TRUE)
    k <- .5*p1*log(tau) - .5*sum(y^2)/phi - .5*n*log(2*pi) - .5*(n-p1)*log(phi) - .5*p1*log(pi)
    if (method=='Laplace') {
      ans <- k + ILaplace$ans
    } else if (method=='MC') {
      Vinv <- solve(ILaplace$Vopt)
      uplim <- ILaplace$thopt + 2*sign(ILaplace$thopt)*sqrt(diag(Vinv))
      sdprop <- diag(abs(uplim)/2,ncol=p1)
      Vprop <- sdprop %*% cov2cor(Vinv) %*% sdprop
      #thsim <- rmvnorm(B,rep(0,p1),Vprop)
      #adj <- - fimomNeg(thsim,XtX=XtX,ytX=ytX,phi=phi,tau=tau) - dmvnorm(thsim,rep(0,p1),Vprop,log=TRUE)
      thsim <- rmvt(B,sigma=Vprop,df=1)
      adj <- - fimomNeg(thsim,XtX=XtX,ytX=ytX,phi=phi,tau=tau) - dmvt(thsim,delta=rep(0,p1),sigma=Vprop,df=1,log=TRUE)
      m <- max(adj)
      adj <- log(mean(exp(adj-m+500))) + m - 500
      ans <- k + adj
    } else {
      stop("Only method=='Laplace' and method=='MC' are currently implemented")
    }
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}



pimomMarginalU <- function(sel, y, x, alpha=0.001, lambda=0.001, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (unknown variance)
# - sel: vector with indexes of variables included in the model
# - y: response variable
# - x: design matrix
# - alpha, lambda: prior for phi is IGamma(alpha/2,lambda/2)
# - tau: prior dispersion parameter
# - method: method to approximate the integral for known phi. Integral wrt phi is performed via integrate. 'Laplace' for Laplace approx which may underestimate true value, 'MC' for exact evaluation which can be very computationally expensive. 'Hybrid' combines Laplace for fixed phi with numerical integration wrt phi. The Laplace error is estimated using a single exact evaluation for a value of phi close to the posterior mode
# - B: number of Monte Carlo samples to use (ignored unless method=='MC')
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
  if (method=='Laplace') method=0 else if (method=='MC') method=1 else if (method=='plugin') method=2 else if (method=='Hybrid') method=3 else stop("Invalid argument 'method'")
  method <- as.integer(method); B <- as.integer(B); logscale <- as.integer(logscale)
  alpha <- as.double(alpha); lambda <- as.double(lambda)
  ngroups= p; nvaringroup= as.integer(rep(1,p))
  ans <- pimomMarginalUI(sel,nsel,n,p,y,sumy2,x,XtX,ytX,tau,method,B,logscale,alpha,lambda,ngroups,nvaringroup)
  return(ans);
}


fimomUNeg <- function(th, XtX, ytX, sumy2, tau, alpha, lambda, n) {
  #Last component in th is eta=log(phi) i.e. log-residual variance
  eta <- th[length(th)]; th <- th[-length(th)]; p <- length(th)
  ss <- lambda + sumy2 - 2*ytX %*% matrix(th,ncol=1) + mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE)
  ans <- tau*exp(eta)*sum(1/th^2) + sum(log(th^2)) + .5*eta*(n-p+alpha) + .5*exp(-eta)*ss
  return(ans)
}

fppimomUNeg <- function(th, XtX, ytX, sumy2, tau, alpha, lambda) {
  eta <- th[length(th)]; th <- th[-length(th)]; p <- length(th)
  ss <- lambda + sumy2 - 2*ytX %*% matrix(th,ncol=1) + mahalanobis(th,center=rep(0,p),cov=XtX,inverted=TRUE)
  ans <- matrix(NA,nrow=p+1,ncol=p+1)
  ans[1:p,1:p] <- exp(-eta) * XtX + diag(6*tau*exp(eta)/th^4 - 2/th^2,ncol=p)
  ans[1:p,p+1] <- ans[p+1,1:p] <- -exp(eta)*2*tau/th^3 - exp(-eta)*as.vector(matrix(th,nrow=1) %*% XtX - ytX)
  ans[p+1,p+1] <- tau*exp(eta)*sum(1/th^2) + .5*exp(-eta)*ss
  return(ans)
}

imomModeU <- function(thini, phiini, XtX, ytX, sumy2, tau, alpha, lambda, n) {
  #piMOM mode when phi is unknown using gradient algorithm
  th <- thini; phi <- phiini
  err <- 1; niter <- 0
  b <- -(n-ncol(XtX)+alpha)/2
  while ((err > 0.001) & (niter<50)) {
    err <- 0; niter <- niter+1
    for (i in 1:length(th)) {
      a <- c(2*tau*phi, 0, -2, (ytX[i]-sum(XtX[i,-i]*th[-i]))/phi, -XtX[i,i]/phi)
      thnew <- polyroot(a)
      thnew <- Re(thnew[abs(Im(thnew))< 1e-7])
      thnew <- thnew[sign(thnew)==sign(th[i])]
      err <- err+abs(th[i]-thnew)
      th[i] <- thnew
    }
    a <- tau*sum(1/th^2)
    b <- .5*(n-ncol(XtX)+alpha)
    c <- -.5*(lambda + sumy2 - 2*ytX %*% matrix(th,ncol=1) + matrix(th,nrow=1) %*% XtX %*% matrix(th,ncol=1))
    d <- sqrt(b^2 - 4*a*c)
    if (-b > d) { phinew <- (-b-d)/(2*a) } else { phinew <- (-b+d)/(2*a) }
    err <- err+abs(phi-phinew)
    phi <- phinew
  }
  return(c(th,log(phi)))
}

imomUIntegralApprox <- function(thini, etaini, XtX, ytX, sumy2, tau, alpha, lambda, n, logscale=TRUE) {
#Laplace approx to product imom marginal (uses gradient search to find mode)
  m <- imomModeU(thini=thini,phiini=exp(etaini),XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda,n=n)
  fopt <- fimomUNeg(th=m,XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda,n=n)
  V <- fppimomUNeg(m,XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda)
  ans <- -fopt - .5*as.numeric(determinant(V,logarithm=TRUE)$modulus) + .5*nrow(XtX)*log(2*tau)
  if (!logscale) ans <- exp(ans)
  return(list(ans=ans,thopt=m,Vopt=V,objective=fopt))
}


imomUIntegralApproxOld <- function(thini, etaini, XtX, ytX, sumy2, tau, alpha, lambda, n, logscale=TRUE) {
#Laplace approx to product imom marginal (uses numerical optimizer to find mode)
  opt <- nlminb(c(thini,etaini), objective=fimomUNeg, XtX=XtX, ytX=ytX, sumy2=sumy2, tau=tau, alpha=alpha, lambda=lambda, n=n)
  V <- fppimomUNeg(opt$par,XtX=XtX,ytX=ytX,sumy2=sumy2,tau=tau,alpha=alpha,lambda=lambda)
  ans <- -opt$objective - .5*as.numeric(determinant(V,logarithm=TRUE)$modulus) + .5*nrow(XtX)*log(2*tau)
  if (!logscale) ans <- exp(ans)
  return(list(ans=ans,thopt=opt$par,Vopt=V,objective=opt$objective))
}


pimomMarginalUR <- function(y, x, alpha=0.001, lambda=0.001, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under a product imom prior (unknown variance)
# - y: response variable
# - x: design matrix
# - alpha, lambda: prior for phi is IGamma(alpha/2,lambda/2)
# - tau: prior dispersion parameter
# - method: method to approximate the integral for known phi. Integral wrt phi is performed via integrate. 'Laplace' for Laplace approx which may underestimate true value, 'MC' for exact evaluation which can be very computationally expensive. 'Hybrid' uses numerical integration to integrate over phi and Laplace to integrate over theta (it corrects the Laplace error with an exact evaluation for a single value of phi close to the posterior mode)
# - B: number of Monte Carlo samples to use (ignored if method=='Laplace')
  #require(actuar)
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  if (ncol(x)==0) {
    n <- length(y)
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    if (missing(XtX)) { XtX <- t(x) %*% x }
    if (missing(ytX)) { ytX <- t(y) %*% x }
    if (method=='Hybrid') {
      f2int <- function(z,method, adj=1) {
        ans <- double(length(z))
        for (i in 1:length(ans)) ans[i] <- pimomMarginalKR(y=y,x=x,phi=z[i],tau=tau,method=method,B=B,logscale=FALSE,XtX=XtX,ytX=ytX)
        ans <- ans * dinvgamma(z, alpha/2, scale=lambda/2) * adj
        return(ans)
      }
      #Compute adjustment factor
      e <- y - x %*% solve(XtX + diag(tau,nrow=nrow(XtX))) %*% t(ytX)
      phiest <- (sum(e^2)+lambda)/(length(y)+alpha)
      intmc <- f2int(phiest,method='MC')
      intlapl <- f2int(phiest,method='Laplace')
      adj <- intmc/intlapl
      ans <- log(integrate(f2int, 0, Inf, method='Laplace', adj=adj)$value)
    } else {
      thini <- as.vector(solve(XtX + tau*diag(nrow(XtX))) %*% t(ytX))
      e <- y - x %*% matrix(thini,ncol=1)
      etaini <- log((sum(e^2)+lambda)/(length(y)+alpha))
      ans <- imomUIntegralApprox(thini=thini,etaini=etaini,XtX=XtX,ytX=ytX,sumy2=sum(y^2),tau=tau,alpha=alpha,lambda=lambda,n=nrow(x),logscale=TRUE)
      ans <- ans$ans + .5*alpha*log(.5*lambda) - .5*nrow(x)*log(2*pi) - lgamma(.5*alpha)
    }
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}



###
### margpemom.R
###

lbesselK <- function(x, nu) log(besselK(x,nu=nu,expon.scaled=TRUE)) - x

femomNeg <- function(th, m, S, phi, tau, logscale=TRUE) .5*mahalanobis(th, center=m, cov=S, inverted=TRUE)/phi + tau*phi*sum(1/th^2)
fpemomNeg <- function(th, m, S, phi, tau) S %*% matrix(th-m, ncol=1)/phi - 2*tau*phi*sum(1/th^3)
fppemomNeg <- function(th, m, S, phi, tau) S/phi + 6*tau*phi*diag(1/th^4,nrow=length(th))

pemomIntegralApproxR <- function(m, S, phi, tau, logscale=TRUE) {
  #Laplace approx to integral N(th; m, phi*solve(S)) prod(exp(-tau*phi/th^2)) wrt th
  opt <- nlminb(m, objective=femomNeg, gradient=fpemomNeg, m=m, S=S, phi=phi, tau=tau)$par
  fopt <- -femomNeg(opt,m=m,S=S,phi=phi,tau=tau)
  hess <- fppemomNeg(opt,m=m,S=S,phi=phi,tau=tau)
  ans <- fopt + .5*log(det(S)) - .5*log(det(hess)) - .5*length(m)*log(phi)
  if (!logscale) ans <- exp(ans)
  return(ans)
}



pemomMarginalKR <- function(y, x, phi, tau, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product eMOM prior (variance phi known)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior proportional to N(th; 0, tau*phi*I) * prod(exp(-tau*phi/th^2)^r
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm)
  n <- length(y); p <- ncol(x)
  if (p==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    ans <- -.5*(sum(y^2) - t(m) %*% S %*% m)/phi + p*sqrt(2)  - .5*n*log(2*pi*phi) - .5*p*log(tau) - log(sqrt(det(S)))
    if (method=='Laplace') {
      I <- pemomIntegralApproxR(m=m, S=S, phi=phi, tau=tau, logscale=TRUE)
    } else if (method=='1storder') {
      I <- - tau*phi*sum(1/m^2)
    } else if (method=='MC') {
      thsim <- rmvnorm(B,m,phi*solve(S))
      eprod <- exp(-tau*phi*rowSums(1/thsim^2))
      I <- log(mean(eprod))
    }
    ans <- ans + I
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}




pemomMarginalUR <- function(y, x, r, alpha=0.001, lambda=0.001, tau, method='Laplace', B=10^5, logscale=TRUE) {
  #Marginal likelihood for product eMOM prior (variance phi unknown)
  # - Likelihood: y ~ N(x %*% th, phi * I)
  # - Prior for th proportional to N(th; 0, tau*phi*I) * prod(exp(-tau*phi/th^2)^r
  # - Prior for phi: IGamma(alpha/2,lambda/2)
  #   i.e. phi is the residual variance; tau the prior dispersion parameter
  #require(mvtnorm)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  n <- length(y); p <- ncol(x)
  if (ncol(x)==0) {
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    S <- t(x) %*% x + diag(p)/tau
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    apost <- n + alpha
    lpost <- as.numeric(lambda + sum(y^2) - t(m) %*% S %*% m)
    #
    if (method=='Laplace') {
      pen <- lpost*tau*sum(1/m^2)
      I <- log(2) - lgamma(.5*apost) + .25*apost*log(.5*pen) + lbesselK(sqrt(2*pen),nu=.5*apost)
    } else if (method=='1storder') {
      phi <- lpost/apost
      I <- - tau*phi*sum(1/m^2)
    } else if (method=='MC') {
      phisim <- 1/rgamma(B, apost, lpost)
      cholV <- t(chol(solve(S)))
      z <- rmvnorm(B,rep(0,p),diag(p))
      thsim <- as.vector(m) + (cholV %*% t(z)) * sqrt(phisim)
      eprod <- -tau*phisim*colSums(1/thsim^2)
      offset <- max(eprod)
      I <- log(mean(exp(eprod-offset))) + offset
    } else {
      stop("Only 'Laplace', '1storder' and 'MC' methods are implemented")
    }
    #
    num <- p*sqrt(2) + .5*alpha*log(lambda/2) + lgamma(apost/2)
    den <- .5*n*log(2*pi) + .5*log(det(S)) + .5*p*log(tau) + lgamma(alpha/2) + .5*apost*log(lpost/2)
    ans <- I + num - den
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
