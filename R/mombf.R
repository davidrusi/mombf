##############################################################################################
##
## PAIRWISE BAYES FACTORS BASED ON QUADRATIC PENALTIES
##
##############################################################################################


###
### mombf.R
###

mombf <- function(lm1,coef,g,prior.mode,baseDensity='normal',nu=3,theta0,logbf=FALSE,B=10^5) {
  if (!(baseDensity %in% c('normal','t'))) stop('baseDensity must either be normal or t')
  if (baseDensity=='t') {
    if (length(g)>1) stop("For baseDensity=='t' g is only allowed to have length 1")
    if (nu <= 2) stop('For t base density nu must be >=3, otherwise prior is improper')
    if (nu >= length(lm1$residuals)-length(coef)+2) warning(paste('Based on the amount of data available, nu should be <',length(lm1$residuals)-length(coef)+2,'in order to guarantee finite sample consistency'))
  }
  UseMethod("mombf")
}

###
### mombf.lm.R
###

mombf.lm <- function(lm1,coef,g,prior.mode,baseDensity='normal',nu=3,theta0,logbf=FALSE,B=10^5) {

if ((!missing(g)) & (!missing(prior.mode))) warning('Both g and prior.mode were specified. g will be ignored')
if ((missing(g)) & (missing(prior.mode))) stop('Either g or prior.mode must be specified')
if (missing(theta0)) theta0 <- rep(0,length(coef)) else if (length(theta0)!=length(coef)) stop('theta0 must have the same length as coef')

  thetahat <- coef(lm1)
  V <- summary(lm1)$cov.unscaled
  n <- length(lm1$residuals); p <- length(thetahat); p1 <- length(coef)
  if ((min(coef)<1) | (max(coef)>p)) stop('Non-valid value for coef. Use only values between 1 and the number of coefficients in lm1')
  ssr <- sum(residuals(lm1)^2); sr <- sqrt(ssr/(n-p))
  if (!missing(prior.mode)) g <- mode2g(prior.mode,prior=paste(baseDensity,'Mom',sep=''),nu=nu)
  if (baseDensity=='normal') {
    bf.mom <- momunknown(thetahat[coef],V[coef,coef],n=n,g=g,ssr=ssr,nuisance.theta=p-p1,theta0=theta0,logbf=logbf)
  } else if (baseDensity=='t') {
    randomg <- 1/rgamma(B,nu/2,nu/2)
    bf.mom <- momunknown(thetahat[coef],V[coef,coef],n=n,g=g*randomg,ssr=ssr,nuisance.theta=p-p1,theta0=theta0,logbf=logbf)
    bf.mom <- mean(bf.mom*randomg)*(nu-2)/nu
  }
  return(bf.mom)
}


###
### momknown.R
###

momknown <- function(theta1hat,V1,n,g=1,theta0,sigma,logbf=FALSE) {
if (missing(sigma)) stop('sigma must be specified')
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat)
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1) * n*g/((1+n*g)*sigma^2) #noncentr param
muk <- p1+l
t1 <- matrix(theta1hat-theta0,nrow=1) %*% solve(V1) %*% matrix(theta1hat-theta0,ncol=1) * n*g/((1+n*g)*sigma^2)
bf <- .5*t1 + log(muk) - log(1+n*g) - (p1/2)*log(1+n*g) - log(p1)
if (!logbf) bf <- exp(bf)
return(bf)
}

###
### momunknown.R
###

momunknown <- function(theta1hat,V1,n,nuisance.theta,g=1,theta0,ssr,logbf=FALSE) {
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat); p <- p1 + nuisance.theta
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1)
sigma2hat <- (ssr + l/(1+n*g))/(n-nuisance.theta)
muk <- p1+ l* n*g/((1+n*g)*sigma2hat)
bf <- (-(n-nuisance.theta)/2)*log(1+n*g*ssr/(ssr+l)) + log(muk) - log(1+n*g) + ((n-p)/2)*log(1+n*g) - log(p1)
if (!logbf) bf <- exp(bf)
return(bf)
}




###
### imombf.R
###

imombf <- function(lm1,
                   coef,
                   g,
                   prior.mode,
                   nu=1,
                   theta0,
                   method='adapt',
                   nquant=100,
                   B=10^5) {
    UseMethod("imombf")
}


###
### imombf.lm.R
###

imombf.lm <- function(lm1,coef,g,prior.mode,nu=1,theta0,method='adapt',nquant=100,B=10^5) {
if ((!missing(g)) & (!missing(prior.mode))) warning('Both g and prior.mode were specified. g will be ignored')
if ((missing(g)) & (missing(prior.mode))) stop('Either g or prior.mode must be specified')
if (missing(theta0)) theta0 <- rep(0,length(coef)) else if (length(theta0)!=length(coef)) stop('theta0 must have the same length as coef')

  thetahat <- coef(lm1)
  V <- summary(lm1)$cov.unscaled
  n <- length(lm1$residuals); p <- length(thetahat); p1 <- length(coef)
  if ((min(coef)<1) | (max(coef)>p)) stop('Non-valid value for coef. Use only values between 1 and the number of coefficients in lm1')
  ssr <- sum(residuals(lm1)^2); sr <- sqrt(ssr/(n-p))
  if (!missing(prior.mode)) g <- mode2g(prior.mode,prior='iMom')
  bf.imom <- imomunknown(thetahat[coef],V[coef,coef],n,nuisance.theta=p-p1,g=g,nu=nu,theta0=theta0,ssr=ssr,method=method,nquant=nquant,B=B)
  return(bf.imom)
}

###
### imomknown.R.R
###

imomknown <- function(theta1hat,V1,n,nuisance.theta,g=1,nu=1,theta0,sigma,method='adapt',B=10^5) {
if (missing(sigma)) stop('sigma must be specified')
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
f <- function(z) { ans <- (n*gi/z)^((nu+p1)/2) * exp(-n*gi/z); ans[z==0] <- 0; return(ans) }

p1 <- length(theta1hat)
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1) / sigma^2 #noncentr param
m <- double(length(g))
if (method=='MC') {
  z <- rchisq(B,df=p1,ncp=l)
  for (i in 1:length(m)) { gi <- g[i]; m[i] <- mean(f(z)) }
} else if (method=='adapt') {
  f2 <- function(z) { return(f(z)*dchisq(z,df=p1,ncp=l)) }
  for (i in 1:length(m)) { gi <- g[i]; m[i] <- integrate(f2,0,Inf)$value }
} else {
  stop('method must be adapt or MC')
}
bf <- exp((p1/2)*log(2/(n*g)) + lgamma(p1/2)-lgamma(nu/2) + .5*l) * m
return(bf)
}


###
### imomunknown.R
###

imomunknown <- function(theta1hat,V1,n,nuisance.theta,g=1,nu=1,theta0,ssr,method='adapt',nquant=100,B=10^5) {

fncp <- function(sigma2) {
  l <- theta1hat-theta0; l <- as.vector(matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1)) / sigma2
  l[l==Inf] <- exp(80)
  return(l)
}
f <- function(z) { ans <- (n*gi/z)^((nu+p1)/2) * exp(-n*gi/z); ans[z==0] <- 0; return(ans) }

if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat)
m <- double(length(g))
if (method=='MC') {
  sigma2 <- 1/rgamma(B,(n-nuisance.theta)/2,ssr/2)
  l <- fncp(sigma2)
  z <- rchisq(B,df=p1,ncp=l)
  for (i in 1:length(m)) { gi <- g[i]; m[i] <- mean(f(z)) }
} else if (method=='adapt') {
  f2 <- function(z,z2) { return(f(z)*dchisq(z,df=p1,ncp=fncp(z2))) }
  qseq <- 1/qgamma((2*(1:nquant)-1)/(2*nquant),(n-nuisance.theta)/2,ssr/2)
  for (i in 1:length(m)) {
    m[i] <- 0; for (j in 1:nquant) { gi <- g[i]; m[i] <- m[i]+integrate(f2,0,Inf,z2=qseq[j])$value }
    m[i] <- m[i]/nquant
  }
} else {
  stop('method must be adapt or MC')
}
t1 <- theta1hat-theta0; t1 <- matrix(t1,nrow=1) %*% solve(V1) %*% matrix(t1,ncol=1) / ssr #noncentr param
bf <- exp((p1/2)*log(2/(n*g)) + lgamma(p1/2) - lgamma(nu/2) + ((n-nuisance.theta)/2)*log(1+t1)) * m
return(bf)
}





