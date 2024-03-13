##############################################################################################
##
## DENSITY, QUANTILES, CDF AND PRIOR ELICITATION FOR NON-LOCAL PRIORS
##
##############################################################################################


###
### dmom.R
###

#Wrapper to call dpmom (product mom) and dqmom (quadratic mom)
dmom <- function(x, tau, a.tau, b.tau, phi=1, r=1, V1, baseDensity='normal', nu=3, logscale=FALSE, penalty='product') {
  if (penalty=='product') {
    ans <- dpmom(x, tau=tau, a.tau=a.tau, b.tau=b.tau, phi=phi, r=r, baseDensity=baseDensity, logscale=logscale)
  } else if (penalty=='quadratic') {
    if (r>1) stop("r>1 not implemented for penalty=='quadratic'. Try penalty=='product' instead")
    ans <- dqmom(x, V1=V1, g=tau, n=1, baseDensity=baseDensity, nu=nu)
  } else {
    stop("Only 'penalty==product' and 'penalty==quadratic' are implemented")
  }
  return(ans)
}

##Product MOM density

setMethod("dpmom", signature(x='vector'), function(x, tau, a.tau, b.tau, phi=1, r=1, baseDensity='normal', logscale=FALSE) {
    if (missing(tau) & (missing(a.tau) | missing(b.tau))) stop("Either tau or (a.tau,b.tau) must be specified")
    if (baseDensity=='normal') {
        if (!missing(tau)) {
            ans <- dnorm(x,0,sd=sqrt(tau*phi),log=TRUE) + r*log(x^2/(tau*phi)) - sum(log(seq(1,2*r-1,by=2)))
        } else {
            ct <- lgamma(r+.5*(a.tau+1)) + r*log(2) - lgamma(.5*a.tau) - .5*log(pi) - (r+.5)*log(b.tau) - .5*log(phi)
            ans <- ct + r*log(x^2/phi) - sum(log(seq(1,2*r-1,by=2))) - (r+.5*(a.tau+1))*log(1+x^2/(b.tau*phi))
        }
    } else if (baseDensity=='laplace') {
        if (r!=1) stop("Only r=1 is implemented for baseDensity='normal'")
        ans= log(x^2) + dalapl(x,0,scale=tau*phi,logscale=TRUE) - log(2*tau*phi)
    } else { stop("Only baseDensity=='normal' or 'laplace' is implemented for the product MOM") }
  if (!logscale) ans <- exp(ans)
  ans
}
)
setMethod("dpmom", signature(x='matrix'), function(x, tau, a.tau, b.tau, phi=1, r=1, baseDensity='normal', logscale=FALSE) {
  if (baseDensity!='normal') stop("Only baseDensity=='normal' is implemented for the product MOM")
  if (missing(tau) & (missing(a.tau) | missing(b.tau))) stop("Either tau or (a.tau,b.tau) must be specified")
  p <- ncol(x)
  normct <- p*sum(log(seq(1,2*r-1,by=2)))
  distval <- rowSums(x^2)
  if (!missing(tau)) {
    ans <- -(p * log(2 * pi) + p*(log(phi)+log(tau))  + distval/(phi*tau))/2 + r*rowSums(log(x^2/(tau*phi))) - normct
  } else {
    anew <- r*p + .5*p + .5*a.tau
    num <- lgamma(anew) - anew*log(1+rowSums(x^2)/(phi*b.tau))  + r*rowSums(log(2*x^2/(phi*b.tau))) - normct
    den <- .5*p*(log(pi)+log(phi)+log(b.tau)) + lgamma(.5*a.tau)
    ans <- num - den
  }
  if (!logscale) ans <- exp(ans)
  ans
}
)

##Quadratic MOM density
dqmom <- function(x,V1=1,g,n=1,theta0,baseDensity='normal',nu=3) {
if (missing(g)) stop("Prior dispersion must be specified for quadratic MOM")
if (!(baseDensity %in% c('normal','t'))) stop("The only implemented baseDensity values are 'normal' and 't'")
if (baseDensity=='t' & nu<3) stop('nu must be >=3, otherwise the prior is improper')
if (missing(V1)) {
  if (is.vector(x)) V1 <- 1 else V1 <- diag(ncol(x))
}
if (missing(theta0)) {
  if (is.vector(V1)) theta0 <- 0 else theta0 <- rep(0,ncol(V1))
}
    
if (is.vector(V1)) {
  qtheta <- (x-theta0)^2/(n*g*V1)
  if (baseDensity=='normal') {
      ans <- qtheta*dnorm(x,theta0,sd=sqrt(n*g*V1))
  } else if (baseDensity=='t') {
    normct <- exp(lgamma(.5*(nu+1))-lgamma(.5*nu)-.5*log(nu*pi*n*g*V1))
    ans <- qtheta*normct*(1+qtheta/nu)^(-.5*(nu+1))*(nu-2)/nu
  }
} else {
  #require(mvtnorm)
  qtheta <- mahalanobis(x,center=theta0,cov=n*g*V1)
  if (baseDensity=='normal') {
    ans <- qtheta * mvtnorm::dmvnorm(x,mean=theta0,sigma=n*g*V1)/ncol(V1)
  } else if (baseDensity=='t') {
    ans <- qtheta * mvtnorm::dmvt(x,delta=theta0,sigma=n*g*V1)*(nu-2)/(nu*length(theta0))
  }
}
return(ans)
}



###
### dimom.R
###

#Wrapper to call dpimom (product iMOM) and dqimom (quadratic imom)
dimom <- function(x, tau=1, phi=1, V1, logscale=FALSE, penalty='product') {
  if (penalty=='product') {
    ans <- dpimom(x, tau=tau, phi=phi, logscale=logscale)
  } else if (penalty=='quadratic') {
    ans <- dqimom(x, V1=V1, g=tau, n=1, nu=1, logscale=logscale)
  } else {
    stop("Only 'penalty==product' and 'penalty==quadratic' are implemented")
  }
  return(ans)
}

#Product iMOM

setMethod("dpimom", signature(x='vector'), function(x, tau=1, phi=1, logscale=FALSE) {
  ans <- .5*(log(tau)+log(phi)) - lgamma(.5) - log(x^2) - tau*phi/x^2
  ans[is.nan(ans)] <- -Inf
  if (!logscale) ans <- exp(ans)
  ans
}
)

setMethod("dpimom", signature(x='matrix'), function(x, tau=1, phi=1, logscale=FALSE) {
  x2 <- x^2
  ans <- ncol(x)*(.5*(log(tau)+log(phi)) - lgamma(.5)) - rowSums(log(x2)) - tau*phi*rowSums(1/x2)
  ans[is.nan(ans)] <- -Inf
  if (!logscale) ans <- exp(ans)
  ans
}
)

setMethod("dpimom", signature(x='data.frame'), function(x, tau=1, phi=1, logscale=FALSE) {
  dpimom(as.matrix(x),tau=tau,phi=phi,logscale=logscale)
}
)


#Quadratic iMOM
dqimom <- function(x,V1=1,g=1,n=1,nu=1,theta0,logscale=FALSE) {
if (is.vector(x)) {
  if (missing(theta0)) theta0 <- 0
  if (missing(V1)) V1 <- 1
  qtheta <- (x-theta0)^2/(n*g*V1)
  k <- -0.5*log(n*g*V1) - lgamma(nu/2)
  p1 <- 1
} else {
  if (missing(theta0)) theta0 <- rep(0,ncol(x))
  if (missing(V1)) V1 <- diag(ncol(x))
  qtheta <- t(matrix(x,nrow=nrow(x))) - theta0
  qtheta <- qtheta %*% solve(n*g*V1) %*% t(qtheta)
  k <- -0.5*log(n*g*det(V1)) - lgamma(nu/2) + lgamma(ncol(x)/2) - .5*ncol(x)*log(pi)
  p1 <- ncol(x)
}
ans <- (k - .5*(nu+p1)*log(qtheta) -1/qtheta)
if (!logscale) { ans <- exp(ans); ans[is.na(ans)] <- 0 }
return(ans)
}




###
### pmom.R
###

pmom <- function(q,V1=1,tau=1) {

  z <- .5-(pnorm(abs(q)/sqrt(V1*tau)) - abs(q)/sqrt(2*pi*V1*tau) * exp(-.5*q^2/(tau*V1)) - .5)
  return(z*(q<=0)+(1-z)*(q>0))
}


###
### qmom.R
###

qmom <- function(p, V1=1, tau=1) {

  e <- function(q) { return((pmom(q, V1, tau)-pneg)^2) }
  ans <- double(length(p))
  for (i in 1:length(p)) {
    pneg <- ifelse(p[i]<=.5, p[i], 1-p[i])
    ans[i] <- nlminb(start=-1, objective=e)$par
    if (p[i]>.5) {
      ans[i] <- -ans[i]
    }
  }
  ans
}


###
### pimom.R
###

pimom <- function(q,V1=1,tau=1,nu=1) {

  z <- (.5*pgamma(tau*V1/q^2,nu/2,1))
  return(z*(q<=0)+(1-z)*(q>0))
}


###
### qimom.R
###

qimom <- function(p, V1=1, tau=1, nu=1) {

  ans <- double(length(p))
  ans[p<=.5] <- -sqrt(tau*V1/qgamma(2*p[p<=.5], nu/2, 1))
  ans[p>.5] <- sqrt(tau*V1/qgamma(2*(1-p[p>.5]), nu/2, 1))
  ans
}



###
### emom.R
###


setMethod("demom",signature(x='vector'),function(x, tau, a.tau, b.tau, phi=1, logscale=FALSE) {
  V1 <- 1
  if (!missing(tau)) {
    pen <- -tau*phi/x^2
    normct <- sqrt(2)
    ans <- pen + dnorm(x,mean=0,sd=sqrt(tau*phi*V1),log=TRUE) + normct
  } else {
    p <- 1; x2phi <- x^2/phi
    anew <- .5*(a.tau+p); bnew <- .5*(b.tau+x2phi)
    num <- sqrt(2)*p + .5*a.tau*log(.5*b.tau)
    den <- lgamma(.5*a.tau) + .5*p*(log(2*pi)+log(phi)) + anew*log(bnew)
    bt <- bnew/x2phi
    ans <- num - den + log(2) + .5*anew*log(bt) + log(besselK(sqrt(4*bt),nu=anew,expon.scaled=TRUE)) - sqrt(4*bt)
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
)

setMethod("demom",signature(x='data.frame'),function(x, tau, a.tau, b.tau, phi=1, logscale=FALSE) {
  demom(as.matrix(x),tau=tau,a.tau=a.tau,b.tau=b.tau,phi=phi,logscale=logscale)
}
)

setMethod("demom",signature(x='matrix'),function(x, tau, a.tau, b.tau, phi=1, logscale=FALSE) {
  p <- ncol(x)
  V1 <- diag(p)
  if (!missing(tau)) {
    pen <- -tau*phi*rowSums(1/x^2)
    normct <- p*sqrt(2)
    ans <- pen + mvtnorm::dmvnorm(x,mean=rep(0,p),sigma=tau*phi*V1,log=TRUE) + normct
  } else {
    anew <- .5*(a.tau+p); bnew <- .5*(b.tau+rowSums(x^2)/phi)
    num <- sqrt(2)*p + .5*a.tau*log(.5*b.tau)
    den <- lgamma(.5*a.tau) + .5*p*(log(2*pi)+log(phi)) + anew*log(bnew)
    bt <- bnew*phi*rowSums(1/x^2)
    ans <- num - den + log(2) + .5*anew*log(bt) + log(besselK(sqrt(4*bt),nu=anew,expon.scaled=TRUE)) - sqrt(4*bt)
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}
)


pemom <- function(q, tau, a.tau, b.tau) integrate(demom,-Inf,q,tau=tau,a.tau=a.tau,b.tau=b.tau)$value



##############################################################################################
## PRIOR ELICITATION
##############################################################################################


###
### g2mode.R
###

g2mode <- function(g,
                   prior=c('iMom', 'normalMom', 'tMom'),
                   nu=1,
                   dim=1) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           2*g
         },
         tMom = {
           g*2*nu/(nu-2+dim)
         },
         iMom = {
           2*g/(nu+dim)
         })
}


###
### mode2g.R
###

mode2g <- function(prior.mode,
                   prior=c('iMom', 'normalMom', 'tMom'),
                   nu=1,
                   dim=1) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           prior.mode / 2
         },
         iMom = {
           prior.mode * (nu+dim) / 2
         },
         tMom = {
           if (nu<3) {
             stop('tMom prior must have nu>2 degrees of freedom')
           }
           prior.mode * (nu-2+dim) / (2*nu)
         })
}


###
### priorp2g.R
###

priorp2g <- function(priorp,
                     q,
                     nu=1,
                     prior=c('iMom', 'normalMom','tMom')) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           e <- function(logg) {
             return((1-2*pmom(-abs(q), tau=exp(logg)) - priorp[i])^2)
           }

           ans <- double(length(priorp))
           for (i in 1:length(priorp)) {
             ans[i] <- exp(nlminb(start=0, objective=e)$par)
           }
           ans
         },
         tMom = {
           stop("prior=='tMom' is not currently implemented")
         },
         iMom = {
           p <- (1-priorp)/2
           qgamma(2*p,nu/2,1)*q^2
         })
}




##############################################################################################
## MARGINAL NLP * IG PRIORS
##############################################################################################


dmomigmarg= function(x,tau,a,b,logscale=FALSE) {
#Marginal MOM density p(x)= int MOM(x;0,tau*phi) IG(phi;a/2,b/2) dphi
    ans= log(2) + lgamma((a+3)/2) - 0.5*log(pi) - 1.5*log(b*tau) - lgamma(a/2) + log(x^2) - ((a+3)/2) * log(1 + x^2/(b*tau))
    if (!logscale) ans= exp(ans)
    return(ans)
}

#Marginal MOM cdf 
pmomigmarg= function(x,tau,a,b) { integrate(dmomigmarg,-Inf,x,tau=tau,a=a,b=b)$value }



#Marginal eMOM density p(x)= int eMOM(x;0,tau*phi) IG(phi;a/2,b/2) dphi
demomigmarg= function(x,tau,a,b,logscale=FALSE) {
    apos= (a+1)/2; bpos= (b+x^2/tau)/2
    ans= log(mgfIG(-tau/x^2,apos,bpos))
    ans= ans + sqrt(2) + lgamma((a+1)/2) - lgamma(a/2) - 0.5*log(pi*tau*b) - 0.5*(a+1) * log(1+x^2/(tau*b))
    if (!logscale) ans= exp(ans)
    return(ans)
}

#Marginal MOM cdf 
pemomigmarg= function(x,tau,a,b) { integrate(demomigmarg,-Inf,x,tau=tau,a=a,b=b)$value }


#Moment-generating function of IG(a,b) evaluated at t
mgfIG= function(t,a,b) {
    if (length(a)==1) a= rep(a,length(t))
    if (length(b)==1) b= rep(b,length(t))
    ans= double(length(t))
    for (i in 1:length(t)) {
        f= function(x) { exp(t[i] * x) * dinvgamma(x,a[i],b[i]) }
        ans[i]= integrate(f,0,Inf)$value
    }
    return(ans)
}



minbeta2tauMOMIGmarg= function(minbeta,a=3,b=3,targetprob=0.99) {
#Find tau such that P(|x| > minbeta)= targetprob, where p(x) is the MOM-IG marginal
    if (minbeta<0) stop("minbeta must be >0")
    f= function(z) abs(targetprob - 2*pmomigmarg(-minbeta,tau=z,a=a,b=b))
    optimize(f,interval=c(.01,5))$minimum
}


minbeta2taueMOMIGmarg= function(minbeta,a=3,b=3,targetprob=0.99) {
#Find tau such that P(|x| > minbeta)= targetprob, where p(x) is the MOM-IG marginal
    if (minbeta<0) stop("minbeta must be >0")
    f= function(z) abs(targetprob - 2*pemomigmarg(-minbeta,tau=z,a=a,b=b))
    tauseq= c(seq(.01,.3,length=10),seq(.4,2,by=.1))
    fseq= sapply(tauseq, f)
    o= order(fseq)
    optimize(f,interval=tauseq[o[1:2]])$minimum
}



##############################################################################################
## CONTINUOUS SPIKE & SLAB PRIORS
##############################################################################################


dmomss= function(x,tau0,tau1,pspike=0.5,baseDensity="normal",logscale=FALSE) {
#Spike & Slab MOM prior
# (1-pspike) * p(x; 0,tau0) + pspike * (x^2/c) * p(x; 0,tau1), where p(x;0,tauj) is Normal for baseDensity="normal" and Laplace for baseDensity="laplace" and c the MOM normalization constant
    if (baseDensity=="normal") {
        d0= dnorm(x,0,tau0,log=TRUE) + log(1-pspike)
        d1= dmom(x,tau=tau1,phi=1,baseDensity="normal",logscale=TRUE) + log(pspike)
    } else if (baseDensity=="laplace") {
        d0= dalapl(x,0,scale=tau0,logscale=TRUE) + log(1-pspike)
        d1= dmom(x,tau=tau1,baseDensity="laplace",logscale=TRUE) + log(pspike)
    }
    ans= exp(d0) + exp(d1)
    if (logscale) ans= log(ans)
    return(ans)       
}

dss= function(x,tau0,tau1,pspike=0.5,baseDensity="normal",logscale=FALSE) {
#Spike & Slab prior
# (1-pspike) * p(x; 0,tau0) + pspike * p(x; 0,tau1), where p(x;0,tauj) is Normal for baseDensity="normal" and Laplace for baseDensity="laplace"
    if (baseDensity=="normal") {
        d0= dnorm(x,0,tau0,log=TRUE) + log(1-pspike)
        d1= dnorm(x,0,tau1,log=TRUE) + log(pspike)
    } else if (baseDensity=="laplace") {
        d0= dalapl(x,0,scale=tau0,logscale=TRUE) + log(1-pspike)
        d1= dalapl(x,0,scale=tau1,logscale=TRUE) + log(pspike)
    }
    ans= exp(d0) + exp(d1)
    if (logscale) ans= log(ans)
    return(ans)       
}
