##################################################################################
## Routines to simulate from the posterior distribution of parameters (BMA, and given one model)
## and to sample from non-local prior distributions
##################################################################################

setMethod("rnlp", signature(y='missing',x='missing',m='missing',V='missing',msfit='msfit',outcometype='missing',family='missing'), function(y, x, m, V, msfit, outcometype, family, priorCoef, priorGroup, priorVar, priorprec, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
  if ((msfit$family=='Continuous') & (msfit$family != 'normal')) stop("Posterior sampling only implemented for Normal residuals")
  y= msfit$ystd; x= msfit$xstd
  outcometype= msfit$outcometype; family= msfit$family
  isgroup= rep(FALSE, length(msfit$groups))
  isgroup[which(table(msfit$groups)>1)]= TRUE
  #Draw model
  pp <- postProb(msfit,method=pp)
  modelid <- strsplit(as.character(pp$modelid), split=',')
  ndraws <- as.numeric(rmultinom(1, size=niter, prob=pp$pp))
  sel <- ndraws>0; modelid <- modelid[sel]; ndraws <- ndraws[sel]
  priorCoef= msfit$priors$priorCoef
  priorGroup= msfit$priors$priorGroup
  priorVar= msfit$priors$priorVar
  priorprec= NULL
  #Draw coefficients
  idx <- c(0,cumsum(ndraws))
  if ((outcometype== 'Continuous') && (family== 'normal')) { ##Linear model
    ans <- matrix(0, nrow=niter, ncol=ncol(x)+1)
    for (i in 1:length(modelid)) {
      b <- min(50, ceiling((burnin/niter) * ndraws[i]))
      colsel <- as.numeric(modelid[[i]])
      if (priorCoef@priorDistr == "icarplus") priorprec = msfit$priorprec[colsel,colsel,drop=FALSE]
      ans[(idx[i]+1):idx[i+1],c(colsel,ncol(ans))] <- rnlp(y=y, x=x[,colsel,drop=FALSE], outcometype=outcometype, family=family, priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, priorprec=priorprec, isgroup=isgroup[colsel], niter=b+idx[i+1]-idx[i], burnin=b)
    }
    if (is.null(colnames(x))) nn <- c(paste('beta',1:ncol(x),sep=''),'phi') else nn <- c(colnames(x),'phi')
    colnames(ans)= nn
  } else {                                                   ##GLM or Survival model
    ans <- matrix(0, nrow=niter, ncol=ncol(x))
    for (i in 1:length(modelid)) {
      b <- min(50, ceiling((burnin/niter) * ndraws[i]))
      colsel <- as.numeric(modelid[[i]])
      if (length(colsel)>0) {
        ans[(idx[i]+1):idx[i+1],colsel] <- rnlp(y=y, x=x[,colsel,drop=FALSE], outcometype=outcometype, family=family, priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, isgroup=isgroup[colsel], niter=b+idx[i+1]-idx[i], burnin=b)
      }
    }
    if (is.null(colnames(x))) nn <- paste('beta',1:ncol(x),sep='') else nn <- colnames(x)
    colnames(ans)= nn
  }
  #Return posterior samples in non-standardized parameterization
  ans= unstdcoef(ans,p=ncol(x),msfit=msfit,coefnames=nn)
  return(ans)
}
)


#Return regression parameter estimates in the parameterization of non-standardized X's (i.e. where X does not have 0 mean, unit variance)
# Input
# - bstd: regression parameters, a matrix with columns corresponding to parameters and rows corresponding to different models / MCMC samples
# - p: columns 1:p in bstd contain regression coefficients for covariate effects, this is to consider that subsequent columns could contain nuisance parameters (phi, the residual variance in linear regression)
# - msfit: object returned by modelSelection.
# - coefnames: names that should be assigned to the output columns
#
# Output: matrix with same dimension as bstd, containing parameter estimates for non-standardized X's
unstdcoef <- function(bstd, p, msfit, coefnames) {
  my= msfit$stdconstants[1,'shift']; mx= msfit$stdconstants[-1,'shift']
  sy= msfit$stdconstants[1,'scale']; sx= msfit$stdconstants[-1,'scale']
  ct= (sx==0)
  b= bstd[,1:p,drop=FALSE]
  if (any(!ct)) b[,!ct]= t(t(b[,!ct])*sy/sx[!ct])  #re-scale regression coefficients
  if (any(ct)) {
      #adjust intercept, if already present
      if (any(!ct)) {
          b[,ct]= my + sy*b[,ct] - colSums(t(b[,!ct,drop=FALSE])*mx[!ct])
      } else {
          b[,ct]= my + sy*b[,ct]
      }
      bstd[,1:p]= b
  } else {
      #add intercept, if not already present
      intercept= my - colSums(t(b[,!ct,drop=FALSE])*mx[!ct]) 
      bstd= cbind(intercept,b,bstd[,-1:-p]); colnames(bstd)= c('intercept',coefnames)
  }
  if ('phi' %in% coefnames) bstd[,'phi']= sy^2*bstd[,'phi'] #re-scale residual variance
  return(bstd)
}



setMethod("rnlp", signature(y='ANY',x='matrix',m='missing',V='missing',msfit='msfit',outcometype='missing',family='missing'), function(y, x, m, V, msfit, outcometype, family, priorCoef, priorGroup, priorVar, priorprec, isgroup, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
    warning("If msfit is specified arguments y and x are ignored")
    rnlp(msfit=msfit, priorprec=priorprec, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm')
}
)


setMethod("rnlp", signature(y='ANY',x='matrix',m='missing',V='missing',msfit='missing',outcometype='missing',family='missing'), function(y, x, m, V, msfit, outcometype, family, priorCoef, priorGroup, priorVar, priorprec, isgroup, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
    rnlp(y=y,x=x,family='Continuous',priorprec=priorprec)
}
)


setMethod("rnlp", signature(y='ANY',x='matrix',m='missing',V='missing',msfit='missing',outcometype='character',family='character'), function(y, x, m, V, msfit, outcometype, family, priorCoef, priorGroup, priorVar, priorprec, isgroup, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
  if ((outcometype== 'Continuous') && (family== 'normal')) {  ##Linear model
    ans <- rnlpLM(y=y,x=x,priorCoef=priorCoef,priorGroup=priorGroup,priorVar=priorVar,priorprec=priorprec,isgroup=isgroup,niter=niter,burnin=burnin,thinning=thinning)
  } else if (outcometype=='glm') { #GLM
    ans <- rnlpGLM(y=y,x=x,family=family,priorCoef=priorCoef,priorGroup=priorGroup,priorVar=priorVar,isgroup=isgroup,niter=niter,burnin=burnin,thinning=thinning)
  } else if ((outcometype=='Survival') && (family=='Cox')) {  #Cox model
    ans <- rnlpCox(y=y,x=x,priorCoef=priorCoef,priorGroup=priorGroup,isgroup=isgroup,niter=niter,burnin=burnin,thinning=thinning)
  } else {
    stop(paste("outcometype",outcometype,"and family=",family,"not implemented",sep=""))
  }
  return(ans)
}
)



rnlpLM <- function(y, x, priorCoef, priorGroup, priorVar, priorprec, isgroup, niter=10^3, burnin=round(niter/10), thinning=1) {
    if (missing(priorGroup)) priorGroup= priorCoef
    p <- ncol(x); n <- length(y)
    tau <- as.double(formatmsPriorsMarg(priorCoef=priorCoef, priorGroup=priorCoef, priorVar=priorVar, priorSkew=momprior(tau=0.348), n=n)$tau)
    if (nrow(x) != n) stop('Dimensions of y and x do not match')
    if (priorVar@priorDistr=='invgamma') {
        a_phi <- as.double(priorVar@priorPars['alpha'])
        b_phi <- as.double(priorVar@priorPars['lambda'])
    } else stop("Only invgamma prior for residual variance is currently implemented")
    if ((priorCoef@priorDistr %in% c('pMOM','peMOM','piMOM')) && (priorCoef@priorDistr == priorGroup@priorDistr)) {
      if (priorCoef@priorDistr=='pMOM') {
        prior <- as.integer(0); r <- as.integer(priorCoef@priorPars['r'])
      } else if (priorCoef@priorDistr=='piMOM') {
        prior <- as.integer(1); r <- as.integer(0)
      } else {
        prior <- as.integer(2); r <- as.integer(0)
      }
      if (p==0) {
        ans <- matrix(1/rgamma((niter-burnin)/thinning, .5*(a_phi+n), .5*(b_phi+sum(y^2))), ncol=1)
        colnames(ans) <- 'phi'
      } else {
        ans <- rnlpPostCI_lm(as.integer(niter),as.integer(burnin),as.integer(thinning),as.double(y),as.double(x),as.integer(p),as.integer(r),tau,a_phi,b_phi,prior)
        #ans <- .Call("rnlpPostCI_lm",as.integer(niter),as.integer(burnin),as.integer(thinning),as.double(y),as.double(x),as.integer(p),as.integer(r),tau,a_phi,b_phi,prior)
        ans <- matrix(ans,ncol=p+1)
        if (is.null(colnames(x))) colnames(ans) <- c(paste('beta',1:ncol(x),sep=''),'phi') else colnames(ans) <- c(colnames(x),'phi')
      }
    } else if (priorCoef@priorDistr %in% c('zellner','bic','normalid','icarplus')) {
      if (priorCoef@priorDistr == 'bic') tau= Inf
      if (p==0) {
        ans <- matrix(1/rgamma((niter-burnin)/thinning, .5*(a_phi+n), .5*(b_phi+sum(y^2))), ncol=1)
        colnames(ans) <- 'phi'
      } else {
        if (priorCoef@priorDistr == 'normalid') {
           S <- solve(t(x) %*% x + diag(p)/tau)
        } else if (priorCoef@priorDistr == 'icarplus') {
           S <- solve(t(x) %*% x + priorprec)
        } else {
           S <- solve((1+1/tau) * t(x) %*% x)
        }
        m <- as.vector(S %*% t(x) %*% matrix(y,ncol=1))
        ssr <- sum(y * (y - (x %*% m)))
        phi <- 1 / rgamma((niter - burnin)/thinning, 0.5*(n+a_phi), 0.5*(ssr + b_phi))
        beta <- matrix(rnorm(ncol(x)*(niter-burnin)/thinning),ncol=ncol(x)) %*% t(chol(S)) * phi
        beta <- t(t(beta)+m)
        ans <- cbind(beta, phi)
        if (is.null(colnames(x))) colnames(ans) <- c(paste('beta',1:ncol(x),sep=''),'phi') else colnames(ans) <- c(colnames(x),'phi')
      }
    } else {
        priorCoef <- zellnerprior(); priorGroup= groupzellnerprior()
        ans <- rnlpLM(y=y, x=x, priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, isgroup=isgroup, niter=niter, burnin=burnin, thinning=thinning)
    }    
    return(ans)
}



setMethod("rnlp", signature(y='missing',x='missing',m='numeric',V='matrix',msfit='missing',outcometype='missing',family='missing'), function(y, x, m, V, msfit, outcometype, family, priorCoef, priorGroup, priorVar, priorprec=priorprec, isgroup, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
  if (missing(priorGroup)) priorGroup= priorCoef
  p <- ncol(V)
  tau <- as.double(priorCoef@priorPars['tau'])
  if ((priorCoef@priorDistr %in% c('pMOM','peMOM','piMOM')) && (priorCoef@priorDistr == priorGroup@priorDistr)) {
    if (priorCoef@priorDistr=='pMOM') {
      prior <- as.integer(0); r <- as.integer(priorCoef@priorPars['r'])
    } else if (priorCoef@priorDistr=='piMOM') {
      prior <- as.integer(1); r <- as.integer(0)
    } else if (priorCoef@priorDistr=='peMOM') {
      prior <- as.integer(2); r <- as.integer(0)
    } else stop("This kind of prior is not implemented")
    ans <- rnlpCI(as.integer(niter),as.integer(burnin),as.integer(thinning),as.double(m),as.double(V),as.integer(p),as.integer(r),tau,prior)
    ans <- matrix(ans,ncol=p)
  } else if (priorCoef@priorDistr %in% c('zellner','bic')) {
    ans <- t(m + t(chol(V)) %*% matrix(rnorm(p*(niter-burnin)/thinning),nrow=p))
  } else stop("This kind of prior is not implemented")
  if (is.null(colnames(ans))) {
      if (is.null(names(m))) colnames(ans) <- paste('beta',1:length(m),sep='') else colnames(ans) <- names(m)
  }
  return(ans)
}
)



glmfamilycode <- function(family) {
    #Convert string indicating GLM family to object of class family, as required by function glm
    if (length(grep("binomial",family))>0) {
        link= strsplit(family,split=" ")[[1]][2]
        if (is.na(link)) link= 'logit'
        ans= binomial(link=link)
    } else if (length(grep("gamma",family))>0) {
        link= strsplit(family,split=" ")[[1]][2]
        ans= Gamma(link=link)
        if (is.na(link)) link='inverse'
    } else if (length(grep("inverse.gaussian",family))>0) {
        link= strsplit(family,split=" ")[[1]][2]
        if (is.na(link)) link="1/mu^2"
        ans= inverse.gaussian(link=link)
    } else if (family=='normal') {
        ans= gaussian()
    } else if (length(grep("poisson",family))>0) {
        link= strsplit(family,split=" ")[[1]][2]
        if (is.na(link)) link="log"
        ans= poisson(link=link)
    } else { stop(paste("family",family,"not recognized",sep=" ")) }
    return(ans)
}



rnlpGLM <- function(y, x, family, priorCoef, priorGroup, priorVar, isgroup=rep(FALSE,ncol(x)), niter=10^3, burnin=round(niter/10), thinning=1) {
#Approximate posterior sampling for GLMs
    p <- ncol(x)
    if (p==0) {
      ans= matrix(double(0),nrow=(niter-burnin)/thinning,ncol=0)
    } else {
      fam= glmfamilycode(family)
      fit= glm(y ~ -1 + ., data=data.frame(x), family=fam)
      pm= postmomentsGLM(fit=fit, priorCoef=priorCoef, priorGroup=priorGroup, isgroup=isgroup)
      ans= rnlp(m=pm$m, V=pm$S, priorCoef=priorCoef, priorGroup=priorGroup, niter=niter, burnin=burnin, thinning=thinning)
    }
    return(ans)
}


postmomentsGLM <- function(fit, priorCoef, priorGroup, isgroup) {
# Extract approximate posterior mean and covariance from GLM fit
# Input
#  - fit: object returned by glm
#  - priorCoef: prior distribution on coefficients
#  - priorGroup: prior on coefficient groups
# Output: list with two elements, m containing the posterior mean and S the posterior covariance
  if (priorCoef@priorDistr != priorGroup@priorDistr) stop("priorGroup != priorCoef not currently implemented")
  tau <- as.double(priorCoef@priorPars['tau'])
  taugroup <- as.double(priorGroup@priorPars['tau'])
  thhat <- matrix(coef(fit),ncol=1)
  if (priorCoef@priorDistr == 'zellner') {
      Vinv <- (t(fit$R) %*% fit$R)
      S <- solve((t(fit$R) %*% fit$R)) / (1+1/tau)
      m <- S %*% Vinv %*% thhat
  } else if (priorCoef@priorDistr == 'bic') {
      S <- vcov(fit)
      m <-  thhat
  } else {
      #if (priorCoef@priorDistr %in% c('pMOM','peMOM','piMOM')) {
      Vinv <- (t(fit$R) %*% fit$R)
      d <- rep(1,nrow(thhat)); d[!isgroup] <- d[!isgroup]/tau; d[isgroup] <- d[isgroup] / taugroup
      D <- diag(nrow(thhat)); diag(D) <- d
      S <- solve(Vinv + D)
      #S <- solve(Vinv + diag(nrow(thhat))/tau)
      m <- S %*% Vinv %*% thhat
  }
  return(list(m=as.vector(m), S=S))
}

rnlpCox <- function(y, x, priorCoef, priorGroup, isgroup=isgroup, niter=10^3, burnin=round(niter/10), thinning=1) {
    tau <- as.double(priorCoef@priorPars['tau'])
    taugroup <- as.double(priorGroup@priorPars['tau'])
    p <- ncol(x)
    if (p==0) {
      ans <- matrix(double(0),nrow=(niter-burnin)/thinning,ncol=0)
    } else {
      fit <- survival::coxph(y ~ ., data=data.frame(x))
      thhat <- matrix(coef(fit),ncol=1); Vinv <- solve(fit$var)
      if (priorCoef@priorDistr == 'zellner') {
        Sinv <- Vinv * (1+1/tau)
      } else if (priorCoef@priorDistr == 'bic') {
        Sinv <- Vinv
      } else {
        #if (priorCoef@priorDistr %in% c('pMOM','peMOM','piMOM')) {
        d <- rep(1,p); d[!isgroup] <- d[!isgroup]/tau; d[isgroup] <- d[isgroup] / taugroup
        D <- diag(p); diag(D) <- d
        Sinv <- Vinv + D
        #Sinv <- Vinv + diag(p)/tau
      }
      S <- solve(Sinv)
      m <- S %*% Vinv %*% thhat
      ans <- rnlp(m=as.vector(m), V=S, priorCoef=priorCoef, priorGroup=priorGroup, niter=niter, burnin=burnin, thinning=thinning)
      if (is.null(colnames(x))) colnames(ans) <- paste('beta',1:ncol(x),sep='') else colnames(ans) <- colnames(x)
    }
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
    ans <- rnorm_truncMultCI(as.integer(n),as.double(lower),as.double(upper),as.double(m),as.double(sd))
    #ans <- .Call("rnorm_truncMultCI",as.integer(n),as.double(lower),as.double(upper),as.double(m),as.double(sd))
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
  ans <- rtmvnormCI(as.integer(n), as.double(m), as.double(Sigma), as.double(lower), as.double(upper), as.integer(within), method)
  #ans <- .Call("rtmvnormCI",as.integer(n), as.double(m), as.double(Sigma), as.double(lower), as.double(upper), as.integer(within), method)
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
  ans= rtmvnormProdCI(as.integer(n), as.double(m), as.double(Sigma), as.integer(k), as.double(lower), upper, lowtrunc, uptrunc, as.integer(burnin));
  #ans= .Call("rtmvnormProdCI",as.integer(n), as.double(m), as.double(Sigma), as.integer(k), as.double(lower), upper, lowtrunc, uptrunc, as.integer(burnin));
  matrix(ans,nrow=n)
}


