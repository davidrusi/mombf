##########################################################################################
##
## BAYESIAN MODEL SELECTION & AVERAGING UNDER DIAGONAL & BLOCK-DIAGONAL X'X
##
##########################################################################################


postModeOrtho <- function(y, x, priorCoef=momprior(tau=0.348), priorDelta=modelbbprior(1,1),priorVar=igprior(0.01,0.01), bma=FALSE, maxvars=100) {
    #Find posterior mode for linear regression under orthogonal X'X
    # - y: outcome
    # - x: predictors
    # - priorCoef: prior on coefficients
    # - priorDelta: prior on model space
    # - priorVar: prior on residual variance
    # - bma: set to TRUE to obtain Bayesian model averaging parameter estimates
    # - integrateMethod: if exactpp=TRUE this is the numerical integration method to obtain logpy. Set to 'CSR' to compound Simpson's rule based on an adaptive grid defined by enumerating highly probable models, to 'AQ' to use adaptive quadrature as implemented in R function 'integrate'
    # - maxvars: the search is restricted to models with up to maxvars variables (note: posterior model prob and BMA are valid regardless of maxvars)
    # Output
    # -
    integrateMethod <- 'AQ'; exactpp <- TRUE #exactpp=FALSE means we only find posterior mode (for Zellner's prior, for pMOM this doesn't work)
    n <- length(y); p <- ncol(x)
    if (priorDelta@priorDistr=='binomial' & ('p' %in% names(priorDelta@priorPars))) {
        rho <- priorDelta@priorPars['p']
        priorModel <- function(nvar) nvar*log(rho) + (p-nvar)*log(1-rho)
        #priorModel <- function(nvar) dbinom(nvar, size=ncol(x),prob=rho,log=TRUE)
    } else if (priorDelta@priorDistr=='binomial' & !('p' %in% names(priorDelta@priorPars))) {
        alpha=priorDelta@priorPars['alpha.p']; beta=priorDelta@priorPars['beta.p']
        priorModel <- function(nvar) lbeta(nvar + alpha, p - nvar + beta) - lbeta(alpha, beta)
        rho <- alpha/(alpha+beta)
    } else if (priorDelta@priorDistr=='uniform') {
        rho <- 0.5
        priorModel <- function(nvar) return(-p*log(2))
    } else { stop("Prior on model space not recognized. Use modelbbprior(), modelunifprior() or modelbinomprior()") }
    if (priorCoef@priorDistr == 'zellner') {
        g <- as.double(priorCoef@priorPars['tau'])
    } else if (priorCoef@priorDistr == 'pMOM') {
        tau <- as.double(priorCoef@priorPars['tau'])
        g <- n * tau
    } else {
        stop("priorCoef must be pMOM or Zellner. Use momprior() or zellnerprior()")
    }
    if (priorVar@priorDistr != 'invgamma') stop('priorVar must be an inverse gamma. Use igprior()')
    a.phi <- as.double(priorVar@priorPars['alpha'])
    l.phi <- as.double(priorVar@priorPars['lambda'])
    if (exactpp & !(integrateMethod %in% c('CSR','AQ'))) {
        stop("integrateMethod should be 'CSR' for compound Simpson's rule or 'AQ' for adaptive quadrature")
    }
    #Pre-compute useful quantities
    sumy2 <- l.phi + sum(y^2)
    apos <- a.phi+n
    shrinkage <- g/(1+g)
    Xty <- t(x) %*% matrix(y,ncol=1); XtX <- colSums(x^2)
    u <- Xty^2 / XtX
    #Avoid numerical overflow (also cases where X'X is not diagonal)
    sumu <- sum(u)
    if (sumu > sumy2-l.phi) u <- u * ((sumy2-l.phi)/sumu)
    #Consider sequence of models
    o <- order(u,decreasing=TRUE)
    uord <- u[o]
    cumsumu <- cumsum(uord)
    gamma <- rep(FALSE,ncol(x))
    modelid <- character(length(o)+1)
    modelusum <- c(0,cumsumu)
    iseq <- 0:length(o)
    lpos <- sumy2 - shrinkage * c(0,cumsumu)
    pp <- priorModel(iseq) - 0.5*iseq*log(1+g) - 0.5*apos * log(lpos)
    for (i in 1:length(o)) { modelid[i+1] <- paste(o[1:i],collapse=',') }
    maxpp <- max(pp)
    variableids <- strsplit(as.character(modelid),split=',')
    nvars <- sapply(variableids,length)
    sel <- nvars<=maxvars; modelid <- modelid[sel]; variableids <- variableids[sel]; nvars <- nvars[sel]
    if (any(!sel)) { modelusum <- modelusum[1:(maxvars+1)]; lpos <- lpos[1:(maxvars+1)] }
    if (exactpp) {
      qnull <- 1/qgamma(.001,shape=.5*apos[1],.5*lpos[1])
      qfull <- 1/qgamma(.999,shape=.5*apos[length(apos)],.5*lpos[length(lpos)])
      phiseq <- c(qnull,.5*unique(lpos) / (.5*apos + 1),qfull)
      #Refine grid for phi by adding further promising models
      goodmodelsizes <- which((maxpp - pp[-1]) < log(1000))
      if (length(goodmodelsizes)>0) {
          goodmodelsizes <- 1:(goodmodelsizes[length(goodmodelsizes)])
          goodmodelsizes <- goodmodelsizes[goodmodelsizes <= maxvars]
          phiseqextra <- modelidextra <- modelusumextra <- ppextra <- vector("list",length(goodmodelsizes))
          for (i in 1:length(goodmodelsizes)) {
              nvars <- goodmodelsizes[i]
              maxmodels <- lchoose(ncol(x),nvars)
              if (maxmodels>0) {
                  if (nvars == 1) {
                      #idx <- list(2,3,4)
                      idx <- list(2,3,4,5,6)
                  } else {
                      if (nvars>2) idx <- 1:(nvars-2) else idx <- integer(0)
                      #idx <- list(c(idx,nvars-1,nvars+1), c(idx,nvars-1,nvars+2), c(idx,nvars,nvars+1))
                      idx <- list(c(idx,nvars-1,nvars+1), c(idx,nvars-1,nvars+2), c(idx,nvars,nvars+1), c(idx,nvars-1,nvars+3), c(idx,nvars,nvars+2))
                  }
                  idx <- idx[1:ifelse(maxmodels>=log(6),5,exp(maxmodels)-1)]  #if <5 extra models available, just take those
                  modelidextra[[i]] <- sapply(idx,function(z) paste(o[z],collapse=','))
                  modelusumextra[[i]] <- sapply(idx, function(z) sum(uord[z]))
                  lposextra <- sumy2 - shrinkage * modelusumextra[[i]]
                  phiseqextra[[i]] <- .5* lposextra / (.5*apos + 1)
                  ppextra[[i]] <- priorModel(nvars) - 0.5*nvars*log(1+g) - 0.5*apos * log(lposextra)
              }
          }
          phiseq <- c(phiseq,unlist(phiseqextra))
          pp <- c(pp,unlist(ppextra))
          modelid <- c(modelid,unlist(modelidextra))
          modelusum <- c(modelusum,unlist(modelusumextra))
          variableidsextra <- strsplit(as.character(modelidextra),split=',')
          variableids <- c(variableids,variableidsextra)
          nvars <- c(nvars,sapply(variableidsextra,length))
      }
      #Evaluate marginal of phi at selected grid points
      phiseq <- phiseq[order(phiseq,decreasing=TRUE)]
      if (priorCoef@priorDistr == 'zellner') {
          phiseqpost <- jointPhiyZellner(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
      } else {  #pMOM
          phiseqpost <- jointPhiyMOM(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
      }
      #Refine grid where target increases > tolf % its max value
      maxphipost <- max(phiseqpost)
      tolf <- 0.01
      sel <- which(abs(diff(exp(phiseqpost-maxphipost))) > tolf)
      while (length(sel)>0) {
          phiseqnew <- unlist(lapply(sel, function(z) seq(phiseq[z],phiseq[z+1],length=5)))
          if (priorCoef@priorDistr == 'zellner') {
              phiseqpostnew <- jointPhiyZellner(phiseqnew, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
          } else {  #pMOM
              phiseqpostnew <- jointPhiyMOM(phiseqnew, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
          }
          phiseq <- c(phiseq,phiseqnew)
          phiseqpost <- c(phiseqpost,phiseqpostnew)
          ophiseq <- order(phiseq,decreasing=TRUE)
          phiseq <- phiseq[ophiseq]; phiseqpost <- phiseqpost[ophiseq]
          sel <- which(abs(diff(exp(phiseqpost-maxphipost))) > tolf)
      }
      #Compute marginal p(y)
      phiseqpost <- exp(phiseqpost-maxphipost)
      if (integrateMethod=='CSR') {
          myintegral <- int.simpson2(phiseq,phiseqpost,equi=FALSE,method="CSR")
      } else {
          if (priorCoef@priorDistr=='zellner') {
            f2int <- function(phi) { exp(jointPhiyZellner(phi,sumy2=sumy2,apos=apos,shrinkage=shrinkage,u=u,g=g,rho=rho,logscale=TRUE) - maxphipost) }
          } else {
            f2int <- function(phi) { exp(jointPhiyMOM(phi,sumy2=sumy2,apos=apos,shrinkage=shrinkage,u=u,g=g,rho=rho,logscale=TRUE) - maxphipost) }
          }
          myintegral <- integrate(f2int,phiseq[length(phiseq)],Inf)$value
      }
      phiseqpost <- phiseqpost / myintegral  # p(phi|y)
      logpy <- log(myintegral) + maxphipost  #log p(y)
      #Format
      phi <- data.frame(phi=phiseq,phipostprob=phiseqpost)
      phi <- phi[order(phi[,1],decreasing=TRUE),]
      #Exact posterior probabilities
      modelidx <- sapply(strsplit(modelid,split=','),as.numeric)
      nvars <- sapply(modelidx,length)
      if (priorCoef@priorDistr=='zellner') {
          pp <- priorModel(nvars) - 0.5*nvars*log(g+1) - logpy + lgamma(0.5*apos) - 0.5*apos * log((sumy2 - shrinkage * modelusum)/2)
          pp <- exp(pp)
      } else {
          pp <- double(length(modelid))
          w <- phi[,2]/sum(phi[,2])
          for (i in 1:length(modelid)) {
              ppcond <- sapply(phi[,1],jointppMomKnown,sel=modelidx[[i]],u=u,g=g,shrinkage=shrinkage,rho=rho,logscale=FALSE)
              pp[i] <- sum(ppcond * w)
          }
      }
    } else {
      phi <- NA
      logpy <- NA
      pp <- exp(pp - maxpp)
      pp <- pp/sum(pp)
    }
    #Remove unnecessary grid points
    dd <- abs(diff(phi[,2]))/max(phi[,2])
    ddsum <- 0; sel <- rep(TRUE,nrow(phi))
    for (i in 2:(nrow(phi)-1)) {
        ddsum <- ddsum + dd[i-1]
        if (ddsum < 0.0001) { sel[i] <- FALSE } else { ddsum <- 0 }
    }
    phi <- phi[sel,]
    #
    models <- data.frame(modelid=modelid,pp=pp)
    models <- models[order(pp,decreasing=TRUE),]
    ans <- list(models=models, phi=phi, logpy=logpy)
    if (bma) {
        margpp <- margpp2 <- double(p)
        postphigrid <- phi[,2]/sum(phi[,2])
        if (priorCoef@priorDistr=='zellner') {
            m <- shrinkage * Xty / XtX
            for (i in 1:length(u)) {
                margppcond <- margppZellnerKnown(phi=phi[,1],u=u[i],g=g,shrinkage=shrinkage,rho=rho,logscale=FALSE)
                margpp[i] <- sum(margppcond * postphigrid)
                #margpp[i] <- int.simpson2(phi[,1], margppcond * phi[,2], equi=FALSE, method="CSR")
            }
            pm <- m * margpp
        } else {
            m <- shrinkage * Xty / XtX
            pm <- double(p)
            for (i in 1:length(u)) {
                mi <- m[i] * (1 + 2 / (1+shrinkage*u[i]/phi[,1]))
                margppcond <- margppMomKnown(phi=phi[,1],u=u[i],g=g,shrinkage=shrinkage,rho=rho,logscale=FALSE)
                margpp[i] <- sum(margppcond * postphigrid)
                pm[i] <- sum(mi * margppcond * postphigrid)
                #margpp[i] <- int.simpson2(phi[,1], margppcond * phi[,2], equi=FALSE, method="CSR")
                #pm[i] <- int.simpson2(phi[,1], mi * margppcond * phi[,2], equi=FALSE, method="CSR")
            }
        }
        ans$bma <- data.frame(margpp=margpp,coef=pm)
    }
    return(ans)
}



## Auxiliary functions ##
#########################

#Marginal inclusion probability for a single variable conditional on a vector of phi values
# Input
# - phi: vector of phi values (residual variance)
# - u: u statistic for a single variable
# - g: prior dispersion parameter
# - shrinkage: g/(1+g)
# - rho: prior inclusion probability
# - logscale: set to TRUE to obtain log inclusion probability
# Output: vector of marginal inclusion probabilities conditional on the given phi's
margppZellnerKnown <- function(phi,u,g,shrinkage=g/(1+g),rho,logscale=TRUE) {
  ans <- -log(1 + sqrt(1+g) * exp(-.5 * shrinkage*u/phi) * (1-rho)/rho)
  if (!logscale) ans <- exp(ans)
  return(ans)
}

margppMomKnown <- function(phi,u,g,shrinkage=g/(1+g),rho,logscale=TRUE) {
  uu <- shrinkage*u/phi
  ans <- -log(1 + (1+g)^(1.5) * exp(-.5 * uu) * (1-rho)/(rho * (1+uu)))
  if (!logscale) ans <- exp(ans)
  return(ans)
}

#Joint model probability conditional on a single phi value
#Input
# - phi: single value for residual variance
# - sel: model indicator, i.e. vector indicating indexes of active variables under current model
# - u: u-scores for all variables (i.e. u[sel] selects those for active variables)
# - g, shrinkage, rho, logscale: as in margppMomKnown
#Output: posterior probability of the model conditional on phi
jointppMomKnown <- function(phi,sel,u,g,shrinkage=g/(1+g),rho,logscale=TRUE) {
    uu <- shrinkage*u/phi
    odds01 <- (1+g)^(1.5) * exp(-.5 * uu) * (1-rho)/(rho * (1+uu))
    if (length(sel)>0 & length(sel)<length(u)) {
        ans <- -sum(log(1+odds01[sel])) -sum(log(1+1/odds01[-sel]))
    } else if (length(sel)==0) {
        ans <- -sum(log(1+1/odds01))
    } else {
        ans <- -sum(log(1+odds01))
    }
    if (!logscale) ans <- exp(ans)
    return(ans)
}


#Evaluate quantity proportional to joint of (phi,y) on a grid
jointPhiyZellner <- function(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE) {
    #Old code: ok but gave overflow
    #fseqnew <- do.call(cbind,lapply(phiseq, function(phi) exp(0.5 * shrinkage * u / phi)))
    #ans <- -0.5*sumy2/phiseq - (.5*apos+1)*log(phiseq) + colSums(log(1 + rho*(fseqnew/sqrt(1+g) -1)))
    #ans <- -0.5*sumy2/phiseq - (.5*apos+1)*log(phiseq) + colSums(log(1-rho + rho*fseqnew/sqrt(1+g)))
    #New code
    ans <- double(length(phiseq))
    ct <- -0.5*sumy2/length(u)
    for (i in 1:length(phiseq)) {
        fseqnew <- exp((ct + 0.5 * shrinkage * u) / phiseq[i])
        ans[i] <- -(.5*apos+1)*log(phiseq[i]) + sum(log((1-rho)*exp(ct/phiseq[i]) + rho*fseqnew/sqrt(1+g)))
    }
    ans[is.infinite(ans)] <- -Inf  #avoid numerical overflow
    if (!logscale) ans <- exp(ans)
    return(ans)
}

jointPhiyMOM <- function(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE) {
    fseqnew <- do.call(cbind,lapply(phiseq, function(phi) exp(0.5 * shrinkage * u / phi + log(1 + shrinkage * u / phi) - log(1+g))))
    ans <- -0.5*sumy2/phiseq - (.5*apos+1)*log(phiseq) + colSums(log(1 + rho*(fseqnew/sqrt(1+g) -1)))
    ans[is.infinite(ans)] <- -Inf  #avoid numerical overflow
    if (!logscale) ans <- exp(ans)
    return(ans)
}


#Function int.simpson2 copied from R package fda.usc
# x: grid of x values
# y: f(x)
# equi: set to TRUE if grid has equally spaced points
# method: "TRAPZ" for trapezoidal rule, "CSR" for composite Simpson's rule, "ESR" for extended Simpson's rule
int.simpson2 <- function (x, y, equi=FALSE, method="CSR") {
    n = length(x)
    ny = length(y)
    if (n != ny) stop("Different length in the input data")
    if (n == 2 || ny == 2) method = "TRAPZ"
    out <- switch(method, TRAPZ = {
        idx = 2:length(x)
        value <- as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
    }, CSR = {
        if (!equi) {
            n = 2 * n - 1
            app = approx(x, y, n = n)
            x = app$x
            y = app$y
        }
        h = (max(x) - min(x))/(n - 1)
        value = (h/3) * (y[n] + y[1] + 2 * sum(y[2 * (1:((n-1)/2)) + 1]) + 4 * sum(y[2 * (1:((n - 1)/2))]))
    }, ESR = {
        if (!equi) {
            n = 2 * n - 1
            app = approx(x, y, n = n)
            x = app$x
            y = app$y
        }
        h = (max(x) - min(x))/(n - 1)
        if (n <= 4) stop("This method needs n>4")
        value = 17 * (y[1] + y[n]) + 59 * (y[2] + y[n - 1]) + 43 * (y[3] + y[n - 2]) + 49 * (y[4] + y[n - 3])
        value = value + 48 * sum(y[5:(n - 4)])
        value = (h/48) * value
    })
    return(out)
}
