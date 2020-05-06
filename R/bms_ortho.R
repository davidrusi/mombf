##########################################################################################
##
## BAYESIAN MODEL SELECTION & AVERAGING UNDER DIAGONAL & BLOCK-DIAGONAL X'X
##
##########################################################################################


postModeOrtho <- function(y, x, priorCoef=momprior(tau=0.348), priorDelta=modelbbprior(1,1),priorVar=igprior(0.01,0.01), bma=FALSE, includeModels, maxvars=100) {
    #Find posterior mode for linear regression under orthogonal X'X
    # - y: outcome
    # - x: predictors
    # - priorCoef: prior on coefficients
    # - priorDelta: prior on model space
    # - priorVar: prior on residual variance
    # - bma: set to TRUE to obtain Bayesian model averaging parameter estimates
    # - includeModels: list of models that should always be included when computing posterior model probabilities. Each element should be a vector z where x[,z] selects the required variables
    # - integrateMethod: if exactpp=TRUE this is the numerical integration method to obtain logpy. Set to 'CSR' to compound Simpson's rule based on an adaptive grid defined by enumerating highly probable models, to 'AQ' to use adaptive quadrature as implemented in R function 'integrate'
    # - maxvars: the search is restricted to models with up to maxvars variables (note: posterior model prob and BMA are valid regardless of maxvars)
    # Output
    # -
    integrateMethod <- 'AQ'; exactpp <- TRUE #exactpp=FALSE means we only find posterior mode (for Zellner's prior, for pMOM this doesn't work)
    n <- length(y); p <- ncol(x)
    if (priorDelta@priorDistr=='binomial' & ('p' %in% names(priorDelta@priorPars))) {
        rho <- priorDelta@priorPars[['p']]
        if (length(rho)>1) stop("postModeOrtho not implemented for vector p in modelbinomprior")
        priorModel <- function(nvar) nvar*log(rho) + (p-nvar)*log(1-rho)
        #priorModel <- function(nvar) dbinom(nvar, size=ncol(x),prob=rho,log=TRUE)
    } else if (priorDelta@priorDistr=='binomial' & !('p' %in% names(priorDelta@priorPars))) {
        alpha=priorDelta@priorPars['alpha.p']; beta=priorDelta@priorPars['beta.p']
        priorModel <- function(nvar) lbeta(nvar + alpha, p - nvar + beta) - lbeta(alpha, beta)
        rho <- alpha/(alpha+beta)
    } else if (priorDelta@priorDistr=='uniform') {
        priorModel <- function(nvar) -p*log(2)
        rho <- 1/2
    } else if (priorDelta@priorDistr=='complexity') {
        priorc= priorDelta@priorPars['c']
        priornorm= log(1-1/p^(priorc*(p+1))) - log(1-1/p^priorc)
        priorModel <- function(nvar) lbeta(nvar + 1, p - nvar + 1) - (priorc*nvar)*log(p) - priornorm
        lseq= 1:(p+1)
        rho= sapply(lseq, function(l) -(priorc*l)*log(p) + lchoose(p-1,l-1) - lchoose(p,l) - priornorm )
        rho= sum(exp(rho[!is.nan(rho)]))
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
    #gamma <- rep(FALSE,ncol(x))
    modelid <- character(length(o)+1)
    modelusum <- c(0,cumsumu)
    iseq <- 0:length(o)
    lpos <- sumy2 - shrinkage * c(0,cumsumu)
    pp <- priorModel(iseq) - 0.5*iseq*log(1+g) - 0.5*apos * log(lpos)
    for (i in 1:length(o)) { modelid[i+1] <- paste((o[1:i])[order(o[1:i])],collapse=',') }
    #Compute lpos, modelid, modelusum for includeModels
    if (!missing(includeModels)) {
        iextra= modelusumextra= double(length(includeModels))
        modelidextra= character(length(includeModels))
        for (i in 1:length(includeModels)) {
            if (is.numeric(includeModels[[i]])) {
                modelidextra[i]= paste(includeModels[[i]][order(includeModels[[i]])],collapse=',')
                iextra[i]= length(includeModels[[i]])
            } else if (is.logical(includeModels[[i]])) {
                modelidextra[i]= paste(which(includeModels[[i]]),collapse=',')
                iextra[i]= sum(includeModels[[i]])
            } else { stop("includeModels must be a list, each element either a logical or numeric vector indicating the variables to be included") }
            modelusumextra[i]= sum(u[includeModels[[i]]])
        }
        sel= !(modelidextra %in% modelid)
        if (any(sel)) {
            modelid= c(modelid,modelidextra[sel])
            modelusum= c(modelusum,modelusumextra[sel])
            lposextra= sumy2 - shrinkage * modelusumextra[sel]
            lpos= c(lpos, lposextra)
            pp= c(pp, priorModel(iextra[sel]) - 0.5*iextra[sel]*log(1+g) - 0.5*apos * log(lposextra[sel]))
        }
    }
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
          if (length(goodmodelsizes)!=1 | goodmodelsizes[1]!=ncol(x)) {  #skip if only 1 model with ncol(x) variables remains
              #goodmodelsizes <- 1:(goodmodelsizes[length(goodmodelsizes)])
              goodmodelsizes <- goodmodelsizes[goodmodelsizes <= min(maxvars,ncol(x)-1)]
              phiseqextra <- modelidextra <- modelusumextra <- ppextra <- vector("list",length(goodmodelsizes))
              for (i in 1:length(goodmodelsizes)) {
                  nvars <- goodmodelsizes[i]
                  maxmodels <- lchoose(ncol(x),nvars)
                  if (maxmodels>0) {
                      if (nvars == 1) {
                          idx <- list(2,3,4,5,6) #idx <- list(2,3,4)
                          idx <- idx[idx<=ncol(x)]
                      } else {
                          if (nvars>2) idx <- 1:(nvars-2) else idx <- integer(0)
                          if (nvars+3<=ncol(x)) {
                              idx <- list(c(idx,nvars-1,nvars+1), c(idx,nvars-1,nvars+2), c(idx,nvars,nvars+1), c(idx,nvars-1,nvars+3), c(idx,nvars,nvars+2))
                          } else if (nvars+2<=ncol(x)) {
                              idx <- list(c(idx,nvars-1,nvars+1), c(idx,nvars-1,nvars+2), c(idx,nvars,nvars+1))
                          } else {
                              idx <- list(c(idx,nvars-1,nvars+1), c(idx,nvars,nvars+1))
                          }
                      }
                      idx <- idx[1:min(length(idx),ifelse(maxmodels>=log(6),5,exp(maxmodels)-1))]  #if <5 extra models available, just take those
                      modelidextra[[i]] <- sapply(idx,function(z) paste(o[z][order(o[z])],collapse=','))
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
              #variableidsextra <- strsplit(as.character(modelidextra),split=',')
              variableidsextra <- strsplit(unlist(modelidextra),split=',')
              variableids <- c(variableids,variableidsextra)
              nvars <- c(nvars,sapply(variableidsextra,length))
          }
      }
      #Evaluate marginal of phi at selected grid points
      phiseq <- phiseq[order(phiseq,decreasing=TRUE)]
      if (priorCoef@priorDistr == 'zellner') {
          phiseqpost <- jointPhiyZellnerOrtho(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
      } else {  #pMOM
          phiseqpost <- jointPhiyMOMOrtho(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
      }
      #Avoid phi values associated to 0 posterior density
      zeropost <- match(-Inf,phiseqpost)
      if (!is.na(zeropost)) {
          zerophi <- phiseq[zeropost]
          phiseq <- phiseq[1:max(1,zeropost-1)]
          phiseqpost <- phiseqpost[1:max(1,zeropost-1)]
      }
      #Refine grid where target increases > tolf % its max value
      maxphipost <- max(phiseqpost)
      tolf <- 0.01
      sel <- which(abs(diff(exp(phiseqpost-maxphipost))) > tolf)
      while (length(sel)>0) {
          phiseqnew <- unlist(lapply(sel, function(z) seq(phiseq[z],phiseq[z+1],length=5)))
          if (priorCoef@priorDistr == 'zellner') {
              phiseqpostnew <- jointPhiyZellnerOrtho(phiseqnew, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
          } else {  #pMOM
              phiseqpostnew <- jointPhiyMOMOrtho(phiseqnew, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE)
          }
          maxphipost <- max(maxphipost, max(phiseqpostnew))
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
            f2int <- function(phi) { exp(jointPhiyZellnerOrtho(phi,sumy2=sumy2,apos=apos,shrinkage=shrinkage,u=u,g=g,rho=rho,logscale=TRUE) - maxphipost) }
          } else {
            f2int <- function(phi) { exp(jointPhiyMOMOrtho(phi,sumy2=sumy2,apos=apos,shrinkage=shrinkage,u=u,g=g,rho=rho,logscale=TRUE) - maxphipost) }
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
          #Old version. It was numerically unstable
          #pp <- priorModel(nvars) - 0.5*nvars*log(g+1) - logpy + lgamma(0.5*apos) - 0.5*apos * log((sumy2 - shrinkage * modelusum)/2)
          #pp <- exp(pp)
          pp <- double(length(modelid))
          w <- phi[,2]/sum(phi[,2])
          for (i in 1:length(modelid)) {
              ppcond <- sapply(phi[,1],jointppZellnerKnownOrtho,sel=modelidx[[i]],u=u,g=g,shrinkage=shrinkage,rho=rho,logscale=FALSE)
              pp[i] <- sum(ppcond * w)
          }
      } else {
          pp <- double(length(modelid))
          w <- phi[,2]/sum(phi[,2])
          for (i in 1:length(modelid)) {
              ppcond <- sapply(phi[,1],jointppMomKnownOrtho,sel=modelidx[[i]],u=u,g=g,shrinkage=shrinkage,rho=rho,logscale=FALSE)
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
    if (sum(models$pp)>1) {
        warning("Posterior model probabilities are unstable, this often signals that t(x) %*% x is not diagonal")
        models$pp <- models$pp/sum(models$pp)
    }
    ans <- list(models=models, phi=phi, logpy=logpy)
    if (bma) {
        margpp <- double(p)
        postphigrid <- phi[,2]/sum(phi[,2])
        if (priorCoef@priorDistr=='zellner') {
            m <- shrinkage * Xty / XtX
            for (i in 1:length(u)) {
                margppcond <- margppZellnerKnownOrtho(phi=phi[,1],u=u[i],g=g,shrinkage=shrinkage,rho=rho,logscale=FALSE)
                margpp[i] <- sum(margppcond * postphigrid)
                #margpp[i] <- int.simpson2(phi[,1], margppcond * phi[,2], equi=FALSE, method="CSR")
            }
            pm <- m * margpp
        } else {
            m <- shrinkage * Xty / XtX
            pm <- double(p)
            for (i in 1:length(u)) {
                mi <- m[i] * (1 + 2 / (1+shrinkage*u[i]/phi[,1]))
                margppcond <- margppMomKnownOrtho(phi=phi[,1],u=u[i],g=g,shrinkage=shrinkage,rho=rho,logscale=FALSE)
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




postModeBlockDiag <- function(y, x, blocks, priorCoef=zellnerprior(tau=nrow(x)), priorDelta=modelbinomprior(p=1/ncol(x)), priorVar=igprior(0.01,0.01), bma=FALSE, maxvars=100, momcoef) {
    #Find posterior mode for linear regression under block-diagonal X'X
    # - y: outcome
    # - x: predictors
    # - blocks: factor or integer of length ncol(x) indicating the block that each column in x belongs to
    # - priorCoef: prior on coefficients
    # - priorDelta: prior on model space
    # - priorVar: prior on residual variance
    # - bma: set to TRUE to obtain Bayesian model averaging parameter estimates
    # - integrateMethod: if exactpp=TRUE this is the numerical integration method to obtain logpy. Set to 'CSR' to compound Simpson's rule based on an adaptive grid defined by enumerating highly probable models, to 'AQ' to use adaptive quadrature as implemented in R function 'integrate'
    # - maxvars: the search is restricted to models with up to maxvars variables (note: posterior model prob and BMA are valid regardless of maxvars)
    # - momcoef: optional argument containing pre-computed coefficients needed to obtain the marginal likelihood under the pMOM prior. A first call to postModeBlockDiag returns these coefficients, thus this argument is useful to speed up successive calls.
    # Output
    # -
    integrateMethod <- 'AQ'; exactpp <- TRUE #exactpp=FALSE means we only find posterior mode (for Zellner's prior, for pMOM this doesn't work)
    n <- length(y); p <- ncol(x)
    if (length(blocks) != p) stop('length(blocks) should be equal to ncol(x)')
    blocks <- as.integer(factor(blocks))
    nblocks <- max(blocks)
    blocksize <- table(blocks)
    nn <- 1:ncol(x)
    varidx <- lapply(as.integer(names(blocksize)), function(k) nn[blocks==k])
    if (priorDelta@priorDistr=='binomial' & ('p' %in% names(priorDelta@priorPars))) {
        rho <- priorDelta@priorPars[['p']]
        if (length(rho)>1) stop("postModeBlockDiag not implemented for vector p in modelbinomprior")
        if (any(rho<0) | any(rho>1)) stop("Specified prior inclusion probability outside [0,1]. Check priorDelta")
        priorModel <- function(nvar) nvar*log(rho) + (p-nvar)*log(1-rho)
        priorModelBlock <- function(nvar,blocksize) nvar*log(rho) + (blocksize-nvar)*log(1-rho)
    } else if (priorDelta@priorDistr=='binomial' & !('p' %in% names(priorDelta@priorPars))) {
        stop("Beta-Binomial prior not currently implemented")
        #alpha=priorDelta@priorPars['alpha.p']; beta=priorDelta@priorPars['beta.p']
        #priorModel <- function(nvar) lbeta(nvar + alpha, p - nvar + beta) - lbeta(alpha, beta)
    } else if (priorDelta@priorDistr=='uniform') {
        rho <- 0.5
        priorModel <- function(nvar) rep(-p*log(2),length(nvar))
        priorModelBlock <- function(nvar,blocksize) rep(-blocksize*log(2),length(nvar))
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
    Xty <- t(x) %*% matrix(y,ncol=1)
    XtX <- lapply(1:nblocks,function(i) { sel <- which(blocks==i); t(x[,sel]) %*% x[,sel] })
    #Compute u-scores, least squares estimates and covariances within blocks
    ubyblocks <- blockuscores(blocks, blocksize=blocksize, y=y, Xty=Xty, XtX=XtX)
    u <- ubyblocks$u; umv <- ubyblocks$umv; ubest <- ubyblocks$ubest
    #Avoid overflow in cases where X'X not exactly block-diagonal
    sumumax <- sum(sapply(u, function(z) z$u[1]))
    if (sumumax > (sumy2-l.phi)) {
        u <- lapply(u, function(z) { z$u <- z$u * ((sumy2-l.phi)/sumumax); return(z) } )
        ubest <- lapply(ubest, function(z) { z$u <- z$u * ((sumy2-l.phi)/sumumax); return(z) } )
        for (i in 1:length(umv)) { umv[[i]] <- lapply(umv[[i]], function(z) { z$u <- z$u * ((sumy2-l.phi)/sumumax); return(z) } ) }
    }
    #Enumerate models
    models <- coolblock(ubest,g=g,priorModelBlock=priorModelBlock,varidx=varidx,maxvars=maxvars)
    models <- rbind(data.frame(nvars=0,u=0,modelid='',phi=Inf,m=NA,block=NA,u.lower=0,u.upper=0,stringsAsFactors=FALSE),models)
    lpos.upper <- sumy2 - shrinkage * models$u.upper
    pp.upper <- priorModel(models$nvars) - 0.5*models$nvars*log(1+g) - 0.5*apos * log(lpos.upper)
    mis <- is.na(models$u) & (pp.upper>max(pp.upper[!is.na(models$u)],na.rm=TRUE))
    if (any(mis)) {
        enumsizes <- models$nvar[mis & (pp.upper>max(pp.upper[!mis]))]
        for (i in 1:length(enumsizes)) {
            #Enumerate all models candidate to be best of variable enumsizes[i]
            candidateModels <- enumblock(enumsizes[i], blocksize=blocksize)
            #Find u-scores
            candidateU <- matrix(0, nrow=nrow(candidateModels), ncol=ncol(candidateModels))
            for (k in 1:ncol(candidateModels)) {
                sel <- candidateModels[,k]>0
                candidateU[sel,k] <- ubest[[k]][candidateModels[sel,k],'u']
            }
            candidateU <- rowSums(candidateU)
            sel <- which.max(candidateU)
            models$u[enumsizes[i]+1] <- models$u.lower[enumsizes[i]+1] <- models$u.upper[enumsizes[i]+1] <- candidateU[sel]
            #Figure out model id
            modelsel <- character(ncol(candidateModels))
            for (k in 1:ncol(candidateModels)) {
                if (candidateModels[sel,k]>0) {
                    modelsel[k] <- paste(varidx[[k]][as.integer(strsplit(ubest[[k]][candidateModels[sel,k],'modelid'], split=',')[[1]])],collapse=',')
                }
            }
            models$modelid[enumsizes[i]+1] <- paste(modelsel[modelsel!=""],collapse=',')
            pp.upper[enumsizes[i]+1] <- priorModel(enumsizes[i]) - 0.5*enumsizes[i]*log(1+g) - 0.5*apos * log(sumy2 - shrinkage * candidateU[sel])
        }
    }
    #Marginal posterior of phi
    phiseq <- models$phi[!is.na(models$phi)][-1]
    qnull <- 1/qgamma(.001,shape=.5*apos,.5*sumy2)
    lfull <- max(c(sumy2 - shrinkage * max(models$u,na.rm=TRUE),l.phi))
    qfull <- 1/qgamma(.999,shape=.5*apos,.5*lfull)
    if (qnull > phiseq[1]) phiseq <- c(qnull,phiseq)
    if (qfull < phiseq[length(phiseq)]) phiseq <- c(phiseq,qfull)
    #Evaluate marginal of phi at selected grid points
    if (priorCoef@priorDistr == 'zellner') {
        phiseqpost <- jointPhiyZellnerBlockDiag(phiseq,sumy2=sumy2,apos=apos,a.phi=a.phi,l.phi=l.phi,shrinkage=shrinkage,u=u,g=g,priorModelBlock=priorModelBlock,blocksize=blocksize,logscale=TRUE)
        maxphipost <- max(phiseqpost)
        while (phiseqpost[length(phiseqpost)]==maxphipost) {
            phiseq <- c(phiseq,phiseq[length(phiseq)]/2)
            phiseqpost <- c(phiseqpost,jointPhiyZellnerBlockDiag(phiseq[length(phiseq)],sumy2=sumy2,apos=apos,a.phi=a.phi,l.phi=l.phi,shrinkage=shrinkage,u=u,g=g,priorModelBlock=priorModelBlock,blocksize=blocksize,logscale=TRUE))
            maxphipost <- max(c(maxphipost,phiseqpost[length(phiseqpost)]))
        }
    } else {  #pMOM
        if (missing(momcoef)) {
            momcoef <- vector("list",length(umv))
            for (k in 1:length(momcoef)) momcoef[[k]] <- lapply(umv[[k]], momPolyCoef, shrinkage=shrinkage, g=g)
        }
        phiseqpost <- jointPhiyMOMBlockDiag(phiseq,sumy2=sumy2,apos=apos,a.phi=a.phi,l.phi=l.phi,shrinkage=shrinkage,u=u,umv=umv,momcoef=momcoef,g=g,priorModelBlock=priorModelBlock,blocksize=blocksize)
        maxphipost <- max(phiseqpost)
        while (phiseqpost[length(phiseqpost)]==maxphipost) {
            phiseq <- c(phiseq,phiseq[length(phiseq)]/2)
            phiseqpost <- c(phiseqpost,jointPhiyMOMBlockDiag(phiseq[length(phiseq)],sumy2=sumy2,apos=apos,a.phi=a.phi,l.phi=l.phi,shrinkage=shrinkage,u=u,umv=umv,momcoef=momcoef,g=g,priorModelBlock=priorModelBlock,blocksize=blocksize,logscale=TRUE))
            maxphipost <- max(c(maxphipost,phiseqpost[length(phiseqpost)]))
        }
    }
    #Refine grid where target increases > tolf % its max value
    tolf <- 0.01
    sel <- which(abs(diff(exp(phiseqpost-maxphipost))) > tolf)
    itrefine <- 0
    while ((length(sel)>0) && (itrefine<10)) {
        phiseqnew <- unlist(lapply(sel, function(z) { ans= seq(phiseq[z],phiseq[z+1],length=5); ans[c(-1,-length(ans))] } ))
        if (priorCoef@priorDistr == 'zellner') {
            phiseqpostnew <- jointPhiyZellnerBlockDiag(phiseqnew, sumy2, apos, a.phi, l.phi, shrinkage, u, g, priorModelBlock, blocksize, logscale=TRUE)
        } else {  #pMOM
            phiseqpostnew <- jointPhiyMOMBlockDiag(phiseqnew, sumy2, apos, a.phi, l.phi, shrinkage, u, umv, momcoef, g, priorModelBlock, blocksize, logscale=TRUE)
        }
        maxphipost <- max(maxphipost, max(phiseqpostnew))
        phiseq <- c(phiseq,phiseqnew)
        phiseqpost <- c(phiseqpost,phiseqpostnew)
        ophiseq <- order(phiseq,decreasing=TRUE)
        phiseq <- phiseq[ophiseq]; phiseqpost <- phiseqpost[ophiseq]
        sel <- which(abs(diff(exp(phiseqpost-maxphipost))) > tolf)
        itrefine <- itrefine+1
    }
    #Compute marginal p(y)
    zeropost <- which(is.infinite(phiseqpost))
    if (length(zeropost)>1) { phiseq <- phiseq[-zeropost[-1]]; phiseqpost <- phiseqpost[-zeropost[-1]] }
    phiseqpost <- exp(phiseqpost-maxphipost)
    if (integrateMethod=='CSR') {
        myintegral <- int.simpson2(phiseq[phiseqpost!=0],phiseqpost[phiseqpost!=0],equi=FALSE,method="CSR")
    } else {
        if (priorCoef@priorDistr=='zellner') {
            f2int <- function(phi) { exp(jointPhiyZellnerBlockDiag(phi,sumy2=sumy2,apos=apos,a.phi=a.phi,l.phi=l.phi,shrinkage=shrinkage,u=u,g=g,priorModelBlock=priorModelBlock,blocksize=blocksize,logscale=TRUE) - maxphipost) }
        } else {
            f2int <- function(phi) { exp(jointPhiyMOMBlockDiag(phi,sumy2=sumy2,apos=apos,a.phi=a.phi,l.phi=l.phi,shrinkage=shrinkage,u=u,umv=umv,momcoef=momcoef,g=g,priorModelBlock=priorModelBlock,blocksize=blocksize,logscale=TRUE) - maxphipost) }
        }
        myintegral <- integrate(f2int,phiseq[length(phiseq)],Inf)$value
    }
    phiseqpost <- phiseqpost / myintegral  # p(phi|y)
    logpy <- log(myintegral) + maxphipost  #log p(y)
    #Format
    phi <- data.frame(phi=phiseq,phipostprob=phiseqpost)
    phi <- phi[order(phi[,1],decreasing=TRUE),]
    #Remove unnecessary grid points
    dd <- abs(diff(phi[,2]))/max(phi[,2])
    ddsum <- 0; sel <- rep(TRUE,nrow(phi))
    for (i in 2:(nrow(phi)-1)) {
        ddsum <- ddsum + dd[i-1]
        if (ddsum < 0.0001) { sel[i] <- FALSE } else { ddsum <- 0 }
    }
    phi <- phi[sel,]
    #Exact posterior probabilities
    nvars <- models$nvars
    if (priorCoef@priorDistr=='zellner') {
        pp.upper <- priorModel(nvars) - 0.5*nvars*log(g+1) - logpy + lgamma(0.5*apos) - 0.5*apos * log((sumy2 - shrinkage*models$u.upper)/2)
        pp.upper <- exp(pp.upper)
        pp <- ifelse(is.na(models$u),NA,pp.upper)
    } else {
        modelidx <- sapply(strsplit(models$modelid,split=','), function(z) { ans= rep(FALSE,p); ans[as.numeric(z)] <- TRUE; return(ans) })
        pp <- double(nrow(models))
        w <- phi[,2]/sum(phi[,2])
        for (i in 1:length(pp)) {
            ppcond <- sapply(phi[,1],jointppMomKnownBlockDiag,sel=modelidx[,i],blocks=blocks,blocksize=blocksize,u=u,g=g,shrinkage=shrinkage,priorModelBlock=priorModelBlock,momcoef=momcoef,logscale=FALSE)
            pp[i] <- sum(ppcond * w)
        }
        pp.upper <- pp
    }
    models <- data.frame(models[,c('modelid','nvars')],pp=pp,pp.upper=pp.upper)
    ans <- list(models=models, phi=phi, logpy=logpy)
    #Coefficient estimates under each visited model
    modelidx <- sapply(strsplit(models$modelid,split=','), function(z) { ans= rep(FALSE,p); ans[as.numeric(z)] <- TRUE; return(ans) })
    pmcond <- matrix(0, nrow=nrow(models), ncol=p)
    for (i in 2:nrow(models)) {
        for (k in 1:length(umv)) {
            includeidx <- which(modelidx[blocks==k,i])
            if (length(includeidx)>0) {
                mincluded <- rep(0,blocksize[k])
                mincluded[includeidx] <- umv[[k]][[paste(includeidx,collapse=',')]]$m
                pmcond[i,blocks==k] <- mincluded
            }
        }
    }
    if (priorCoef@priorDistr=='zellner') pmcond <- shrinkage * pmcond
    ans$postmean.model <- data.frame(modelid=models$modelid,pmcond)
    #Marginal inclusion probabilities and BMA estimates
    if (bma) {
        margpp <- pm <- double(p)
        postphigrid <- phi[,2]/sum(phi[,2])
        if (priorCoef@priorDistr=='zellner') {
            for (k in 1:length(u)) {
                margppcond <- margppZellnerKnownBlockDiag(phi=phi[,1],u=u[[k]],umv=umv[[k]],g=g,shrinkage=shrinkage,priorModelBlock=priorModelBlock,logscale=FALSE)
                margpp[blocks==k] <- colSums(margppcond$margpp * postphigrid)
                pm[blocks==k] <- colSums(margppcond$coef * postphigrid)
            }
        } else {
            margpp <- pm <- double(p)
            for (k in 1:length(u)) {
                margppcond <- margppMOMKnownBlockDiag(phi=phi[,1],u=u[[k]],umv=umv[[k]],g=g,shrinkage=shrinkage,priorModelBlock=priorModelBlock,momcoef=momcoef[[k]],logscale=FALSE)
                margpp[blocks==k] <- colSums(margppcond$margpp * postphigrid)
                pm[blocks==k] <- colSums(margppcond$coef * postphigrid)
            }
        }
        ans$bma <- data.frame(margpp=margpp,coef=pm)
    }
    if (priorCoef@priorDistr=='pMOM') ans$momcoef <- momcoef
    return(ans)
}


blockuscores <- function(blocks, blocksize, y, Xty, XtX) {
    u <- umv <- ubest <- vector("list",length(blocksize))
    for (i in 1:length(u)) {
        Xtyblock <- Xty[blocks==i]
        modelsblock <- expand.grid(lapply(1:blocksize[i], function(z) c(FALSE,TRUE)))
        umv[[i]] <- apply(modelsblock[-1,,drop=FALSE], 1, function(sel) uscore(y,XtX[[i]][sel,sel],Xtyblock[sel]))
        ublock <- c(0,sapply(umv[[i]],'[[',1)) #null model has u-score=0
        nvarsblock <- rowSums(modelsblock)
        o <- order(nvarsblock,ublock,decreasing=c(TRUE,FALSE))
        modelid <- apply(modelsblock,1,function(z) paste(which(z),collapse=','))
        u[[i]] <- data.frame(modelid=modelid[o],u=ublock[o],nvars=nvarsblock[o],stringsAsFactors=FALSE)
        names(umv[[i]]) <- modelid[-1]
        ubest[[i]] <- do.call(rbind,by(u[[i]], INDICES=u[[i]][,'nvars'], FUN=function(z) z[1,]))[-1,] #best model of each size
    }
    return(list(u=u,umv=umv,ubest=ubest))
}


coolblock <- function(ubest, g, priorModelBlock, varidx, maxvars) {
    #Coolblock algorithm. Returns sequence of models with highest post prob under Zellner's prior conditional on any value for phi (residual variance)
    #Input
    # - ubest: list with as many elements as blocks. For each block it contains a data.frame with the best model (largest u-score) of each size
    # - g: Zellner prior dispersion parameter
    # - priorModelBlock: function taking a single integer argument k and returning the log-prior probability of including k variables within a block
    # - varidx: list with as many elements as blocks. For each block it indicates the indexes of the variables in that block
    # - maxvars: stop when >maxvars included into the model
    #Output: data.frame with info about best model of each size, containing the following columns
    # - nvars: number of variables
    # - u: u-score for best model of that size (may be missing for models not visited by coolblock)
    # - modelid: model identifier (comma-separated variable indexes for variables included in the model)
    # - phi: decreasing thresholds on residual variance. modelid[i] has highest posterior prob conditional on any residual variance <=phi[i] and >=phi[i+1]
    # - m: the posterior mode given phi[i] is obtained by setting active variables in block[i] to m[i]
    # - block: the posterior mode given phi[i] is obtained by setting active variables in block[i] to m[i]
    # - u.upper: upper bound for u (equals the exact u if this model size was visited by the algorithm)
    # - u.lower: lower bound for u (equals the exact u if this model size was visited by the algorithm)
    K <- length(ubest)
    shrinkage <- g/(1+g)
    log1g <- log(1+g)
    nvars <- sapply(ubest,nrow)
    #Store prior odds. The code takes advantage that for uniform/binom the odds don't depend on the block size
    #priorodds[[l]][j] has prior odds of l vars vs j vars
    priorodds <- lapply(1:max(nvars), function(l) { priorModelBlock(l,blocksize=max(nvars)) - sapply(0:(l-1), function(j) priorModelBlock(j,blocksize=max(nvars))) })
    #Step 1. Within-blocks optimal configuration
    #
    phiseq <- r <- m <- vector("list",K)
    for (k in 1:length(phiseq)) {
        #Obtain all pairwise r scores
        r[[k]] <- matrix(NA,nrow=nvars[k]+1,ncol=nvars[k]+1)
        rownames(r[[k]]) <- colnames(r[[k]]) <- 0:nvars[k]
        r[[k]][-1,1] <- r[[k]][1,-1] <- shrinkage * ubest[[k]][,'u'] / (ubest[[k]]$nvars * log1g - 2 * sapply(priorodds[ubest[[k]]$nvars],'[[',1))
        #r[[k]][-1,1] <- r[[k]][1,-1] <- shrinkage * ubest[[k]][,'u'] / (ubest[[k]]$nvars * log1g - 2 * priorodds[[nvars[k]]][ubest[[k]]$nvars])
        if (nvars[k]>=2) {
          for (l in 2:nvars[k]) {
              for (j in 1:(l-1)) {
                  r[[k]][l+1,j+1] <- r[[k]][j+1,l+1] <- shrinkage * (ubest[[k]][l,'u'] - ubest[[k]][j,'u']) / ((l-j) * log1g - 2 * priorodds[[l]][j+1])
              }
          }
        }
        phiseq[[k]] <- double(nvars[k])
        m[[k]] <- character(nvars[k])
        ll <- 1
        if (nvars[k]>=2) {
            for (l in 1:(nvars[k]-1)) {
                rmin <- min(r[[k]][l+1,1:l])
                rmax <- max(r[[k]][(l+2):(nvars[k]+1),l+1])
                if (rmax < rmin) {
                    phiseq[[k]][ll] <- rmin
                    m[[k]][ll] <- as.character(ubest[[k]][l,'modelid'])
                    ll <- ll+1
                }
            }
        }
        l <- nvars[k]  #treat last case separately (no larger models to compare with)
        phiseq[[k]][ll] <- min(r[[k]][l+1,1:l])
        m[[k]][ll] <- as.character(ubest[[k]][l,'modelid'])
        phiseq[[k]] <- phiseq[[k]][1:ll]
        m[[k]] <- m[[k]][1:ll]
    }
    #
    #Step 2. Combine blocks
    #
    phiseq <- data.frame(phi=unlist(phiseq),m=unlist(m),block=rep(1:length(phiseq),sapply(phiseq,length)),stringsAsFactors=FALSE)
    phiseq <- phiseq[order(phiseq[,1],decreasing=TRUE),]
    #
    #Step 3. Sequence of models
    #
    modelid <- character(nrow(phiseq))
    nvarsmodel <- integer(nrow(phiseq))
    useq <- double(nrow(phiseq))
    varidxWithinBlock <- strsplit(phiseq$m,split=',')
    curmodel <- curvars <- lapply(1:length(ubest), function(k) character(0))
    curu <- rep(0,length(ubest))
    i <- 1; stopvars <- FALSE
    while ((i <= nrow(phiseq) & (!stopvars))) {
        curblock <- phiseq$block[[i]]
        curmodel[[curblock]] <- phiseq$m[i]
        curvars[[curblock]] <- varidx[[curblock]][as.numeric(varidxWithinBlock[[i]])]
        nn <- unlist(curvars)
        modelid[i] <- paste(nn,collapse=',')
        nvarsmodel[i] <- length(nn)
        curu[curblock] <- ubest[[curblock]][length(varidxWithinBlock[[i]]),'u']
        useq[i] <- sum(curu)
        if (nvarsmodel[i]>maxvars) stopvars <- TRUE
        i <- i+1
    }
    if (!stopvars) {
        modelseq <- data.frame(nvars=nvarsmodel,u=useq,modelid=modelid,phiseq,stringsAsFactors=FALSE)
    } else {
        sel <- 1:(i-1)
        modelseq <- data.frame(nvars=nvarsmodel[sel],u=useq[sel],modelid=modelid[sel],phiseq[sel,],stringsAsFactors=FALSE)
    }
    rownames(modelseq) <- modelseq$nvars
    ans <- data.frame(nvars=1:max(modelseq$nvars),u=NA,modelid=NA,phi=NA,m=NA,block=NA)
    ans[modelseq$nvars,] <- modelseq
    #Best models of size 1 & 2 are always returned
    if (any(is.na(ans$u[1:2]))) onevar <- data.frame(block=1:length(ubest),do.call(rbind,lapply(ubest,function(z) z[1,])), stringAsFactors=FALSE)
    if (is.na(ans$u[1])) {
        sel <- which.max(onevar$u)
        ans$u[1] <- onevar$u[sel]
        ans$block[1] <- onevar$block[sel]
        ans$modelid[1] <- ans$m[1] <- varidx[[onevar$block[sel]]][as.integer(onevar$modelid[sel])]
    }
    if (is.na(ans$u[2])) {
        sel1 <- ((1:nrow(onevar))[order(onevar$u,decreasing=TRUE)])
        sel1 <- sel1[1:min(2,length(sel1))]
        u1 <- sum(onevar$u[sel1])
        if (any(nvars>1)) {
            twovar <- data.frame(block=1:length(ubest),do.call(rbind,lapply(ubest,function(z) z[2,])), stringAsFactors=FALSE)
            sel2 <- which.max(twovar$u)
            u2 <- twovar$u[sel2]
        }
        if (all(nvars==1) | u1>u2) {
            ans$u[2] <- u1
            ans$block[2] <- onevar$block[sel1][2]
            ans$modelid[2] <- paste(varidx[[onevar$block[sel1[1]]]][as.integer(onevar$modelid[sel1[1]])],varidx[[onevar$block[sel1[2]]]][as.integer(onevar$modelid[sel1[2]])],sep=',')
            ans$m[2] <- onevar$modelid[sel1[2]]
        } else {
            ans$u[2] <- u2
            ans$block[2] <- twovar$block[sel2]
            ans$modelid[2] <- paste(varidx[[twovar$block[sel2]]][as.integer(strsplit(twovar$modelid[sel2],split=',')[[1]])],collapse=',')
            ans$m[2] <- twovar$modelid[sel2]
        }
    }
    #
    #Step 4. Bound u-score for non-visited model sizes
    #
    mis <- which(is.na(ans$u[-nrow(ans)]))
    ans$u.upper <- ans$u.lower <- ans$u
    if (length(mis)>0) {
        for (i in 1:length(mis)) {
            l <- mis[i]+1
            while (is.na(ans$u[l]) & (l<=nrow(ans))) l <- l+1
            ans$u.upper[mis[i]] <- ans$u[l]
            l <- mis[i]-1
            while (is.na(ans$u[l]) & (l>=1)) l <- l-1
            ans$u.lower[mis[i]] <- ans$u[l]
        }
    }
    if (is.na(ans$u[nrow(ans)])) { ans$u.lower[nrow(ans)] <- ans$u[nrow(ans)-1] }
    return(ans)
}



## Auxiliary functions for block-diagonal design ##
###################################################

#Marginal inclusion probability & coefficient estimates for a block of variables variable conditional on a vector of phi values
margppZellnerKnownBlockDiag <- function(phi,u,umv,g,shrinkage,priorModelBlock,logscale=FALSE) {
    maxvars <- max(u$nvars)
    modelidx <- sapply(strsplit(u$modelid,split=','), function(z) { ans= rep(FALSE,maxvars); ans[as.numeric(z)] <- TRUE; return(ans) })
    priorpp <- priorModelBlock(0:maxvars,blocksize=maxvars)
    gpow <- (-(0:maxvars)/2) * log(g+1)
    margpp <- coefest <- matrix(NA,nrow=length(phi),ncol=maxvars)
    modelcoef <- shrinkage * sapply(u$modelid, function(z) { ans= rep(0,maxvars); ans[as.numeric(strsplit(z,split=',')[[1]])]= umv[[z]]$m; return(ans) })
    for (i in 1:nrow(margpp)) {
        h <- gpow[u$nvars+1] + .5*(shrinkage*u$u)/phi[i] + priorpp[u$nvars+1]
        hmax <- max(h)
        pp <- exp(h-hmax)/sum(exp(h-hmax))
        margpp[i,] <- as.vector(modelidx %*% pp)
        coefest[i,] <- as.vector(modelcoef %*% pp)
    }
    return(list(margpp=margpp,coef=coefest))
}


margppMOMKnownBlockDiag <- function(phi,u,umv,g,shrinkage,priorModelBlock,momcoef,logscale=FALSE) {
    #Note: posterior mean under MOM is approximated using MLE
    maxvars <- max(u$nvars)
    modelidx <- sapply(strsplit(u$modelid,split=','), function(z) { ans= rep(FALSE,maxvars); ans[as.numeric(z)] <- TRUE; return(ans) })
    priorpp <- priorModelBlock(0:maxvars,blocksize=maxvars)
    gpow <- (-(0:maxvars)/2) * log(g+1)
    margpp <- coefest <- matrix(NA,nrow=length(phi),ncol=maxvars)
    modelcoef <- sapply(u$modelid, function(z) { ans= rep(0,maxvars); ans[as.numeric(strsplit(z,split=',')[[1]])]= umv[[z]]$m; return(ans) })
    for (i in 1:length(phi)) {
        phipow <- 1/(phi[i]^(0:maxvars))
        logf <- log(sapply(momcoef, function(z) sum(z*phipow[1:length(z)])))
        logf <- c(logf[u$modelid[u$modelid!='']],0) #last model is null model
        h <- gpow[u$nvars+1] + .5*(shrinkage*u$u)/phi[i] + priorpp[u$nvars+1] + logf
        hmax <- max(h)
        pp <- exp(h-hmax)/sum(exp(h-hmax))
        margpp[i,] <- as.vector(modelidx %*% pp)
        coefest[i,] <- as.vector(modelcoef %*% pp)
    }
    return(list(margpp=margpp,coef=coefest))
}


#Joint model probability conditional on a single phi value
#Input
# - phi: single value for residual variance
# - sel: logical vector of length p indicating the active variables
# - block: vector of length p indicating the block that each variable belongs to
# - u: u-scores for all variables (i.e. u[sel] selects those for active variables)
# - g, shrinkage, rho, logscale: as in margppMomKnown
#Output: posterior probability of the model conditional on phi
jointppMomKnownBlockDiag <- function(phi,sel,blocks,blocksize,u,g,shrinkage=g/(1+g),priorModelBlock,momcoef,logscale=TRUE) {
    priorpp <- lapply(1:length(u), function(k) priorModelBlock(0:blocksize[k],blocksize=blocksize[k]))
    gpow <- (-(0:max(blocksize))/2) * log(g+1)
    phipow <- 1/(phi^(0:max(blocksize)))
    ans <- 0
    for (k in 1:length(u)) {
        nactive <- sum(sel[blocks==k])
        idx <- paste(which(sel[blocks==k]),collapse=',')
        logf <- log(sapply(momcoef[[k]], function(z) sum(z*phipow[1:length(z)])))
        logf <- c(logf[u[[k]]$modelid[u[[k]]$modelid!='']],0) #last model is null model
        h <- gpow[u[[k]]$nvars+1] + .5*(shrinkage*u[[k]]$u)/phi + priorpp[[k]][u[[k]]$nvars+1] + logf
        names(h) <- u[[k]]$modelid
        hmax <- max(h)
        ans <- ans + (ifelse(idx=="",h[length(h)],h[idx])-hmax) - log(sum(exp(h-hmax)))
        #exp(h[idx]) / sum(exp(h))
    }
    if (!logscale) ans <- exp(ans)
    return(ans)
}


#Evaluate quantity proportional to joint of (phi,y) on a grid
jointPhiyZellnerBlockDiag <- function(phiseq, sumy2, apos, a.phi, l.phi, shrinkage, u, g, priorModelBlock, blocksize, logscale=TRUE) {
    umax <- sapply(u, function(z) z$u[1])
    sumumax <- sum(umax)
    #Avoid overflow in cases where X'X not exactly block-diagonal
    if (sumumax > (sumy2-l.phi)) {
        umax <- umax/sumumax * (sumy2-l.phi)
        sumumax <- sum(umax)
    }
    ans <- -0.5*(sumy2-sumumax)/phiseq - (.5*apos+1)*log(phiseq)
    priorpp <- lapply(1:length(u), function(k) priorModelBlock(0:blocksize[k],blocksize=blocksize[k]))
    gpow <- (-(0:max(blocksize))/2) * log(g+1)
    for (i in 1:length(phiseq)) {
        hsum <- 0
        for (k in 1:length(u)) {
            hsum <- sum(exp(gpow[u[[k]]$nvars+1] + .5*(shrinkage*u[[k]]$u - umax[k])/phiseq[i] + priorpp[[k]][u[[k]]$nvars+1]))
            ans[i] <- ans[i] + sum(log(hsum))
        }
    }
    ans[is.infinite(ans)] <- -Inf  #avoid numerical overflow
    if (!logscale) ans <- exp(ans)
    return(ans)
}


jointPhiyMOMBlockDiag <- function(phiseq, sumy2, apos, a.phi, l.phi, shrinkage, u, umv, momcoef, g, priorModelBlock, blocksize, logscale=TRUE) {
    umax <- sapply(u, function(z) z$u[1])
    sumumax <- sum(umax)
    #Avoid overflow in cases where X'X not exactly block-diagonal
    if (sumumax > (sumy2-l.phi)) {
        umax <- umax/sumumax * (sumy2-l.phi)
        sumumax <- sum(umax)
    }
    ans <- -0.5*(sumy2-sumumax)/phiseq - (.5*apos+1)*log(phiseq)
    priorpp <- lapply(1:length(u), function(k) priorModelBlock(0:blocksize[k],blocksize=blocksize[k]))
    gpow <- (-(0:max(blocksize))/2) * log(g+1)
    for (i in 1:length(phiseq)) {
        phipow <- 1/(phiseq[i]^(0:max(blocksize)))
        for (k in 1:length(u)) {
            logf <- log(sapply(momcoef[[k]], function(z) sum(z*phipow[1:length(z)])))
            logf <- logf[u[[k]]$modelid]
            logf[length(logf)] <- 0 #last model is null model
            #logf <- sapply(umv[[k]][u[[k]]$modelid], function(z) log(eprod(shrinkage*z$m, z$V * (shrinkage*phiseq[i]), power=2)) - sum(log(diag(z$V)))) - u[[k]]$nvars * log(phiseq[i]*g) #same but slower
            hsum <- sum(exp(gpow[u[[k]]$nvars+1] + .5*(shrinkage*u[[k]]$u - umax[k])/phiseq[i] + priorpp[[k]][u[[k]]$nvars+1] + logf))
            ans[i] <- ans[i] + sum(log(hsum))
        }
    }
    ans[is.infinite(ans)] <- -Inf  #avoid numerical overflow
    if (!logscale) ans <- exp(ans)
    return(ans)
}


momPolyCoef <- function(z, shrinkage, g) {
    #Pre-compute polynomial coefficients for fast evaluation of pMOM marginal likelihood given phi
    # - z: z$m contains posterior mean, z$V the inverse of X'X, e.g. as stored in umv[[k]][[j]] and returned by uscore
    nvars <- nrow(z$m)
    ans <- double(nvars+1)
    v <- as.matrix(expand.grid(lapply(1:nvars, function(zz) 0:2)))
    h <- 1-v
    unitv <- (-2)^rowSums(v==1)
    R <- cov2cor(z$V)
    mstd <- shrinkage * z$m / sqrt(diag(z$V))
    #r=0
    ans[1]= sum(exp(sapply(1:nrow(h), function(i) { nvars*log((matrix(h[i,],nrow=1) %*% mstd)^2) }) - lfactorial(2*nvars) - nvars*log(g)) * unitv)
    if (nvars>1) {
        nvlog <- nvars*log(1+g)
        for (r in 1:(nvars-1)) {
            ss <- sapply(1:nrow(h), function(i) { r * (log(matrix(h[i,],nrow=1) %*% R %*% matrix(h[i,],ncol=1)) - log(2)) + (nvars-r)*log((matrix(h[i,],nrow=1) %*% mstd)^2) })
            ans[r+1] <- sum(exp(ss -lfactorial(r) - lfactorial(2*(nvars-r)) -nvlog + (nvars-r)*log(shrinkage)) * unitv)
        }
    }
    ans[nvars+1]= sum(exp(sapply(1:nrow(h), function(i) { nvars * (log(matrix(h[i,],nrow=1) %*% R %*% matrix(h[i,],ncol=1)) - log(2)) }) -lfactorial(r) - nvars*log(1+g)) * unitv)
    #cat('.')
    return(rev(ans))
}


enumblock <- function(l, blocksize, usedvars=0, curblocks=NULL, blockstart=1) {
    #Enumerate all models of given size l. Needless enumeration is avoided by using the within-block rankings defined by the u-statistics.
    if (blockstart==length(blocksize)) {
        ans <- c(curblocks,l)
    } else {
        minv <- max(c(0, l-sum(blocksize[(blockstart+1):length(blocksize)])))
        maxv <- min(c(blocksize[blockstart], l))
        ans <- vector("list", maxv-minv+1)
        for (i in minv:maxv) {
            ans[[i+1]] <- enumblock(l-i, blocksize=blocksize, usedvars=usedvars+i, curblocks=c(curblocks,i), blockstart=blockstart+1)
        }
        ans <- do.call(rbind,ans)
    }
    return(ans)
}


uscore <- function(y,XtX,Xty) {
    #Compute u-score, least squares estimate and its covariance for a given model
    V <- solve(XtX)
    m <- V %*% Xty
    u <- t(Xty) %*% m
    return(list(u=u,m=m,V=V))
}




## Auxiliary functions for orthogonal design ##
###############################################

#Marginal inclusion probability for a single variable conditional on a vector of phi values
# Input
# - phi: vector of phi values (residual variance)
# - u: u statistic for a single variable
# - g: prior dispersion parameter
# - shrinkage: g/(1+g)
# - rho: prior inclusion probability
# - logscale: set to TRUE to obtain log inclusion probability
# Output: vector of marginal inclusion probabilities conditional on the given phi's
margppZellnerKnownOrtho <- function(phi,u,g,shrinkage=g/(1+g),rho,logscale=TRUE) {
  ans <- -log(1 + sqrt(1+g) * exp(-.5 * shrinkage*u/phi) * (1-rho)/rho)
  if (!logscale) ans <- exp(ans)
  return(ans)
}

margppMomKnownOrtho <- function(phi,u,g,shrinkage=g/(1+g),rho,logscale=TRUE) {
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
jointppZellnerKnownOrtho <- function(phi,sel,u,g,shrinkage=g/(1+g),rho,logscale=TRUE) {
    uu <- shrinkage*u/phi
    odds01 <- sqrt(1+g) * exp(-.5 * uu) * (1-rho)/rho
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

#Joint model probability conditional on a single phi value
#Input
# - phi: single value for residual variance
# - sel: model indicator, i.e. vector indicating indexes of active variables under current model
# - u: u-scores for all variables (i.e. u[sel] selects those for active variables)
# - g, shrinkage, rho, logscale: as in margppMomKnown
#Output: posterior probability of the model conditional on phi
jointppMomKnownOrtho <- function(phi,sel,u,g,shrinkage=g/(1+g),rho,logscale=TRUE) {
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
jointPhiyZellnerOrtho <- function(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE) {
    #Code version 1
    fseqnew <- do.call(cbind,lapply(phiseq, function(phi) exp(0.5 * shrinkage * u / phi)))
    ans <- -0.5*sumy2/phiseq - (.5*apos+1)*log(phiseq) + colSums(log(1 + rho*(fseqnew/sqrt(1+g) -1)))
    #Code version 2. Smarter algebra but often causes overflow
    #ans <- double(length(phiseq))
    #ct <- -0.5*sumy2/length(u)
    #for (i in 1:length(phiseq)) {
    #    fseqnew <- exp((ct + 0.5 * shrinkage * u) / phiseq[i])
    #    ans[i] <- -(.5*apos+1)*log(phiseq[i]) + sum(log((1-rho)*exp(ct/phiseq[i]) + rho*fseqnew/sqrt(1+g)))
    #}
    ans[is.infinite(ans)] <- -Inf  #avoid numerical overflow
    if (!logscale) ans <- exp(ans)
    return(ans)
}

jointPhiyMOMOrtho <- function(phiseq, sumy2, apos, shrinkage, u, g, rho, logscale=TRUE) {
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
