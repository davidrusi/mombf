modelSelectionGGM= function(y, priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), sampler='Gibbs', niter=10^3, burnin= round(niter/10), initialize='null') {
  #Check input args
  if (!is.matrix(y)) y = as.matrix(y)
  if (ncol(y) <=1) stop("y must have at least 2 columns")
  if (!is.numeric(y)) stop("y must be numeric")
  #Format prior parameters
  prCoef= formatmsPriorsMarg(priorCoef=priorCoef, priorVar=priorDiag)
  prCoef= as.list(c(priorlabel=prCoef$priorCoef@priorDistr, prCoef[c('prior','tau','lambda')]))
  prModel= as.list(c(priorlabel=priorModel@priorDistr, priorPars= priorModel@priorPars))
  #Format posterior sampler parameters
  samplerPars= list(sampler, as.integer(niter), as.integer(burnin))
  names(samplerPars)= c('sampler','niter','burnin')
  #Initial value for sampler
  initialEstimate(y, initialize)
    
  #Call C++ function
  ans= modelSelectionGGMC(y, prCoef, prModel, samplerPars)

}


initialEstimateGGM= function(y, initialize) {

  if (initialize=='null') {
        
    ans= Matrix(nrow=ncol(y), ncol=ncol(y), sparse=TRUE)
        
  } else if (initialize=='glasso') {
        
    sfit= glassopath(cov(y), trace=0)  #from package glasso
    #sfit= glassopath(cov(y),rho=seq(.005,.05,length=20),trace=0)
     
    logl= double(length(sfit$rholist))
    for (i in 1:length(sfit$rholist)) {
        logl[i]= sum(dmvnorm(y, sigma=sfit$w[,,i], log=TRUE))
        npar= sum(sfit$w[,,i] != 0)
    }
    bic= -2*logl + npar * log(nrow(y))
    ans= sfit$wi[,,which.min(bic)]
  } else {
    stop("initialize must be 'null' or 'glasso'")
  }
    
  return(ans)

}
