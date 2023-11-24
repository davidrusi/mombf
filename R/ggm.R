modelSelectionGGM= function(y, priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), sampler='Gibbs', niter=10^3, burnin= round(niter/10)) {
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
  #Call C++ function
  ans= modelSelectionGGMC(y, prCoef, prModel, samplerPars)

}
