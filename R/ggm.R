modelSelectionGGM= function(y, priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), sampler='Gibbs', niter=10^3, burnin= round(niter/10), initialize='null') {
  #Check input args
  if (!is.matrix(y)) y = as.matrix(y)
  if (ncol(y) <=1) stop("y must have at least 2 columns")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!(sampler %in% c('Gibbs','zigzag'))) stop("sampler must be 'Gibbs' or 'zigzag'")
  #Format prior parameters
  prCoef= formatmsPriorsMarg(priorCoef=priorCoef, priorVar=priorDiag)
  prCoef= as.list(c(priorlabel=prCoef$priorCoef@priorDistr, prCoef[c('prior','tau','lambda')]))
  prModel= as.list(c(priorlabel=priorModel@priorDistr, priorPars= priorModel@priorPars))
  #Format posterior sampler parameters
  samplerPars= list(sampler, as.integer(niter), as.integer(burnin))
  names(samplerPars)= c('sampler','niter','burnin')
  #Initial value for sampler
  Omegaini= initialEstimateGGM(y, initialize)
    
  #Call C++ function
  ans= modelSelectionGGMC(y, prCoef, prModel, samplerPars, Omegaini)

}

#Multivariate normal density given the precision matrix
dmvnorm_prec <- function(x, sigmainv, logdet.sigmainv, mu = rep(0, ncol(sigmainv)), log = FALSE) {
  if (missing(logdet.sigmainv)) logdet.sigmainv= determinant(sigmainv,log=TRUE)$modulus
  d= mahalanobis(x, center=mu, cov=sigmainv, inverted=TRUE)
  ans= -0.5 * (d + ncol(sigmainv) * log(2*pi) - logdet.sigmainv)
  if(!log) ans= exp(ans)
  return(ans)
}

initialEstimateGGM= function(y, initialize) {

  if (initialize=='null') {
        
    ans= Matrix(nrow=ncol(y), ncol=ncol(y), sparse=TRUE)
        
  } else if (initialize=='glasso') {
        
    sfit= glassopath(cov(y), trace=0)  #from package glasso
    #sfit= glassopath(cov(y),rho=seq(.005,.05,length=20),trace=0)
     
    logl= double(length(sfit$rholist))
    for (i in 1:length(sfit$rholist)) {
        logl[i]= sum(dmvnorm_prec(y, sigmainv=sfit$w[,,i], log=TRUE))
        npar= sum(sfit$w[,,i] != 0)
    }
    bic= -2*logl + npar * log(nrow(y))
    ans= Matrix(sfit$w[,,which.min(bic)], sparse=TRUE)
  } else {
    stop("initialize must be 'null' or 'glasso'")
  }
    
  return(ans)

}
