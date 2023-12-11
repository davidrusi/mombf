modelSelectionGGM= function(y, priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), center=TRUE, scale=TRUE, sampler='Gibbs', niter=10^3, burnin= round(niter/10), initialize='null') {
  #Check input args
  if (!is.matrix(y)) y = as.matrix(y)
  if (ncol(y) <=1) stop("y must have at least 2 columns")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!(sampler %in% c('Gibbs','zigzag'))) stop("sampler must be 'Gibbs' or 'zigzag'")
  y = scale(y, center=center, scale=scale)
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


#Model log-joint (log marginal likelihood + log model prior) for the non-zero entries in colsel of Omega
#This is an R version of the C++ function GGM_rowmarg, build for debugging purposes
GGM_rowmargR= function(y, colsel, Omega, model= (Omega[,colsel]!=0), priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), center=TRUE, scale=TRUE) {
  #Format input
  y = scale(y, center=center, scale=scale)
  prCoef= formatmsPriorsMarg(priorCoef=priorCoef, priorVar=priorDiag)
  prCoef= as.list(c(priorlabel=prCoef$priorCoef@priorDistr, prCoef[c('prior','tau','lambda')]))
  prModel= as.list(c(priorlabel=priorModel@priorDistr, priorPars= priorModel@priorPars))
  lambda= prCoef[["lambda"]]
  tau= prCoef[["tau"]]
  if (sum(model)>1) {
    S= t(y) %*% y
    #Obtain Omegainv, select submatrix for model
    Omegainv= solve(Omega[-colsel,-colsel,drop=FALSE])
    Omegainv_model= Omegainv[model[-colsel], model[-colsel], drop=FALSE]
    #Obtain posterior mean m and covariance Uinv
    U= (S[colsel,colsel] + lambda) * Omegainv_model + diag(1/tau, nrow=nrow(Omegainv_model))
    Uinv= solve(U)
    s= S[which(model)[-colsel], colsel, drop=FALSE]
    m= Uinv %*% s
    #Log marginal likelihood
    logmarg= 0.5 * t(m) %*% U %*% m - 0.5 * sum(model) * log(tau) - 0.5 * log(det(U))
  } else {
    logmarg= - 0.5 * sum(model) * log(tau)
    Omegainv= Omegainv_model= Uinv= m= NULL
  }
  #Log model prior
  p= 1/ncol(y); nsel= sum(model[-colsel])
  logprior= nsel * log(p) + (length(model) - 1 - nsel) * log(1-p)
  #Returns
  ans= list(logjoint=logmarg+logprior, Omegainv=Omegainv, Omegainv_model=Omegainv_model, m=m, Uinv=Uinv)
  return(ans)
}
