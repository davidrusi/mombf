###
### ggm.R
###

### Methods for msfit_ggm objects

plot.msfit_ggm= function(x, y, ...) {
  postSample= x$postSample[,x$indexes[1,] != x$indexes[2,],drop=FALSE]
  margppcum= apply(postSample !=0, 2, cumsum) / (1:nrow(postSample))
  plot(margppcum[,1], type='l', ylim=c(0,1), xlab='Iteration', ylab='Marginal posterior inclusion probabilities')
  if (ncol(margppcum)>1) for (i in 2:min(5000,ncol(margppcum))) lines(margppcum[,i])
}


setMethod("show", signature(object='msfit_ggm'), function(object) {
  cat('Gaussian graphical model (msfit_ggm object) with ',object$p,'variables\n')
  cat("Use coef() to get BMA estimates, posterior intervals and posterior marginal prob of entries being non-zero\n")
}
)

coef.msfit_ggm <- function(object,...) {
  if (object$almost_parallel) stop("coef not yet implemented for almost_parallel")
  m= Matrix::colMeans(object$postSample)
  ci= sparseMatrixStats::colQuantiles(object$postSample, prob=c(0.025,0.975))
  if (object$samplerPars['sampler'] == 'Gibbs') {
    ans= cbind(t(object$indexes), m, ci, object$margpp)
  } else {
    ans= cbind(t(object$indexes), m, ci, Matrix::colMeans(object$postSample != 0))
  }
  colnames(ans)[-1:-2]= c('estimate','2.5%','97.5%','margpp')
  return(ans)
}

icov <- function(fit, threshold) {
  if (fit$almost_parallel) stop("coef not yet implemented for almost_parallel")
  if (!inherits(fit, 'msfit_ggm')) stop("Argument fit must be of class msfit_ggm")
  m= Matrix::colMeans(fit$postSample)
  if (!missing(threshold)) {
      margpp= Matrix::colMeans(fit$postSample != 0)
      sel= (margpp >= threshold)
      ans= Matrix::sparseMatrix(i= fit$indexes[1,sel], j= fit$indexes[2,sel], x=m[sel], dims=c(fit$p,fit$p), symmetric=TRUE)
  } else {
      ans= matrix(nrow=fit$p,ncol=fit$p)
      ans[upper.tri(ans,diag=TRUE)]= m
      ans[lower.tri(ans)]= t(ans)[lower.tri(ans)]
  }
  return(ans)
}


### Model selection routines

modelSelectionGGM= function(y, priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), center=TRUE, scale=TRUE, almost_parallel= FALSE, sampler='Gibbs', niter=10^3, burnin= round(niter/10), pbirth=0.5, nbirth, Omegaini='glasso-ebic', verbose=TRUE) {
  #Check input args
  if (!is.matrix(y)) y = as.matrix(y)
  if (ncol(y) <=1) stop("y must have at least 2 columns")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!(sampler %in% c('Gibbs','birthdeath','zigzag'))) stop("sampler must be 'Gibbs', 'birthdeath' or 'zigzag'")
  y = scale(y, center=center, scale=scale)
    
  #Format prior parameters
  prCoef= formatmsPriorsMarg(priorCoef=priorCoef, priorVar=priorDiag)
  prCoef= as.list(c(priorlabel=prCoef$priorCoef@priorDistr, prCoef[c('prior','tau','lambda')]))
  prModel= as.list(c(priorlabel=priorModel@priorDistr, priorPars= priorModel@priorPars))
    
  #Format posterior sampler parameters
  samplerPars= format_GGM_samplerPars(sampler, p=ncol(y), niter, burnin, pbirth, nbirth, verbose)
    
  #Initial value for sampler
  Omegaini= initialEstimateGGM(y, Omegaini)
    
  #Call C++ function
  if (!almost_parallel) {
    ans= modelSelectionGGMC(y, prCoef, prModel, samplerPars, Omegaini)
    postSample= Matrix::t(ans$postSample)
  } else {
    ans= GGM_Gibbs_parallelC(y, prCoef, prModel, samplerPars, Omegaini)
    postSample= lapply(ans, Matrix::t)
  }
    

  #Return output
  priors= list(priorCoef=priorCoef, priorModel=priorModel, priorDiag=priorDiag)

  A= diag(ncol(y))
  indexes= rbind(row(A)[upper.tri(row(A),diag=TRUE)], col(A)[upper.tri(row(A),diag=TRUE)])
  rownames(indexes)= c('row','column')

  ans= list(postSample=postSample, margpp=ans$margpp, priors=priors, p=ncol(y), indexes=indexes, samplerPars=samplerPars, almost_parallel=almost_parallel)

  new("msfit_ggm",ans)
}



#Format posterior sampling parameters to pass onto C++
format_GGM_samplerPars= function(sampler, p, niter, burnin, pbirth, nbirth, verbose) {
  if (missing(nbirth)) {
      nbirth= as.integer(min(p, max(log(p), 10)))
  } else {
      nbirth= as.integer(nbirth)
  }
  samplerPars= list(sampler, as.integer(niter), as.integer(burnin), pbirth, nbirth, as.integer(ifelse(verbose,1,0)))
  names(samplerPars)= c('sampler','niter','burnin','pbirth','nbirth','verbose')
  return(samplerPars)
}

#Multivariate normal density given the precision matrix
dmvnorm_prec <- function(x, sigmainv, logdet.sigmainv, mu = rep(0, ncol(sigmainv)), log = FALSE) {
  if (missing(logdet.sigmainv)) logdet.sigmainv= determinant(sigmainv,logarithm=TRUE)$modulus
  d= mahalanobis(x, center=mu, cov=sigmainv, inverted=TRUE)
  ans= -0.5 * (d + ncol(sigmainv) * log(2*pi) - logdet.sigmainv)
  if(!log) ans= exp(ans)
  return(ans)
}

#Initial estimate of precision matrix
initialEstimateGGM= function(y, Omegaini) {
    
  if (is.character(Omegaini)) {
    ans= initGGM(y, Omegaini)
  } else if (inherits(Omegaini, "matrix")) {
    ans= Matrix(Omegaini, sparse=TRUE)
  } else if (inherits(Omegaini, "dgCMatrix", "ddiMatrix")) {
    ans= Omegaini
  } else stop("Invalid Omegaini. It must be of class matrix, dgCMatrix or ddiMatrix")
  return(ans)
    
}

initGGM= function(y, Omegaini) {
    
  if (Omegaini=='null') {
        
    ans= sparseMatrix(1:ncol(y), 1:ncol(y), x=rep(1,ncol(y)), dims=c(ncol(y),ncol(y)))
        
  } else if (Omegaini %in% c('glasso-bic','glasso-ebic')) {

    method= ifelse(Omegaini == 'glasso-bic', 'BIC', 'EBIC')

    sfit= glassopath(cov(y), trace=0)  #from package glasso
    bic= glasso_getBIC(y=y, sfit=sfit, method=method)

    topbic= which.min(bic)
    found= (topbic != 1) & (topbic != length(bic))

    niter= 1
    while ((!found) & (niter<=5)) {
        
      if (topbic==1) {
          rholist= seq(sfit$rholist[1]/10, sfit$rholist[1], length=10)
      } else {
          l= sfit$rholist[length(sfit$rholist)]
          rholist= seq(l, 10*l, length=10)
      }

      sfit= glassopath(cov(y), rholist=rholist, trace=0)
      bic= glasso_getBIC(y=y, sfit=sfit, method=method)
      
      topbic= which.min(bic)
      found= (topbic != 1) & (topbic != length(bic))
      niter= niter+1

    }
        
    ans= Matrix(sfit$wi[,,which.min(bic)], sparse=TRUE)
      
  } else {
    stop("Omegaini must be 'null', 'glasso-ebic' or 'glasso-bic'")
  }
    
  return(ans)

}


glasso_getBIC= function(y, sfit, method='EBIC') {
  logl= npar= double(length(sfit$rholist))
  for (i in 1:length(sfit$rholist)) {
    logl[i]= sum(dmvnorm_prec(y, sigmainv=sfit$wi[,,i], log=TRUE))
    npar[i]= sum(sfit$wi[,,i] != 0)
  }
  if (method == 'BIC') {
    ans= -2*logl + npar * log(nrow(y))
  } else if (method == 'EBIC') {
    ans= -2*logl + npar * (log(nrow(y)) + log(ncol(y)^2))
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
    idx= which(model)
    idx= idx[idx != colsel]
    s= S[idx, colsel, drop=FALSE]
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
