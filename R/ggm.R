###
### ggm.R
###

### Methods for msfit_ggm objects

plot.msfit_ggm= function(x, y, ...) {
  postSample= x$postSample[,x$indexes[1,] != x$indexes[2,],drop=FALSE]
  margppcum= apply(postSample !=0, 2, cumsum) / seq_len(nrow(postSample))
  plot(margppcum[,1], type='l', ylim=c(0,1), xlab='Iteration', ylab='Marginal posterior inclusion probabilities')
  if (ncol(margppcum)>1) for (i in 2:min(5000,ncol(margppcum))) lines(margppcum[,i])
}


setMethod("show", signature(object='msfit_ggm'), function(object) {
  cat('Gaussian graphical model (msfit_ggm object) with ',object$p,'variables\n')
  cat("Use coef() to get BMA estimates, posterior intervals and posterior marginal prob of entries being non-zero\n")
  cat("use postProb() for posterior model probabilities")
}
)

coef.msfit_ggm <- function(object,...) {
  m= Matrix::colMeans(object$postSample)
  ci= sparseMatrixStats::colQuantiles(object$postSample, prob=c(0.025,0.975))
  #if (object$almost_parallel=='none') {
  ans= cbind(t(object$indexes), m, ci, object$margpp) #use Rao-Blackwellized edge inclusion probabilities
  #} else {
  #ans= cbind(t(object$indexes), m, ci, Matrix::colMeans(object$postSample != 0)) #MCMC edge inclusion proportions
  #}
  colnames(ans)[-1:-2]= c('estimate','2.5%','97.5%','margpp')
  return(ans)
}

icov <- function(fit, threshold) {
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


proportion_visited_models= function(postSample) {
  modelpp <- apply(postSample != 0, 1, function(z) paste(which(z),collapse=','))
  modelpp <- table(modelpp)/length(modelpp)
  modelpp <- data.frame(modelid=names(modelpp), pp=as.numeric(modelpp))
  return(modelpp)
}

setMethod("postProb", signature(object='msfit_ggm'), function(object, nmax, method='norm') {
if (!is.null(object$models)) {
  ans= object$models
} else {
  param_ids= paste('(',object$indexes[1,], ',', object$indexes[2,], ')',sep='')
  npar= ncol(object$indexes)
  #Compute posterior probabilities
  modelpp = proportion_visited_models(object$postSample)
  #Compute model identifiers
  tmp= strsplit(modelpp$modelid, ',')
  modelid= t(sapply(tmp, function(z) { ans= rep(FALSE,npar); ans[as.integer(z)]= TRUE; return(ans) }))
  colnames(modelid)= param_ids
  modelid= Matrix::Matrix(modelid, sparse=TRUE)
  #Sort models in decreasing probability
  o= order(modelpp$pp,decreasing=TRUE)
  modelpp= modelpp[o,]
  modelid= modelid[o,]
  #Return nmax models    
  if (!missing(nmax)) {
      modelpp <- modelpp[1:nmax,]
      modelid = modelid[1:nmax,]
  }
  #Return output
  ans= list(modelid= modelid, pp=modelpp$pp)
}
return(ans)
}
)



### Model selection routines

modelSelectionGGM= function(y, priorCoef=normalidprior(tau=1), priorModel=modelbinomprior(1/ncol(y)), priorDiag=exponentialprior(lambda=1), center=TRUE, scale=TRUE, almost_parallel= "regression", prob_parallel=0.5, tempering=0.5, truncratio= 100, save_proposal=FALSE, niter=10^3, burnin= round(niter/10), updates_per_iter= ceiling(sqrt(ncol(y))), updates_per_column= 10, sampler='birthdeath', pbirth=0.75, pdeath=0.5*(1-pbirth), bounds_LIT, Omegaini='glasso-ebic', verbose=TRUE) {
  #Check input args
  if (!is.matrix(y)) y = as.matrix(y)
  p= ncol(y);
  if (p <=1) stop("y must have at least 2 columns")
  if (!is.numeric(y)) stop("y must be numeric")
  if (!(sampler %in% c('Gibbs','birthdeath','LIT'))) stop("sampler must be 'Gibbs', 'birthdeath', or 'LIT'")
  if (tempering < 0) stop("tempering cannot be negative")
  y = scale(y, center=center, scale=scale)
    
  #Format prior parameters
  prCoef= formatmsPriorsMarg(priorCoef=priorCoef, priorVar=priorDiag)
  prCoef= as.list(c(priorlabel=prCoef$priorCoef@priorDistr, prCoef[c('prior','tau','lambda')]))
  prModel= as.list(c(priorlabel=priorModel@priorDistr, priorPars= priorModel@priorPars))
    
  #Format posterior sampler parameters
  if (missing(bounds_LIT)) bounds_LIT= log(c("lbound_death"=1/p, "ubound_death"= 1, "lbound_birth"=1/p, "ubound_birth"=p))
  samplerPars= format_GGM_samplerPars(sampler, p=p, niter=niter, burnin=burnin, updates_per_iter=updates_per_iter, updates_per_column = updates_per_column, pbirth=pbirth, pdeath=pdeath, prob_parallel=prob_parallel, tempering=tempering, truncratio=truncratio, almost_parallel=almost_parallel, bounds_LIT=bounds_LIT, verbose=verbose)
    
  #Initial value for sampler
  if (verbose) cat(" Obtaining initial parameter estimate...")
  Omegaini= initialEstimateGGM(y, Omegaini)
  if (verbose) cat(" Done\n")
    
  #Call C++ function
  proposal= proposaldensity= NULL
  if (almost_parallel == 'none') {
    ans= modelSelectionGGMC(y, prCoef, prModel, samplerPars, Omegaini)
    postSample= Matrix::t(ans$postSample)
    margpp= ans$margpp
    prop_accept= ans$prop_accept
  } else {
    ans= modelSelectionGGM_parallelC(y, prCoef, prModel, samplerPars, Omegaini)
    postSample= Matrix::t(ans[[1]])
    margpp= ans[[2]]
    prop_accept= ans[[3]]
    if (save_proposal) {
      proposal= lapply(ans[4:(4+p-1)], Matrix::t)
      for (i in seq_len(length(proposal))) proposal[[i]]@x = as.double(proposal[[i]]@x)
      proposaldensity= Matrix::t(ans[[length(ans)]])
    }
  }
    

  #Return output
  priors= list(priorCoef=priorCoef, priorModel=priorModel, priorDiag=priorDiag)

  A= diag(p)
  indexes= rbind(row(A)[upper.tri(row(A),diag=TRUE)], col(A)[upper.tri(row(A),diag=TRUE)])
  rownames(indexes)= c('row','column')

  ans= list(postSample=postSample, prop_accept=prop_accept, proposal=proposal, proposaldensity=proposaldensity, margpp=margpp, priors=priors, p=p, indexes=indexes, samplerPars=samplerPars, almost_parallel=almost_parallel)

  new("msfit_ggm",ans)
}



#Format posterior sampling parameters to pass onto C++
format_GGM_samplerPars= function(sampler, p, niter, burnin, updates_per_iter, updates_per_column, pbirth, pdeath, prob_parallel, tempering, truncratio, almost_parallel, bounds_LIT, verbose) {
  if (missing(updates_per_column)) {
      updates_per_column= as.integer(min(p, max(log(p), 10)))
  } else {
      updates_per_column= as.integer(updates_per_column)
  }
  if (!(almost_parallel %in% c('none','regression','in-sample'))) stop("almost_parallel must be 'none', 'regression' or 'in-sample'")
  if ((updates_per_iter < 1) || (updates_per_iter > p)) stop("updates_per_iter must be between 1 and ncol(y)")
  if ((pbirth + pdeath) > 1) stop("pbirth + pdeath must be <=1")
  if (sampler == "LIT") {
    psum= (pbirth + pdeath)
    pbirth = pbirth / psum
    pdeath = pdeath / psum
  }

  if (!is.vector(bounds_LIT)) stop("bounds_LIT must be a vector")
  if (!all(c("lbound_death", "ubound_death", "lbound_birth", "ubound_birth") %in% names(bounds_LIT))) stop("bounds_LIT must be a named vector containing 'lbound_death', 'ubound_death', 'lbound_birth', 'ubound_birth'")
  
  samplerPars= list(sampler, as.integer(niter), as.integer(burnin), as.integer(updates_per_iter), as.integer(updates_per_column), as.double(pbirth), as.double(pdeath), as.double(prob_parallel), as.double(tempering), as.double(truncratio), almost_parallel, as.double(bounds_LIT["lbound_death"]), as.double(bounds_LIT["ubound_death"]), as.double(bounds_LIT["lbound_birth"]), as.double(bounds_LIT["ubound_birth"]), as.integer(ifelse(verbose,1,0)))
  names(samplerPars)= c('sampler','niter','burnin','updates_per_iter','updates_per_column','pbirth','pdeath','prob_parallel','tempering','truncratio','almost_parallel','lbound_death','ubound_death','lbound_birth','ubound_birth','verbose')
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
  } else {
    if (ncol(Omegaini) != nrow(Omegaini)) stop("Omegaini must be a square matrix")
    if (ncol(y) != ncol(Omegaini)) stop("ncol(Omegaini) must be equal to ncol(y)")
    if (inherits(Omegaini, "matrix")) {
    ans= Matrix::Matrix(Omegaini, sparse=TRUE)
  } else if (inherits(Omegaini, c("dgCMatrix", "ddiMatrix", "dsCMatrix"))) {
    ans= Omegaini
  } else stop("Invalid Omegaini. It must be of class matrix, dgCMatrix, dsCMatrix or ddiMatrix")
  }
  return(ans)
    
}



initGGM= function(y, Omegaini) {
  
  maxiter= ifelse (ncol(y) <=100, 5, 1)  #use 1 iteration for large p, to avoid excessive init time
  
  if (Omegaini=='null') {
    
    ans= Matrix::sparseMatrix(seq_len(ncol(y)), seq_len(ncol(y)), x=rep(1,ncol(y)), dims=c(ncol(y),ncol(y)))
    
  } else if (Omegaini %in% c('glasso-bic','glasso-ebic')) {
    if (Omegaini == 'glasso-bic') {
      method <- "BIC"
      #gamma = 0
    } else if (Omegaini == 'glasso-ebic') {
      #gamma = 0.5
      method <- "EBIC"
    }
    
    sfit= huge::huge(y, method="glasso", scr=TRUE, verbose=FALSE, nlambda = 10)
    #sfit2= huge::huge.select(sfit, criterion="ebic", ebic.gamma = gamma, verbose=FALSE)
    ebic= glasso_getBIC(y=y, sfit=sfit, method=method)
    
    #topbic= which.min(sfit2$ebic.score)
    topbic= which.min(ebic)
    found= (topbic != 1) & (topbic != 10)
    
    niter= 1
    while ((!found) && (niter<=maxiter)) {
      ## huge wants a decreasing sequence
      if((niter == 1) && (topbic == 1)){## No need for lambda bigger
        break
      }
      if (topbic==10) {
        rholist= seq(sfit$lambda[10], sfit$lambda[10]/10, length=10)
      } else {
        rholist= seq(10*sfit$lambda[1], sfit$lambda[1], length=10)
      }
      
      sfit= huge::huge(y, lambda=rholist, method="glasso", scr=TRUE, verbose=FALSE)
      #sfit2= huge::huge.select(sfit, criterion="ebic", ebic.gamma = gamma, verbose=FALSE)
      ebic= glasso_getBIC(y=y, sfit=sfit, method=method)
      
      #topbic= which.min(sfit2$ebic.score)
      topbic= which.min(ebic)
      found= (topbic != 1) & (topbic != 10)
      niter= niter+1
      
    }
    
    ans= Matrix::Matrix(sfit$icov[[topbic]], sparse=TRUE)
    
  } else {
    stop("Omegaini must be 'null', 'glasso-ebic' or 'glasso-bic'")
  }
  
  return(ans)
  
}



glasso_getBIC = function(y, sfit, method='EBIC') {
  #Extraction of BIC/EBIC, copied from huge.select in package huge
  n= nrow(y); d= ncol(y)
  if (method == 'BIC') {
    gamma = 0
  } else if (method == 'EBIC') {
    gamma = 0.5
  }
  #ans = -n * sfit$loglik + log(n) * sfit$df + 4 * gamma * log(d) * sfit$df
  #Manual implementation
  logl= npar= double(length(sfit$icov))
  for (i in seq_len(length(logl))) {
    logl[i]= sum(dmvnorm_prec(y, sigmainv=sfit$icov[[i]], log=TRUE))
    npar[i]= (sum(sfit$icov[[i]] != 0) - d)/2
  }
  if (method == 'BIC') {
    ans= -2*logl + npar * log(n)
  } else if (method == 'EBIC') {
    ans= -2*logl + npar * (log(n) + 4 * gamma * log(d))
  }
  return(ans)
}


#For objects from package huge (deprecated)
#glasso_getBIC = function(y, sfit, method='EBIC') {
#  #Extraction of BIC/EBIC, copied from huge.select in package huge
#  n= nrow(y); d= ncol(y)
#  if (method == 'BIC') {
#    gamma = 0
#  } else if (method == 'EBIC') {
#    gamma = 0.5
#  }
#  #ans = -n * sfit$loglik + log(n) * sfit$df + 4 * gamma * log(d) * sfit$df
#  #Manual implementation (deprecated)    
#  logl= npar= double(length(sfit$icov))
#  for (i in 1:length(logl)) {
#    logl[i]= sum(dmvnorm_prec(y, sigmainv=sfit$icov[[i]], log=TRUE))
#    npar[i]= (sum(sfit$icov[[i]] != 0) - d)/2
#  }
#  if (method == 'BIC') {
#    ans= -2*logl + npar * log(n)
#  } else if (method == 'EBIC') {
#    ans= -2*logl + npar * (log(n) + 4 * gamma * log(d))
#  }
#  return(ans)
#}




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
