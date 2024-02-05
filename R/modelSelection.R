###
### modelSelection.R
###

### Methods for msfit objects


plot.msfit= function(x, y, ...) {
  if (!x$enumerate) {
    margppcum= apply(x$postSample, 2, cumsum) / (1:nrow(x$postSample))
    plot(margppcum[,1], type='l', ylim=c(0,1), xlab='Iteration', ylab='Marginal posterior inclusion probabilities')
    if (ncol(margppcum)>1) for (i in 2:ncol(margppcum)) lines(margppcum[,i])
  } else {
    stop("plot.msfit produces an MCMC convergence diagnostic plot, but MCMC output is not available (you probably set enumerate=FALSE)")
  }
}


setMethod("show", signature(object='msfit'), function(object) {
  cat('msfit object with outcome of type',object$outcometype,',',object$p,'covariates and',object$family,'error distribution\n')
  ifelse(any(object$postMode!=0), paste('  Posterior mode: covariate',which(object$postMode==1)), '  Posterior mode: null model')
  cat("Use postProb() to get posterior model probabilities\n")
  cat("Use coef() or predict() to get BMA estimates and intervals for parameters or given covariate values\n")
  cat("Elements $margpp, $postMode, $postSample and $coef contain further information (see help('msfit') and help('modelSelection') for details)\n")
}
)


hasPostSampling <- function(object) {
#Sends an error message if posterior sampling is not implemented for the priors and outcome type of the msFit object
#
  #List combinations for which posterior sampling is implemented
  hassamples= data.frame(matrix(NA,nrow=10,ncol=4))
  names(hassamples)= c('outcometype','family','priorCoef','priorGroup')
  hassamples[1,] =    c('Continuous','normal',   'pMOM',   'pMOM')
  hassamples[2,] =    c('Continuous','normal',  'peMOM',  'peMOM')
  hassamples[3,] =    c('Continuous','normal',  'piMOM',  'piMOM')
  hassamples[4,] =    c('Continuous','normal','zellner','zellner')
  hassamples[5,] =    c('Continuous','normal','normalid','normalid')
  hassamples[6,] =    c('glm','binomial','bic','bic')
  hassamples[7,] =    c('glm','binomial logit','bic','bic')
  hassamples[8,] =    c('glm','binomial probit','bic','bic')
  hassamples[9,] =    c('glm','gamma inverse','bic','bic')
  hassamples[10,]=    c('glm','inverse.gaussian 1/mu^2','bic','bic')
  hassamples[11,]=    c('glm','poisson','bic','bic')
  hassamples[12,]=    c('glm','poisson log','bic','bic')
  hassamples[13,]=    c('Survival','Cox','bic','bic')
  #hassamples[14,]=    c('Survival','normal','bic','bic') #to be added
  #Check if there's variable groups
  hasgroups= (length(object$groups) > length(unique(object$groups)))
  outcometype= object$outcometype; family= object$family; priorCoef= object$prior$priorCoef@priorDistr; priorGroup= object$prior$priorGroup@priorDistr
  outcomefam= paste(outcometype,family,sep=',')
  if (hasgroups) {
      outcomefamprior= paste(outcomefam,priorCoef,priorGroup,sep=',')
      avail_outcomefamprior= apply(hassamples,1,paste,collapse=',')
  } else {
      outcomefamprior= paste(outcomefam,priorCoef,sep=',')
      avail_outcomefamprior= apply(hassamples[,1:3],1,paste,collapse=',')
  }
  found= outcomefam %in% apply(hassamples[,1:2],1,paste,collapse=',')
  exactsampling= outcomefamprior  %in% avail_outcomefamprior
  if (!found) {
    cat("Inference on parameters currently only available for the following settings: \n\n")
    print(hassamples)
    cat("\n")
    stop("Inference on parameters not implemented for outcometype= ",outcometype,", family=",family,", priorCoef=",priorCoef,", priorGroup=",priorGroup)
  } else {
    if (!exactsampling) warning("Exact posterior sampling not implemented, using Normal approx instead")
  }
}




coef.msfit <- function(object,...) {
    hasPostSampling(object)
    th= rnlp(msfit=object,niter=10^4)
    ct= (object$stdconstants[-1,'scale']==0)
    if (!is.null(names(object$margpp))) {
        nn= names(object$margpp)
    } else if (!is.null(colnames(th))) {
        nn= colnames(th)
    } else { nn= paste('beta',1:ncol(th))  }
    if (any(ct)) {
        if (ncol(th) > length(object$margpp)) {
            if (object$family %in% c('binomial','poisson')) {
              margpp= object$margpp
            } else {
              margpp= c(object$margpp,1)
              nn= c(nn,'phi')
            }
        } else { margpp= object$margpp }
    } else {
        if (ncol(th) > length(object$margpp)) {
            if (object$family %in% c('binomial','poisson')) {
              margpp= c(mean(th[,1]!=0),object$margpp)
              nn= c('intercept',nn)
            } else {
              margpp= c(mean(th[,1]!=0),object$margpp,1)
              nn= c('intercept',nn,'phi')
            }
        } else {
            margpp= c(mean(th[,1]!=0),object$margpp)
            nn= c('intercept',nn)
        }
    }
    ans= cbind(colMeans(th),t(apply(th,2,quantile,probs=c(.025,0.975))),margpp=margpp)
    colnames(ans)= c('estimate','2.5%','97.5%','margpp')
    rownames(ans)= nn
    return(ans)
}


#Return point estimate, 95% interval and posterior model prob for nmax top models
setMethod("coefByModel", signature(object='msfit'), function(object, maxmodels, alpha=0.05, niter=10^3, burnin=round(niter/10)) {

  if (!is.null(object$postmean) & (object$prior$priorCoef@priorDistr %in% c('bic','zellner'))) {

      postmean= object$postmean[1:min(maxmodels,nrow(object$postmean)),]
      postsd= sqrt(object$postvar[1:min(maxmodels,nrow(object$postvar)),])
      ans= list(postmean= postmean, ci.low= postmean + qnorm(alpha/2) * postsd, ci.up= postmean - qnorm(alpha/2) * postsd)
      nn= colnames(ans[[1]])
      
  } else {

    hasPostSampling(object)
    y= object$ystd; x= object$xstd
    outcometype= object$outcometype; family= object$family
    b= min(50, ceiling((burnin/niter) * niter))
    #List models for which estimates are to be obtained
    pp= postProb(object,method="norm")
    modelid= strsplit(as.character(pp$modelid), split=',')
    modelid= modelid[1:min(maxmodels,length(modelid))]
    priorCoef= object$priors$priorCoef
    priorGroup= object$priors$priorGroup
    priorVar= object$priors$priorVar
    #Obtain point estimates and posterior intervals
    ans= vector("list",3); names(ans)= c('postmean','ci.low','ci.up')
    if ((outcometype== 'Continuous') && (family== 'normal')) { ##Linear model
      ans[[1]]= ans[[2]]= ans[[3]]= matrix(0,nrow=length(modelid),ncol=ncol(x)+1)
      for (i in 1:length(modelid)) {  #for each model
        colsel= as.numeric(modelid[[i]])
        bm= coefOneModel(y=y, x=x[,colsel,drop=FALSE], outcometype=outcometype, family=family, priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, alpha=alpha, niter=niter, burnin=b)
        colselphi= c(colsel,ncol(ans[[1]]))
        ans[[1]][i,colselphi]= bm[,1]
        ans[[2]][i,colselphi]= bm[,2]
        ans[[3]][i,colselphi]= bm[,3]
      }
      if (is.null(colnames(x))) nn= c(paste('beta',1:ncol(x),sep=''),'phi') else nn= c(colnames(x),'phi')
    } else {                                                   ##GLM or Survival model
      ans[[1]]= ans[[2]]= ans[[3]]= matrix(0, nrow=length(modelid), ncol=ncol(x))
      for (i in 1:length(modelid)) {   #for each model
        colsel= as.numeric(modelid[[i]])
        if (length(colsel)>0) {
          bm= coefOneModel(y=y, x=x[,colsel,drop=FALSE], outcometype=outcometype, family=family, priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, alpha=alpha, niter=niter, burnin=b)
          ans[[1]][i,colsel]= bm[,1]
          ans[[2]][i,colsel]= bm[,2]
          ans[[3]][i,colsel]= bm[,3]
        }
      }
      if (is.null(colnames(x))) nn= paste('beta',1:ncol(x),sep='') else nn= colnames(x)
    }
    colnames(ans[[1]])= colnames(ans[[2]])= colnames(ans[[3]]) = nn
    rownames(ans[[1]])= rownames(ans[[2]])= rownames(ans[[3]]) = pp$modelid[1:nrow(ans[[1]])]

  }

  #Return parameter estimates in non-standardized parameterization
  ans[[1]]= unstdcoef(ans[[1]],p=ncol(ans[[1]]),msfit=object,coefnames=nn)
  ans[[2]]= unstdcoef(ans[[2]],p=ncol(ans[[2]]),msfit=object,coefnames=nn)
  ans[[3]]= unstdcoef(ans[[3]],p=ncol(ans[[3]]),msfit=object,coefnames=nn)
  return(ans)
}
)



setMethod("coefOneModel", signature(y='ANY',x='matrix',m='missing',V='missing',outcometype='character',family='character'), function(y, x, m, V, outcometype, family, priorCoef, priorGroup, priorVar, alpha=0.05, niter=10^3, burnin=round(niter/10)) {
  if ((outcometype== 'Continuous') && (family== 'normal')) {  ##Linear model
    b= rnlpLM(y=y, x=x, priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, niter=niter, burnin=burnin)
  } else if (outcometype=='glm') { #GLM
    b= rnlpGLM(y=y, x=x, family=family, priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, niter=niter, burnin=burnin)
  } else if ((outcometype=='Survival') && (family=='Cox')) {  #Cox model
    b= rnlpCox(y=y, x=x, priorCoef=priorCoef, priorGroup=priorGroup, niter=niter, burnin=burnin)
  } else {
    stop(paste("outcometype",outcometype,"and family=",family,"not implemented",sep=""))
  }
  ans= cbind(colMeans(b), t(apply(b,2,quantile,probs=c(alpha/2,1-alpha/2))))
  return(ans)
}
)



predict.msfit <- function(object, newdata, data, level=0.95, ...) {
    hasPostSampling(object)
    if (is.null(colnames(object$xstd))) colnames(object$xstd)= paste('x',1:ncol(object$xstd),sep='')
    th= rnlp(msfit=object,niter=10^4)
    mx= object$stdconstants[-1,'shift']; sx= object$stdconstants[-1,'scale']
    if (!missing(newdata)) {
        f= object$call$formula
        if ('formula' %in% class(f)) {
            if (missing(data)) stop("You must specify the argument 'data'")
            alldata= rbind(data,newdata)
            alldata[,as.character(f)[2]]= 0  #ensure there's no NAs in the response, so createDesign doesn't drop those rows from newdata
            nn= rownames(alldata)[(nrow(data)+1):(nrow(data)+nrow(newdata))]
            if (is.null(object$call$smoothterms)) {
                des= createDesign(f, data=alldata)
            } else {
                des= createDesign(f, data=alldata, smoothterms=object$call$smoothterms, splineDegree=object$call$splineDegree, nknots=object$call$nknots)
            }
            newdata= des$x[nn,,drop=FALSE]
        }
        #ct= (sx==0)
        #newdata[,!ct]= t((t(newdata[,!ct]) - mx[!ct])/sx[!ct])
        if (is.null(colnames(newdata))) colnames(newdata)= paste('x',1:ncol(newdata),sep='')
    } else {
        newdata= t(t(object$xstd) * sx + mx)
    }
    sel= colnames(th) %in% colnames(newdata)
    ypred= th[,sel] %*% t(newdata)
    if ("intercept" %in% colnames(th)[!sel]) ypred= ypred + th[,'intercept']
    ans= cbind(mean=colMeans(ypred), t(apply(ypred,2,quantile,probs=c((1-level)/2,1-(1-level)/2))))
    return(ans)
}


getmodelid= function(object) {
  if (!inherits(object, 'msfit')) stop("Function modelid requires an argument of type msfit")
  if (!is.null(object$models)) {
    ans= object$models[,c('modelid','family')]
  } else {
    modelid= apply(object$postSample==1, 1, function(z) paste(which(z),collapse=','))
    if (object$family=='auto') {
      twopiece <- laplace <- logical(length(modelid))
      twopiece[grep(as.character(object$p+1),modelid)] <- TRUE
      laplace[grep(as.character(object$p+2),modelid)] <- TRUE
      family <- character(length(modelid))
      family[(!twopiece) & (!laplace)] <- 'normal'
      family[twopiece & (!laplace)] <- 'twopiecenormal'
      family[(!twopiece) & laplace] <- 'laplace'
      family[twopiece & laplace] <- 'twopiecelaplace'
      modelid <- sub(paste(',',object$p+1,sep=''),'',modelid)
      modelid <- sub(as.character(object$p+1),'',modelid)  #for null model
      modelid <- sub(paste(',',object$p+2,sep=''),'',modelid)
      modelid <- sub(as.character(object$p+2),'',modelid)  #for null model
    } else {
      family= object$family
    }
    ans= data.frame(modelid=modelid, family=family)
  }
  return(ans)
}



setMethod("logjoint", signature(object='msfit'), function(object, return_models=TRUE) {
  if (object$enumerate) {
      if (return_models) {
        ans= data.frame(object$models, object$postProb)
      } else {
        ans= object$postProb
      }
  } else {
    ans= unique(data.frame(object$postSample==1, logpp=object$postProb))
    if (!return_models) ans= ans[,ncol(ans)]
  }
return(ans)
}
)


setMethod("postProb", signature(object='msfit'), function(object, nmax, method='norm') {
if (!is.null(object$models)) {
    ans= object$models
} else {
  if (method=='norm') {
    modelpp <- unique(data.frame(object$postSample==1, logpp=object$postProb))
    modelpp <- data.frame(modelid= apply(modelpp[,1:(ncol(modelpp)-1)], 1, function(z) paste(which(z),collapse=',')), logpp=modelpp$logpp)
    modelpp$logpp <- modelpp$logpp - max(modelpp$logpp)
    modelpp$pp <- exp(modelpp$logpp)/sum(exp(modelpp$logpp))
  } else if (method=='exact') {
    modelpp <- apply(object$postSample==1, 1, function(z) paste(which(z),collapse=','))
    modelpp <- table(modelpp)/length(modelpp)
    modelpp <- data.frame(modelid=names(modelpp), pp=as.numeric(modelpp))
  } else {
    stop("Argument 'method' not recognized")
  }
  modelpp <- modelpp[order(modelpp$pp,decreasing=TRUE),]
  if (!missing(nmax)) modelpp <- modelpp[1:nmax,]
  if (object$family=='auto') {
    modelid <- as.character(modelpp[,'modelid'])
    twopiece <- laplace <- logical(nrow(modelpp))
    twopiece[grep(as.character(object$p+1),modelid)] <- TRUE
    laplace[grep(as.character(object$p+2),modelid)] <- TRUE
    family <- character(nrow(modelpp))
    family[(!twopiece) & (!laplace)] <- 'normal'
    family[twopiece & (!laplace)] <- 'twopiecenormal'
    family[(!twopiece) & laplace] <- 'laplace'
    family[twopiece & laplace] <- 'twopiecelaplace'
    modelid <- sub(paste(',',object$p+1,sep=''),'',modelid)
    modelid <- sub(as.character(object$p+1),'',modelid)  #for null model
    modelid <- sub(paste(',',object$p+2,sep=''),'',modelid)
    modelid <- sub(as.character(object$p+2),'',modelid)  #for null model
    modelpp <- data.frame(modelid=modelid,family=family,pp=modelpp[,'pp'])
  } else {
    modelpp <- data.frame(modelid=modelpp[,'modelid'],family=object$family,pp=modelpp[,'pp'])
  }
  ans= modelpp[,c('modelid','family','pp')]
}
return(ans)
}
)


#Convert text model identifier (e.g. "1,3,4") into logical identifier (e.g. c(TRUE,FALSE,TRUE,TRUE))
modelid2logical= function(modelid, nvars) {
  modelid= strsplit(modelid, split=',')
  ans= matrix(FALSE, nrow=length(modelid), ncol=nvars)
  for (i in 1:nrow(ans)) {
      sel= as.numeric(modelid[[i]])
      ans[i, sel[sel<=nvars]]= TRUE
  }
  return(ans)
}


#Obtain posterior probability for subsets of variables indicated in varsubset (the remaining parameters are passed on to postProb)
# varsubset can take one of the following 4 formats
# - A logical vector indicating which variables are in the subset and length equal to the number of variables, e.g. c(TRUE, FALSE, TRUE, TRUE)
# - A logical matrix where each row indicates a subset, as in the previous entry
# - A numeric vector indicating the indices of the variables in the subset, e.g. c(1,3,4)
# - A list where each entry is a numeric vector as in the previous entry
setMethod("postProbSubset", signature(object='msfit'), function(object, varsubset, nmax, method='norm') {
    if (!is.list(varsubset) & !is.matrix(varsubset)) {
      if (is.logical(varsubset)) varsubset= matrix(varsubset, nrow=1) else varsubset= list(varsubset)
    }
    if (is.matrix(varsubset)) nsubsets= nrow(varsubset) else if (is.list(varsubset)) nsubsets= length(varsubset) else stop("varsubset has the wrong format")
    modelpp= postProb(object, nmax=nmax, method=method)
    modelid= modelid2logical(modelpp$modelid, nvars= ncol(object$xstd))
    colnames(modelid)= colnames(object$xstd)
    ans= double(nsubsets)
    if (is.matrix(varsubset)) {
      for (i in 1:nsubsets) {
          nselvars= rowSums(modelid[,varsubset[i,],drop=FALSE]) #number of variables in the subset selected by each model
          selmodel= (nselvars== sum(varsubset[i,]))  #models selecting all variables in the subset
          ans[i]= sum(modelpp$pp[selmodel])
      }
    } else {
      for (i in 1:nsubsets) {
          nselvars= rowSums(modelid[,varsubset[[i]],drop=FALSE])    #number of variables in the subset selected by each model
          selmodel= (nselvars== length(varsubset[[i]]))  #models selecting all variables in the subset
          ans[i]= sum(modelpp$pp[selmodel])
      }
    }
    return(ans)
}
)



defaultmom= function(outcometype,family) {
    if (outcometype=='Continuous') {
        cat("Using default prior for continuous outcomes priorCoef=momprior(tau=0.348), priorVar=igprior(.01,.01)\n")
        priorCoef= momprior(tau=0.348)
        priorVar= igprior(alpha=.01,lambda=.01)
    } else if (outcometype=='Survival') {
        cat("Using default prior for Normal AFT survival outcomes priorCoef=momprior(tau=0.192), priorVar=igprior(3,3)\n")
        priorCoef= momprior(tau=0.192)
        priorVar= igprior(alpha=3,lambda=3)
    } else if (outcometype=='glm') {
        cat("Using default prior for GLMs priorCoef=momprior(tau=1/3), priorVar=igprior(.01,.01)\n")
        priorCoef= momprior(tau=1/3)
        priorVar= igprior(alpha=.01,lambda=.01)
    } else {
      stop("There is not default priorCoef for this outcome type")
    }
    ans= list(priorCoef=priorCoef, priorVar=priorVar)
    return(ans)
}



#### General model selection routines
modelSelection <- function(y, x, data, smoothterms, nknots=9, groups=1:ncol(x), constraints, center=TRUE, scale=TRUE, enumerate, includevars=rep(FALSE,ncol(x)), models, maxvars, niter=5000, thinning=1, burnin=round(niter/10), family='normal', priorCoef, priorGroup, priorDelta=modelbbprior(1,1), priorConstraints, priorVar=igprior(.01,.01), priorSkew=momprior(tau=0.348), phi, deltaini=rep(FALSE,ncol(x)), initSearch='greedy', method='auto', adj.overdisp='intercept', hess='asymp', optimMethod, optim_maxit, initpar='none', B=10^5, XtXprecomp= ifelse(ncol(x)<10^4,TRUE,FALSE), verbose=TRUE) {
# Input
# - y: either formula with the regression equation or vector with response variable. If a formula arguments x, groups & constraints are ignored
# - x: design matrix with all potential predictors
# - data: data frame where the variables indicated in y (if it's a formula) and smoothterms can be found
# - smoothterms: formula indicating variables for which non-linear effects (splines) should be considered
# - nknots: number of knots
# - groups: vector indicating groups for columns in x (defaults to each variable in a separate group)
# - constraints: constraints on the model space. List with length equal to the number of groups; if group[[i]]=c(j,k) then group i can only be in the model if groups j and k are also in the model
# - center: if center==TRUE y and x are centered to have zero mean, therefore eliminating the need to include an intercept term in x.
# - scale: if scale==TRUE y and columns in x are scaled to have standard deviation 1
# - enumerate: if TRUE all models with up to maxvars are enumerated, else Gibbs sampling is used to explore the model space
# - includevars: set to TRUE for variables that you want to force into the model (for grouped variables, TRUE/FALSE is taken from 1st variable in each group)
# - maxvars: maximum number of variables in models to be enumerated (ignored if enumerate==FALSE)
# - niter: number of Gibbs sampling iterations
# - thinning: MCMC thinning factor, i.e. only one out of each thinning iterations are reported. Defaults to thinning=1, i.e. no thinning
# - burnin: number of burn-in MCMC iterations. Defaults to 10% of niter. Set to 0 for no burn-in.
# - family: assumed residual distribution ('normal','twopiecenormal','laplace','twopiecelaplace')
# - priorCoef: prior distribution for the coefficients. Must be object of class 'msPriorSpec' with slot priorType set to 'coefficients'. Possible values for slot priorDistr are 'pMOM', 'piMOM' and 'peMOM'.
# - priorGroup: prior on grouped coefficients, as indicated by groups
# - priorDelta: prior on model indicator space. Must be object of class 'msPriorSpec' with slot priorType set to 'modelIndicator'. Possible values for slot priorDistr are 'uniform' and 'binomial'
# - priorVar: prior on residual variance. Must be object of class 'msPriorSpec' with slot priorType set to 'nuisancePars'. Slot priorDistr must be equal to 'invgamma'.
# - priorSkew: prior on residual skewness parameter. Ignored unless family=='twopiecenormal' or 'twopiecelaplace'
# - phi: residual variance. Typically this is unknown and therefore left missing. If specified argument priorVar is ignored.
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# - initSearch: algorithm to refine deltaini. initSearch=='greedy' uses a greedy Gibbs sampling search. initSearch=='SCAD' sets deltaini to the non-zero elements in a SCAD fit with cross-validated regularization parameter. initSearch=='none' leaves deltaini unmodified.
# - method: method to compute marginal densities. method=='Laplace' for Laplace approx, method=='MC' for Importance Sampling, method=='Hybrid' for Hybrid Laplace-IS (the latter method is only used for piMOM prior with unknown residual variance phi), method='ALA' (former method=='plugin')
# - adj.overdisp: for method=='ALA' it indicates the over-dispersion adjustment to be made in models where the dispersion parameter is fixed, as in logistic and Poisson regression. adj.overdisp='none' for no adjustment (not recommended, particularly for Poisson models). adj.overdisp='intercept' to estimate over-dispersion from the intercept-only model. ad.overdisp='residuals' from the Pearson residuals of each model (slightly higher computational cost)
# - hess: only used for asymmetric Laplace errors. When hess=='asymp' the asymptotic hessian is used to compute the Laplace approximation to the marginal likelihood, when hess=='asympDiagAdj' a diagonal adjustment to the asymptotic Hessian is used
# - optimMethod: method to maximize objective function when method=='Laplace' or method=='MC'. Only used for family=='twopiecenormal'. optimMethod=='LMA' uses modified Newton-Raphson algorithm, 'CDA' coordinate descent algorithm
# - optim_maxit: maximum number of iterations in optimization method
# - initpar: initial value for regression parameters when finding the posterior mode to approximate the integrated likelihood. Either 'MLE', 'MLE-aisgd', 'L1', 'L2-aisgd', or a numeric vector specifying the initial values. If p<n/2 MLE is used, else L1 (regularization parameter set via BIC). If n>10,000 or p>200, then MLE-aisgd or L2-aisgd are used.
# - B: number of samples to use in Importance Sampling scheme. Ignored if method=='Laplace'.
# - verbose: set verbose==TRUE to print iteration progress
# Output: list
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)

  #Check input
  tmp <- formatInputdata(y=y,x=x,data=data,smoothterms=smoothterms,nknots=nknots,family=family)
  x <- tmp$x; y <- tmp$y; formula <- tmp$formula;
  splineDegree <- tmp$splineDegree
  if (!is.null(tmp$groups)) groups <- tmp$groups
  if (length(groups) != ncol(x)) stop(paste("groups has the wrong length. It should have length",ncol(x)))
  hasgroups <- tmp$hasgroups; isgroups <- tmp$isgroups
  if (!is.null(tmp$constraints)) constraints <- tmp$constraints
  outcometype <- tmp$outcometype; uncens <- tmp$uncens; ordery <- tmp$ordery
  typeofvar <- tmp$typeofvar
  call <- list(formula=formula, smoothterms= NULL, splineDegree=splineDegree, nknots=nknots)
  if (!missing(smoothterms)) call$smoothterms <- smoothterms
  p= ncol(x); n= length(y)
      if (is.numeric(includevars)) {
      tmp= rep(FALSE,p)
      if (max(includevars) > p) stop(paste("includevars contains index ",max(includevars)," but the design matrix only has ",p," columns",sep=""))
      tmp[includevars]= TRUE
      includevars= tmp
  }
  if (length(includevars)!=ncol(x) | (!is.logical(includevars))) stop("includevars must be a logical vector of length ncol(x)")
  if (missing(maxvars)) maxvars= ifelse(family=='auto', p+2, p)
  if (maxvars <= sum(includevars)) stop("maxvars must be >= sum(includevars)")

  #If there are variable groups, count variables in each group, indicate 1st variable in each group, convert group and constraint labels to integers 0,1,...
  if (missing(priorCoef)) {
      defaultprior= defaultmom(outcometype=outcometype,family=family)
      priorCoef= defaultprior$priorCoef; priorVar= defaultprior$priorVar
  }
  if (missing(priorGroup)) { if (length(groups)==length(unique(groups))) { priorGroup= priorCoef } else { priorGroup= groupzellnerprior(tau=n) } }
  tmp= codeGroupsAndConstraints(p=p,groups=groups,constraints=constraints)
  ngroups= tmp$ngroups; constraints= tmp$constraints; invconstraints= tmp$invconstraints; nvaringroup=tmp$nvaringroup; groups=tmp$groups
  if (missing(models)) {
    if (missing(enumerate)) enumerate= ifelse(ngroups<15,TRUE,FALSE)
  } else { enumerate= TRUE }

  #Standardize (y,x) to mean 0 and variance 1 (for continuous or log-survival time outcomes only)
  if (!is.vector(y)) { y <- as.double(as.vector(y)) } else { y <- as.double(y) }
  if (!is.matrix(x)) x <- as.matrix(x)
  mx= colMeans(x); sx= sqrt(colMeans(x^2) - mx^2) * sqrt(n/(n-1))
  ct= (sx==0)
  if (any(is.na(ct))) stop('x contains NAs, this is currently not supported, please remove the NAs')
  if (sum(ct)>1) stop('There are >1 constant columns in x (e.g. two intercepts)')
  if (!center) { my=0; mx= rep(0,p) } else { my= mean(y) }
  if (!scale) { sy=1; sx= rep(1,p) } else { sy= sd(y) }
  mx[typeofvar=='factor']=0; sx[typeofvar=='factor']= 1
  if (!(outcometype %in% c('Continuous','Survival'))) { my=0; sy= 1 }
  ystd= (y-my)/sy; xstd= x; xstd[,!ct]= t((t(x[,!ct]) - mx[!ct])/sx[!ct])
  if (missing(phi)) { knownphi <- as.integer(0); phi <- double(0) } else { knownphi <- as.integer(1); phi <- as.double(phi) }
  stdconstants= rbind(c(my,sy),cbind(mx,sx)); colnames(stdconstants)= c('shift','scale')

  #Format arguments for .Call
  if (missing(deltaini)) {
    deltaini= as.integer(which(includevars)-1); ndeltaini= as.integer(length(deltaini))
  } else {
    if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
    if (!is.logical(deltaini)) { stop('deltaini must be of type logical') } else { ndeltaini <- as.integer(sum(deltaini | includevars)); deltaini <- as.integer(which(deltaini | includevars)-1) }
  }

  thinit= getthinit(y=y, x=xstd, family=family, initpar=initpar, enumerate=enumerate)
  usethinit= thinit$usethinit; thinit= thinit$thinit

  method <-  formatmsMethod(method=method, usethinit=usethinit, optimMethod=optimMethod, optim_maxit=optim_maxit, priorCoef=priorCoef, priorGroup=priorGroup, knownphi=knownphi, outcometype=outcometype, family=family, hasgroups=hasgroups, adj.overdisp=adj.overdisp, hess=hess)
  optimMethod <- method$optimMethod; optim_maxit <- method$optim_maxit; adj.overdisp <- method$adj.overdisp; hesstype <- method$hesstype; method <- method$method

  niter <- as.integer(niter); burnin <- as.integer(burnin); thinning <- as.integer(thinning); B <- as.integer(B)
  sumy2 <- as.double(sum(ystd^2)); sumy <- as.double(sum(ystd)); ytX <- as.vector(matrix(ystd,nrow=1) %*% xstd); colsumsx <- as.double(colSums(xstd))
  if (XtXprecomp) {
      XtX= t(xstd) %*% xstd
      hasXtX= as.logical(TRUE)
  } else {
      XtX= double(0)
      hasXtX= as.logical(FALSE)
  }

  ffamily= formatFamily(family, issurvival= length(uncens)>0)
  familyint= ffamily$familyint; familygreedy= ffamily$familygreedy
  if (familyint == 22) { sumlogyfact= as.double(sum(lgamma(ystd+1))) } else { sumlogyfact= as.double(0) } #Poisson regression
    
  if (!is.null(colnames(xstd))) { nn <- colnames(xstd) } else { nn <- paste('x',1:ncol(xstd),sep='') }

  tmp= formatmsPriorsMarg(priorCoef=priorCoef, priorGroup=priorGroup, priorVar=priorVar, priorSkew=priorSkew, n=n)
  r= tmp$r; prior= tmp$prior; priorgr= tmp$priorgr; tau=tmp$tau; taugroup=tmp$taugroup; alpha=tmp$alpha; lambda=tmp$lambda; taualpha=tmp$taualpha; fixatanhalpha=tmp$fixatanhalpha
  priorCoef= tmp$priorCoef; priorGroup= tmp$priorGroup
    
  priorConstraints <- defaultpriorConstraints(priorDelta, priorConstraints)
  tmp= formatmsPriorsModel(priorDelta=priorDelta, priorConstraints=priorConstraints, constraints=constraints)
  prDelta=tmp$prDelta; prDeltap=tmp$prDeltap; parprDeltap=tmp$parprDeltap
  prConstr=tmp$prConstr; prConstrp= tmp$prConstrp; parprConstrp= tmp$parprConstrp

  #Run model selection
  if (!enumerate) {
    #Initialize
    includevars <- as.integer(includevars)
    if (familyint==0) { postMode <- rep(as.integer(0),p+2) } else { postMode <- rep(as.integer(0),p) }
    postModeProb <- double(1)
    if (initSearch=='greedy') {
      niterGreed <- as.integer(100)
      ans= greedyVarSelCI(knownphi,familygreedy,prior,priorgr,niterGreed,ndeltaini,deltaini,includevars,n,p,ystd,uncens,sumy2,sumy,sumlogyfact,xstd,colsumsx,hasXtX,XtX,ytX,method,adj.overdisp,hesstype,optimMethod,optim_maxit,thinit,usethinit,B,alpha,lambda,phi,tau,taugroup,taualpha,fixatanhalpha,r,prDelta,prDeltap,parprDeltap,prConstr,prConstrp,parprConstrp,groups,ngroups,nvaringroup,constraints,invconstraints,as.integer(verbose))
      postMode <- ans[[1]]; postModeProb <- ans[[2]]
      if (familyint==0) { postMode <- as.integer(c(postMode,0,0)); postModeProb <- as.double(postModeProb - 2*log(2)) }
      postMode[includevars==1] <- TRUE
      ndeltaini <- as.integer(sum(postMode)); deltaini <- as.integer(which(as.logical(postMode))-1)
    } else if (initSearch=='SCAD') {
      if (verbose) cat("Initializing via SCAD cross-validation...")
      deltaini <- rep(TRUE,ncol(xstd))
      cvscad <- cv.ncvreg(X=xstd[,!ct],y=ystd-mean(ystd),family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
      deltaini[!ct] <- ncvreg(X=xstd[,!ct],y=ystd-mean(ystd),penalty='SCAD',dfmax=1000,lambda=rep(cvscad$lambda[cvscad$cv],2))$beta[-1,1]!=0
      deltaini[includevars==1] <- TRUE
      ndeltaini <- as.integer(sum(deltaini)); deltaini <- as.integer(which(deltaini)-1)
      if (verbose) cat(" Done\n")
    }

    #Run MCMC
    ans <- modelSelectionGibbsCI(postMode,postModeProb,knownphi,familyint,prior,priorgr,niter,thinning,burnin,ndeltaini,deltaini,includevars,n,p,ystd,uncens,sumy2,sumy,sumlogyfact,as.double(xstd),colsumsx,hasXtX,XtX,ytX,method,adj.overdisp,hesstype,optimMethod,optim_maxit,thinit,usethinit,B,alpha,lambda,phi,tau,taugroup,taualpha,fixatanhalpha,r,prDelta,prDeltap,parprDeltap,prConstr,prConstrp,parprConstrp,groups,ngroups,nvaringroup,constraints,invconstraints,as.integer(verbose))
    postSample <- matrix(ans[[1]],ncol=ifelse(familyint!=0,p,p+2))
    margpp <- ans[[2]]; postMode <- ans[[3]]; postModeProb <- ans[[4]]; postProb <- ans[[5]]
    margpp[includevars==1]= 1
    postmean= postvar= NULL
    modelid= apply(postSample[,1:ncol(xstd),drop=FALSE]==1, 1, function(z) paste(which(z),collapse=','))

  } else {

    #Model enumeration
    if (verbose) cat("Enumerating models...\n")
    nincludevars= sum(includevars)
    nvars= ifelse(familyint==0,ncol(xstd)+2-nincludevars,ncol(xstd)-nincludevars)
    if (familyint==0) { includeenum= c(includevars[groups+1],FALSE,FALSE) } else { includeenum= includevars[groups+1] }
    if (missing(models)) {
      models= listmodels(vars2list=1:ngroups, includevars=includeenum, constraints=sapply(constraints,function(z) z+1), nvaringroup=nvaringroup, maxvars=maxvars) #listmodels expects group indexes 1,2,...
    } else {
      if (!is.logical(models)) stop("models must be a logical matrix")
      if (ncol(models) != ncol(xstd)) stop(paste("models has",ncol(models),"but x has",ncol(xstd),"columns"))
    }
    if (familyint==0) models= rbind(cbind(models,FALSE,FALSE),cbind(models,FALSE,TRUE),cbind(models,TRUE,FALSE),cbind(models,TRUE,TRUE))
    nmodels= as.integer(nrow(models))
    models= as.integer(models)
    includevars= as.integer(includevars)
    ans= modelSelectionEnumCI(nmodels,models,knownphi,familyint,prior,priorgr,n,p,ystd,uncens,sumy2,sumy,sumlogyfact,as.double(xstd),colsumsx,hasXtX,XtX,ytX,method,adj.overdisp,hesstype,optimMethod,optim_maxit,thinit,usethinit,B,alpha,lambda,phi,tau,taugroup,taualpha,fixatanhalpha,r,prDelta,prDeltap,parprDeltap,prConstr,prConstrp,parprConstrp,groups,ngroups,nvaringroup,constraints,invconstraints,as.integer(verbose))
    postMode <- ans[[1]]; postModeProb <- ans[[2]]; postProb <- ans[[3]]
    postSample <- matrix(nrow=0,ncol=ifelse(familyint!=0,p,p+2))
    models <- matrix(models,nrow=nmodels)
    pp <- exp(postProb-postModeProb); pp <- matrix(pp/sum(pp),ncol=1)
    margpp <- as.vector(t(models) %*% pp)
    modelid= apply(models[,1:ncol(xstd),drop=FALSE]==1, 1, function(z) paste(which(z),collapse=','))
    if (familyint==0) {
        modelfam= models[,ncol(xstd)+1] + 2*models[,ncol(xstd)+2]
        margpp= c(margpp[1:ncol(xstd)],sum(pp[modelfam==0]),sum(pp[modelfam==1]),sum(pp[modelfam==2]),sum(pp[modelfam==3]))
        modeltxt= ifelse(modelfam==0,'normal',ifelse(modelfam==1,'twopiecenormal',ifelse(modelfam==2,'laplace','twopiecelaplace')))
        models= data.frame(modelid=modelid,family=modeltxt,pp=pp)
        postmean= postvar= NULL
    } else {
        models= data.frame(modelid=modelid,family=family,pp=pp)
        postmean= postvar= NULL
    }
    modelid= models$modelid
    models= models[order(models$pp,decreasing=TRUE),]
  }

  #Post-process output
  if (familyint!=0) {
    colnames(postSample) <- names(postMode) <- names(margpp) <- nn
  } else {
    colnames(postSample) <- names(postMode)<- c(nn,'asymmetry','laplace')
    names(margpp) <- c(nn,'family.normal','family.tpnormal','family.laplace','family.tplaplace')
  }

  priors= list(priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, priorConstraints=priorConstraints, priorVar=priorVar, priorSkew=priorSkew)
  if (length(uncens)>0) { ystd[ordery]= ystd; uncens[ordery]= uncens; ystd= Surv(time=ystd, event= uncens); xstd[ordery,]= xstd }
  names(constraints)= paste('group',0:(length(constraints)-1))
  ans <- list(postSample=postSample,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb,modelid=modelid,postmean=postmean,postvar=postvar,family=family,p=ncol(xstd),enumerate=enumerate,priors=priors,ystd=ystd,xstd=xstd,groups=groups,constraints=constraints,stdconstants=stdconstants,outcometype=outcometype,call=call)
  if (enumerate) { ans$models= models }
  new("msfit",ans)
}

# format input data from either formula (y), formula and data.frame (y,data) or matrix and vector (y, x)
# it accepts smoothterms, groups and survival data
formatInputdata <- function(y,x,data,smoothterms,nknots,family) {
  valid_families <- c('normal','twopiecenormal','laplace','twopiecelaplace','auto','binomial','binomial logit','poisson','poisson log')
  if (!(family %in% valid_families)) stop(paste("Invalid family. Valid values are", valid_families))
  call <- match.call()
  groups <- NULL; constraints <- NULL; ordery <- NULL
  if ('formula' %in% class(y)) {
      formula= y; is_formula=TRUE; splineDegree= 3
      des= createDesign(y, data=data, smoothterms=smoothterms, splineDegree=splineDegree, nknots=nknots)
      x= des$x; groups= des$groups; constraints= des$constraints; typeofvar= des$typeofvar
      if ('Surv' %in% class(des$y)) {
          if (all(des$y[,1] >0)) {
              cat("Response type is survival and all its values are >0. Remember that you should log-transform the response prior to running modelSelection\n")
          }
          outcometype= 'Survival'; uncens= as.integer(des$y[,2]); y= des$y[,1]
          ordery= c(which(uncens==1),which(uncens!=1)); y= y[ordery]; x= x[ordery,,drop=FALSE]; uncens= uncens[ordery]
          if (family !="normal") stop("For survival outcomes only family='normal' is currently implemented")
      } else {
          if (family %in% c('normal','twopiecenormal','laplace','twopiecelaplace','auto')) {
            outcometype= 'Continuous'
          } else {
            outcometype= 'glm'
          }
          y= des$y; uncens= integer(0)
      }
      nlevels <- apply(x,2,function(z) length(unique(z)))
      typeofvar[nlevels==2]= 'factor'
  } else {
      if ('Surv' %in% class(y)) {
          outcometype= 'Survival'; uncens= as.integer(y[,2]); y= y[,1]
          ordery= c(which(uncens==1),which(uncens!=1)); y= y[ordery]; x= x[ordery,,drop=FALSE]; uncens= uncens[ordery]
          if (family !="normal") stop("For survival outcomes only family='normal' is currently implemented")
      } else {
        if (family %in% c('normal','twopiecenormal','laplace','twopiecelaplace','auto')) {
          outcometype= 'Continuous'
        } else {
          outcometype= 'glm'
        }
        uncens= integer(0)
      }
      formula= splineDegree= NA; is_formula=FALSE; typeofvar= rep('numeric',ncol(x))
  }
  if (nrow(x)!=length(y)) stop('nrow(x) must be equal to length(y)')
  if (any(is.na(y))) stop('y contains NAs, this is currently not supported, please remove the NAs')
  hasgroups <-  (length(groups) > length(unique(groups)))
  ningroup <- table(groups)
  isgroup <- (groups %in% as.numeric(names(ningroup)[ningroup>1]))
  y <- as.double(y)
  #Check that support of y is valid for the specified family
  if (family %in% c('binomial','binomial logit')) {
      if (any(!(y %in% c(0,1)))) stop("Invalid value for the response. For logistic regression it must be 0 or 1")
  } else if (family %in% c('poisson','poisson logit')) {
      if (any(y < 0) || any((y %% 1) != 0)) stop("Invalid value for the response. For Poisson regression it must be a natural number")
  }
  if (any(is.na(y)) | any(is.infinite(y))) stop("y contains NAs or infinite values")
  if (any(is.na(x)) | any(is.infinite(x))) stop("x contains NAs or infinite values")
  ans <- list(
    x=x, y=y, formula=formula, is_formula=is_formula, splineDegree=splineDegree,
    groups=groups, hasgroups=hasgroups, isgroup=isgroup, constraints=constraints, outcometype=outcometype, uncens=uncens,
    ordery=ordery, typeofvar=typeofvar
  )
  return(ans)
}



#Create a design matrix for the given formula. Return also covariate groups (e.g. from factors), covariate type (factor/numeric) and hierarchical constraints (e.g. from interaction terms), these are the parameters "groups" and "constraints" in modelSelection
createDesign <- function(formula, data, smoothterms, subset, na.action, splineDegree=3, nknots=14) {
    call <- match.call()
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    if (!missing(subset)) {
        m <- match(c("formula", "data", "subset"), names(mf), 0L)
    } else {
        m <- match(c("formula", "data"), names(mf), 0L)
    }
    mf <- mf[c(1L, m)]
    mf$na.action = quote(na.pass)
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    #gam.slist <- gam.smoothers()$slist
    mt <- if (missing(data)) terms(formula) else terms(formula, data = data)
    #mt <- if (missing(data)) terms(formula, gam.slist) else terms(formula, gam.slist, data = data)
    mf$formula <- mt
    mf <- eval(mf, parent.frame())
    if (missing(na.action)) {
        naa = getOption("na.action", "na.fail")
        na.action = get(naa)
    }
    mf = na.action(mf)
    mt = attributes(mf)[["terms"]]  #mt is an object of class "terms" storing info about the model, see help(terms.object) for a description
    y <- model.response(mf, "any")
    x <- if (!is.empty.model(mt)) model.matrix(mt, mf) else matrix(, NROW(y), 0)
    #x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(, NROW(y), 0)
    groups <- attr(x,"assign") #group that each variable belongs to, e.g. for factors
    tab= table(groups);
    typeofvar= ifelse(groups %in% as.numeric(names(tab)[tab>1]),'factor','numeric')
    intercept= ifelse(min(groups)==0,1,0)
    groups2vars= attr(mt,"factors")[-1,] #for each variable group, hierarchical dependence on other groups
    nn= colnames(groups2vars)[!(colnames(groups2vars) %in% rownames(groups2vars))]
    if (length(nn)>0) { #there's interaction terms
        tmp= matrix(0,nrow=length(nn),ncol=ncol(groups2vars))
        rownames(tmp)= nn
        colnames(tmp)= colnames(groups2vars)
        for (i in 1:nrow(tmp)) {
            nni= paste(strsplit(nn[i],split=":")[[1]], collapse=":.*") #regular expression, e.g. if nn[i]="Xj:Xl" it checks for "Xj:.*Xl"
            tmp[i,grep(nni,colnames(tmp))]= 1
        }
        groups2vars= rbind(groups2vars,tmp)
        constraints= lapply(1:ncol(groups2vars), function(i) { ans= as.integer(which(groups2vars[,i]>0)); ans[ans!=i] + intercept })
        if (intercept==1) constraints= c(list(integer(0)),constraints)
    } else { #there are no interactions
        constraints= lapply(1:(max(groups)+intercept), function(i) integer(0))
    }
    groups= groups+intercept
    #Add spline terms
    if (!missing(smoothterms)) {
        if (!any(c('formula','matrix','data.frame') %in% class(smoothterms))) stop("smoothterms should be of class 'formula', 'matrix' or 'data.frame'")
        Snested= nestedSplines(x=x, groups=groups, smoothterms=smoothterms, data=data, subset=subset, na.action=na.action, splineDegree=splineDegree, nknots=nknots)
        x= cbind(x, Snested$L, Snested$W)
        groups= c(groups, Snested$groups, Snested$groupsW)
        typeofvar= c(typeofvar, Snested$typeofvar)
        constraints= c(constraints, Snested$constraintsL, Snested$constraintsW)
        #old code, before structuring nestedSplines as a separate function
        #x= cbind(x,L[,selL],W)
        #groups= c(groups, groupsL[selL], groupsW)
        #typeofvar= c(typeofvar, rep('numeric',sum(selL)+length(groupsW)))
        #constraints= c(constraints, constraintsL[selL], constraintsW)
    }
    return(list(y=y,x=x,groups=groups,constraints=constraints,typeofvar=typeofvar))
}


#Create design matrix for nested splines: linear within knots[1], knots[1] within knots[2], etc.
nestedSplines= function(x, groups, smoothterms, data, subset, na.action, splineDegree, nknots) {
    if (length(nknots)>1) nknots= nknots[order(nknots)]
    maxgroups= max(groups)
    if ('formula' %in% class(smoothterms)) {
        smoothterms= formula(paste("~ ",-1,"+",as.character(smoothterms)[2])) #remove intercept
        L= createDesign(smoothterms, data=data, subset=subset, na.action=na.action)$x
    } else {
        L= as.matrix(smoothterms)
        if (is.null(colnames(L))) colnames(L)= paste("L",1:ncol(L),sep="")
    }

    W= constraintsW= constraintsL= groupsW= groupsL= vector("list", length(nknots))
    for (kk in 1:length(nknots)) {
        W[[kk]] <- matrix(NA,nrow=nrow(L),ncol=(nknots[kk]-4)*ncol(L))
        namesW <- character(ncol(W[[kk]]))
        groupsW[[kk]] <- integer(ncol(W[[kk]])); constraintsW[[kk]]= vector("list",ncol(L))
        groupsL[[kk]] <- integer(ncol(L)); constraintsL[[kk]]= lapply(1:ncol(L), function(i) integer(0))
        selL= rep(FALSE,ncol(L)); nselL= 0
        for (j in 1:ncol(L)) {
            m= range(L[,j])
            tmp= bspline(L[,j], degree=splineDegree, knots=seq(m[1],m[2],length=nknots[kk])) #equally-spaced knots
            Lj= cbind(1,L[,j]); b= solve(t(Lj) %*% Lj) %*% (t(Lj) %*% tmp)
            idx= (1+(j-1)*(nknots[kk]-4)):(j*(nknots[kk]-4))
            W[[kk]][,idx]= tmp - Lj %*% b #project splines onto space orthogonal to linear term
            namesW[idx]= paste(colnames(L)[j],'.s',1:length(idx),sep='')
            repeated= which(colnames(x)==colnames(L)[j])
            if (length(repeated)==1) {  #linear term was already in x
                constraintsW[[kk]][[j]]= groups[repeated]
            } else {                    #linear term wasn't in x
                selL[j]= TRUE
                nselL= nselL+1
                groupsL[[kk]][j]= maxgroups + nselL
                constraintsW[[kk]][[j]]= groupsL[[kk]][j]
            }
            groupsW[[kk]][idx]= j
        }
        colnames(W[[kk]])= namesW
        groupsW[[kk]]= maxgroups + nselL + groupsW[[kk]]
        varW= (colMeans(W[[kk]]^2) - colMeans(W[[kk]])^2) * nrow(W[[kk]])/(nrow(W[[kk]])-1)
        if (any(varW<1.0e-4)) {
            W[[kk]]= W[[kk]][,varW>1.0e-4]; groupsW[[kk]]= groupsW[[kk]][varW>1.0e-4]; constraintsW[[kk]]= constraintsW[[kk]][unique(groupsW[[kk]]-min(groupsW[[kk]])+1)]
        }
    }

    #Combine design matrices hierarchically into single matrix with hierarchical constraints
    Wall= W[[1]]; groupsWall= groupsW[[1]]; constraintsWall= constraintsW[[1]]
    maxgroups= max(groupsWall)
    if (length(nknots)==1) {
        ans= list(L=L[,selL], W=Wall, groupsL=groupsL[selL], groupsW=groupsWall, typeofvar= rep('numeric',sum(selL)+length(groupsW)), constraintsL=constraintsL[selL], constraintsW=constraintsWall)
    } else {
        groupsWidx= unique(groupsW[[1]])
        for (kk in 2:length(nknots)) {
            keep= vector("list",length(groupsWidx))
            for (g in 1:length(keep)) {
                sel1= which(groupsW[[kk]] == groupsWidx[g])
                sel2= which(groupsWall == groupsWidx[g])
                r= cor(W[[kk]][,sel1], Wall[,sel2])  #correlation with basis with previous number of knots
                o= order(abs(apply(r, 1, max)))      #sort by max absolute correlation
                keep[[g]]= sel1[o[1:(length(sel1) - sum(groupsW[[kk-1]] == groupsWidx[g]))]]
                keep[[g]]= keep[[g]][order(keep[[g]])]
            }
            newgroups= groupsW[[kk]][unlist(keep)]; newgroups= newgroups + maxgroups - min(newgroups) + 1
            maxgroups= max(newgroups)
            newconstraintsW= as.list(unique(newgroups) - length(groupsWidx))
            Wall= cbind(Wall, W[[kk]][,unlist(keep)])
            groupsWall= c(groupsWall, newgroups)
            constraintsWall= c(constraintsWall, newconstraintsW)
        }
        for (j in 1:ncol(L)) {
            pattern= sub("\\[","\\\\[",colnames(L)[j]); pattern= sub("\\]","\\\\]",pattern)
            selj= grep(pattern, colnames(Wall))
            colnames(Wall)[selj]= paste(colnames(L)[j],'.s',1:length(selj), sep='')
        }
        ans= list(L=L[,selL], W=Wall, groupsL=groupsL[selL], groupsW=groupsWall, typeofvar= rep('numeric',sum(selL)+length(groupsWall)), constraintsL=constraintsL[selL], constraintsW=constraintsWall)
    }

    return(ans)
}





#Count variables in each group, indicate 1st variable in each group, convert group and constraint labels to integers 1,2,...
codeGroupsAndConstraints= function(p,groups,constraints) {
    groupsnum= as.numeric(factor(groups)); groupsnum= cumsum(c(TRUE,groupsnum[-1]!=groupsnum[-length(groupsnum)])) #re-code groups (in order of appearance)
    ngroups= max(groupsnum)
    if (ngroups>p) stop("There cannot be more groups than variables (columns in x)")
    if (missing(constraints)) {
        constraints= sapply(1:ngroups, function(i) integer(0))
    } else {
        if (length(constraints) != ngroups) stop("length(constraints) must be equal to number of variable groups")
        #Ensure that constraints match the group order in groupsnum
        g2code= sapply(names(constraints), function(nn) groupsnum[match(TRUE,groups==nn)])
        if (any(names(constraints) != as.character(g2code))) {
            constraints= constraints[g2code]
            names(constraints)= g2code
            for (i in 1:length(constraints)) { if (length(constraints[[i]]>0)) constraints[[i]]= as.integer(g2code[constraints[[i]]]) -as.integer(1) } #group codes start at 0
        } else {
            constraints= lapply(constraints,function(z) as.integer(z)-as.integer(1)) #group codes start at 0
        }
    }
    if (ngroups==p) {
        nvaringroup= as.integer(rep(1,p))
        groups= as.integer(0:(p-1))
    } else {
        nvaringroup= as.integer(table(groupsnum)) #number of variables in each group
        groups= as.integer(groupsnum-1) #group id that each variable belongs to
        #groups= as.integer(c(0,as.numeric(cumsum(nvaringroup[-length(nvaringroup)])))) #1st variable in each group (0-indexed)
    }
    #Determine inverse constraints
    invconstraints= vector("list",ngroups)
    tabconstr= cbind(group=rep(0:(ngroups-1), sapply(constraints,length)), requires= unlist(constraints))
    for (i in 1:ngroups) { invconstraints[[i]]= as.integer(tabconstr[tabconstr[,'requires']==(i-1), 'group']) }
    ans= list(ngroups=ngroups,constraints=constraints,invconstraints=invconstraints,nvaringroup=nvaringroup,groups=groups)
    return(ans)
}


#Routine to enumerate all models satisfying hierarchical constraints, e.g. x[i] can only be in model if x[j] and x[k] are in model
#
# - vars2list: vector indicating variable groups, e.g. 1:10 means there's 10 variable groups named 1-10
# - includevars: logical vector of length= length(vars2list). TRUE indicates that all models should include that variable group
# - constraints: list with length= length(vars2list). Each element indicates hierarchical restrictions. Restrictions must be given in order, i.e. a restriction on group i can only depend on groups < i
# - nvaringroup: number of variables in each group
# - fixedvars: for internal use only. When calling the function recursively, fixedvars are the variables that are currently included in the model
#
# OUTPUT: matrix with models in rows and variables in columns. If the (i,j) entry is TRUE, model i includes variable j
#
# EXAMPLE 1: list models under restriction that x3 included only when (x1,x2) also included
#
# listmodels(vars2list=1:3, constraints=list(integer(0),integer(0),c(1,2)))
#
# EXAMPLE 2: same but forcing inclusion of x2
#
# listmodels(vars2list=1:3, includevars= c(FALSE,TRUE,FALSE), constraints=list(integer(0),integer(0),c(1,2)))
#
# EXAMPLE 3: list models under restriction that group 3 included only when (group 1, group 2) also included
#
# listmodels(vars2list=1:3, constraints=list(integer(0),integer(0),c(1,2)), nvaringroup= c(1,3,2))
#
# EXAMPLE 4: list models under restrictions x4 requires (x1,x2), x5 requires (x1,x3), x6 requires (x2,x3), x7 requires (x1,...,x6)
#
# listmodels(vars2list=1:7, constraints=list(integer(0),integer(0),integer(0),c(1,2),c(1,3),c(2,3),1:6))

listmodels= function(vars2list, includevars=rep(FALSE,length(vars2list)), fixedvars=integer(0), constraints, nvaringroup=rep(1,length(vars2list)), maxvars) {
    var1= vars2list[1]
    if (includevars[var1]) { #forcing inclusion of this variable group
        if (length(vars2list)>1) {
            ans= listmodels(vars2list[-1], includevars=includevars, fixedvars=c(fixedvars,var1), constraints=constraints, nvaringroup=nvaringroup, maxvars=maxvars)
        } else {
            fixedLogical= rep(FALSE,length(constraints)-1)
            fixedLogical[fixedvars]= TRUE
            fixedLogical= rep(fixedLogical,nvaringroup[-length(nvaringroup)])
            ans= c(fixedLogical,rep(TRUE,nvaringroup[length(nvaringroup)]))
        }
    } else { #not forcing inclusion of this variable group
        nfixedvars= length(fixedvars)
        if (length(constraints[[var1]])>0) {
            var1constraint= all(constraints[[var1]] %in% fixedvars) #is constraint on var1 satisfied by fixedvars?
        } else { var1constraint= TRUE }
        if (length(vars2list)>1) { #if this is not the last variable, call function recursively
            ans1= listmodels(vars2list[-1], includevars=includevars, fixedvars=fixedvars, constraints=constraints, nvaringroup=nvaringroup, maxvars=maxvars)
            if (var1constraint & (nfixedvars<maxvars)) {
                ans2= listmodels(vars2list[-1], includevars=includevars, fixedvars=c(fixedvars,var1), constraints=constraints, nvaringroup=nvaringroup, maxvars=maxvars)
            } else { ans2= NULL }
            ans= rbind(ans1,ans2)
        } else {  #if this is the last variable, return entire variable inclusion vector
            fixedLogical= rep(FALSE,length(constraints)-1)
            fixedLogical[fixedvars]= TRUE
            fixedLogical= rep(fixedLogical,nvaringroup[-length(nvaringroup)])
            if (var1constraint & (nfixedvars<maxvars)) {
                ans= rbind(c(fixedLogical,rep(FALSE,nvaringroup[length(nvaringroup)])),c(fixedLogical,rep(TRUE,nvaringroup[length(nvaringroup)])))
            } else {
                ans= c(fixedLogical, rep(FALSE,nvaringroup[length(nvaringroup)]))
            }
        }
    }
    return(ans)
}


#Routine to format method indicating how integrated likelihoods should be computed in modelSelection
#
#For any method other than 'auto', it sets a numeric code corresponding to that method to pass onto C.
#for method=='auto', the integration method is decided automatically according to the following rules
#
# 1. Continuous outcomes
# - if family normal, no groups: for pMOM exact / ALA; for other NLPs & LPs: Laplace
# - if family normal, groups:  for pMOMgMOM and pMOMgZell ALA; for other NLPs & LPs: Laplace
# - if family != normal: Laplace approximation (since ALA not currently implemented)
#
# 2. Survival outcomes. Use ALA when available, LA otherwise. Currently ALA available for pmom/groupMOM/groupzellner + pmom/groupMOM/groupzellner
#
# 3. GLMs. Laplace approximation
#
# Note: usethinit==3 indicates to use the initial parameter values in LA and ALA.
#
# OUTPUT
#
#   method
#   0 means Laplace approximation (note: for normal outcomes + normal priors this means the calculation is exact)
#   1 means Monte Carlo (Importance Sampling)
#   2 means ALA (approximate Laplace approximation)
#   -1 only available for pMOM, it means to use exact calculation for small models (<=3 parameters) and ALA for larger models
#  
#  optimMethod: 1 means Newton-Raphson type algorithm, 2 means Coordinate Descent Algorith, 0 means use the default. This argument is used in twopiecenormal models and GLMs, else it's ignored
#
#  optim_maxit: maximum number of iterations for optimization method. If not specified on input, a -1 is returned (so the C code uses a default maxit)
#
#  hesstype: what type of hessian to use in twopiecelaplace models, where the observed hessian is not defined
#
#  adj.overdisp
#  0: no over-dispersion adjustment
#  1: estimate over-dispersion from intercept-only model
#  2: estimate over-dispersion from Pearson residuals under each model
formatmsMethod= function(method, usethinit, initpar, optimMethod, optim_maxit=optim_maxit, priorCoef, priorGroup, knownphi, outcometype, family, hasgroups, adj.overdisp, hess) {
  hesstype <- as.integer(ifelse(hess=='asympDiagAdj',2,1))

  #Obtain code for the optimization method
  if (missing(optimMethod)) optimMethod <- 'auto'
  if (outcometype == 'glm') {
    if (optimMethod=='auto') {
        optimMethod <- as.integer(0)
    } else if (optimMethod %in% c('Newton','LMA')) {
        optimMethod <- as.integer(1)
    } else if (optimMethod == 'CDA') {
        optimMethod <- as.integer(2)
    } else { stop("Invalid optimMethod. For this family, only 'auto', 'Newton', 'LMA' or 'CDA' are implemented") }
  } else {
    if (optimMethod %in% c('Newton','LMA')) {
        optimMethod <- as.integer(1)
    } else {
        optimMethod <- as.integer(2)
    }
  }

  #Obtain code to be passed to C for the method to compute the integrated likelihood
  if ((method=='ALA') && (usethinit==3)) {
      method= 'Laplace'  #ALA init at theta != 0 is implemented in C as Laplace with limited optim_maxit
      if (missing(optim_maxit)) { optim_maxit= as.integer(0) } else { optim_maxit= as.integer(optim_maxit) } #maxit>0 tells C to use specified maxit
  } else {
      if (missing(optim_maxit)) { optim_maxit= as.integer(-1) } else { optim_maxit= as.integer(optim_maxit) } #maxit= -1 tells C to use its own defaults
  }
  
  if (method=='Laplace') {
    method <- as.integer(0)
  } else if (method=='MC') {
    method <- as.integer(1)
  } else if (method=='Hybrid') {
    if ((priorCoef@priorDistr!='piMOM') | (knownphi==1)) {
      warning("method=='Hybrid' is only available for 'piMOM' priors with unknown phi. Using method=='Laplace' instead")
      method <- as.integer(0)
    } else {
      method <- as.integer(2)
    }
  } else if (method=='ALA') {
      method <- as.integer(2)
  } else if (method=='auto') {
    if (outcometype=='Continuous') {
      if (family=='normal') {
          if (!hasgroups) {
            if (priorCoef@priorDistr=='pMOM') { method <- as.integer(-1) } else { method <- as.integer(0) }
          } else {
            if ((priorCoef@priorDistr=='pMOM') & (priorGroup@priorDistr %in% c('groupMOM','zellner','groupzellner'))) {
                method <- as.integer(2)
            } else { method <- as.integer(0) }
          }
      } else {
          method <- as.integer(0)
      }
    } else if (outcometype=='Survival') {
      if (family=='normal') {
        if ((priorCoef@priorDistr %in% c('pMOM','groupMOM','groupzellner')) & (priorGroup@priorDistr %in% c('pMOM','groupMOM','groupzellner'))) {
           method <- as.integer(2)
        } else {
           method <- as.integer(0)
        }
      } else { stop("For survival outcomes, only family=='normal' currently implemented") }
    } else if (outcometype=='glm') {
        method <- as.integer(0)
    } else { stop("outcometype must be 'Continuous', 'Survival' or 'glm'") }
  } else if ((method=='ALA') | (method=='plugin')) {
    method <- as.integer(2)
  } else {
    stop("Invalid 'method'")
  }

  adj.overdisp <- as.integer(ifelse(adj.overdisp=='none',0,ifelse(adj.overdisp=='intercept',1,2)))
  ans <- list(method=method, optimMethod=optimMethod, optim_maxit=optim_maxit, adj.overdisp=adj.overdisp, hesstype= hesstype)
  return(ans)
}

#Routine to format modelSelection prior distribution parameters for marginal likelihood
#Input: priorCoef, priorVar, priorGroup, priorSkew
#Output: parameters for prior on coefficients (r, prior, tau), prior on variance parameter (alpha, lambda), skewness parameter (taualpha, fixatanhalpha)
formatmsPriorsMarg <- function(priorCoef, priorGroup, priorVar, priorSkew, n) {
  r= as.integer(1)

  if (priorCoef@priorDistr=='bic') {
      
    prior <- priorgr <- as.integer(100)
    tau <- as.double(priorCoef@priorPars['tau'])
    taugroup <- alpha <- lambda <- taualpha <- fixatanhalpha <- as.double(-1)
      
  } else {
  
    if (missing(priorGroup)) priorGroup= priorCoef
    if (missing(priorSkew)) priorSkew= priorCoef
    has_taustd <- "taustd" %in% names(priorCoef@priorPars)
    has_taugroupstd <- "taustd" %in% names(priorGroup@priorPars)
    if (has_taustd) {
      taustd <- as.double(priorCoef@priorPars['taustd'])
    } else {
      tau <- as.double(priorCoef@priorPars['tau'])
    }
    if (has_taugroupstd) {
      taugroupstd <- as.double(priorGroup@priorPars['taustd'])
    } else {
      taugroup <- as.double(priorGroup@priorPars['tau'])
    }
    if (priorCoef@priorDistr=='pMOM') {
      r <- as.integer(priorCoef@priorPars['r'])
      prior <- as.integer(0)
      if (has_taustd) tau <- taustd * 1/3
    } else if (priorCoef@priorDistr=='piMOM') {
      prior <- as.integer(1)
    } else if (priorCoef@priorDistr=='peMOM') {
      prior <- as.integer(2)
    } else if (priorCoef@priorDistr=='zellner') {
      prior <- as.integer(3)
      if (has_taustd) tau <- taustd * n
    } else if (priorCoef@priorDistr=='normalid') {
      prior <- as.integer(4)
      if (has_taustd) tau <- taustd
    } else if (priorCoef@priorDistr=='groupMOM') {
      prior <- as.integer(10)
      if (has_taustd) tau <- taustd
    } else if (priorCoef@priorDistr=='groupzellner') {
      prior <- as.integer(13)
      if (has_taustd) tau <- taustd * n
    } else {
      stop('Prior specified in priorDistr not recognized')
    }
    if (priorGroup@priorDistr=='pMOM') {
      priorgr= as.integer(0)
      if (has_taugroupstd) taugroup <- taugroupstd * 1/3
    } else if (priorGroup@priorDistr=='piMOM') {
      priorgr= as.integer(1)
    } else if (priorGroup@priorDistr=='peMOM') {
      priorgr= as.integer(2)
    } else if (priorGroup@priorDistr=='zellner') {
      priorgr= as.integer(3)
      if (has_taugroupstd) taugroup <- taugroupstd * n
    } else if (priorGroup@priorDistr=='normalid') {
      priorgr= as.integer(4)
      if (has_taugroupstd) taugroup <- taugroupstd
    } else if (priorGroup@priorDistr=='groupMOM') {
      priorgr= as.integer(10)
      if (has_taugroupstd) taugroup <- taugroupstd
    } else if (priorGroup@priorDistr=='groupiMOM') {
      priorgr= as.integer(11)
    } else if (priorGroup@priorDistr=='groupeMOM') {
      priorgr= as.integer(12)
    } else if (priorGroup@priorDistr=='groupzellner') {
      priorgr= as.integer(13)
      if (has_taugroupstd) taugroup <- taugroupstd * n
    } else {
      stop('Prior in priorGroup not recognized')
    }
    alpha <- as.double(priorVar@priorPars['alpha']); lambda <- as.double(priorVar@priorPars['lambda'])
    #
    if ('msPriorSpec' %in% class(priorSkew)) {
        taualpha <- as.double(priorSkew@priorPars['tau'])
        fixatanhalpha <- as.double(-10000)
    } else {
        taualpha <- 0.358
        fixatanhalpha <- as.double(priorSkew)
    }
    if (has_taustd) priorCoef@priorPars['tau']= tau
    if (has_taugroupstd) priorGroup@priorPars['tau']= taugroup
  }
  
  ans= list(r=r,prior=prior,priorgr=priorgr,tau=tau,taugroup=taugroup,alpha=alpha,lambda=lambda,taualpha=taualpha,fixatanhalpha=fixatanhalpha,priorCoef=priorCoef,priorGroup=priorGroup)
  return(ans)
}

defaultpriorConstraints <- function(priorDelta, priorConstraints) {
  if (missing(priorConstraints)) {
    if ((priorDelta@priorDistr=='binomial') && ('p' %in% names(priorDelta@priorPars)) && (length(priorDelta@priorPars[['p']]) > 1)) {
      priorConstraints <- modelbinomprior(p=0.5)
    } else {
      priorConstraints <- priorDelta
    }
  }
  return(priorConstraints)
}

#Routine to format modelSelection prior distribution parameters in model space
#Input: priorDelta, priorConstraints, constraints
#Output: model space prior (prDelta, prDeltap, parprDeltap) and constraints (prConstr,prConstrp,parprConstrp)
formatmsPriorsModel <- function(priorDelta, priorConstraints, constraints) {
  #Prior on model space (parameters not subject to hierarchical constraints)
  n_unconstrained <- sum(sapply(constraints, function(x) length(x) == 0))
  n_constrained <- length(constraints) - n_unconstrained
  if (priorDelta@priorDistr=='uniform') {
    prDelta <- as.integer(0)
    prDeltap <- as.double(0)
    parprDeltap <- double(2)
  } else if (priorDelta@priorDistr=='binomial') {
    if ('p' %in% names(priorDelta@priorPars)) {
      prDelta <- as.integer(1)
      prDeltap <- as.double(priorDelta@priorPars[['p']])
      if (any(prDeltap<=0) | any(prDeltap>1)) stop("For the binomial model prior the inclusion probabilities p must lie in (0,1]")
      if ((length(prDeltap) != 1) & (length(prDeltap) != n_unconstrained)) stop("p in priorDelta must be a scalar or have length=number of unconstrained variables")
      parprDeltap <- as.double(length(prDeltap))
    } else {
      prDelta <- as.integer(2)
      prDeltap <- as.double(.5)
      parprDeltap <- as.double(priorDelta@priorPars[c('alpha.p','beta.p')])
    }
  } else if (priorDelta@priorDistr=='complexity') {
      prDelta <- as.integer(3)
      prDeltap <- as.double(priorDelta@priorPars['c'])
      if (prDeltap<0) stop("c must be >0 for priorDelta@priorDistr=='complexity'")
      parprDeltap <- double(2)
  } else {
    stop('Prior specified in priorDelta not recognized')
  }
  #Prior on model space (parameters subject to hierarchical constraints)
  if (priorConstraints@priorDistr=='uniform') {
    prConstr <- as.integer(0)
    prConstrp <- as.double(0)
    parprConstrp <- double(2)
  } else if (priorConstraints@priorDistr=='binomial') {
    if ('p' %in% names(priorConstraints@priorPars)) {
      prConstr <- as.integer(1)
      prConstrp <- as.double(priorConstraints@priorPars[['p']])
      if (any(prConstrp<=0) | any(prConstrp>=1)) stop("p must be between 0 and 1 for priorConstraints@priorDistr=='binomial'")
      if ((length(prConstrp) != 1) & (length(prConstrp) != n_constrained)) stop("p in priorConstraints must be a scalar or have length=number of constrained variables")
      parprConstrp <- as.double(length(prConstrp))
    } else {
      prConstr <- as.integer(2)
      prConstrp <- as.double(.5)
      parprConstrp <- as.double(priorConstraints@priorPars[c('alpha.p','beta.p')])
    }
  } else if (priorConstraints@priorDistr=='complexity') {
      prConstr <- as.integer(3)
      prConstrp <- as.double(priorConstraints@priorPars['c'])
      if (prConstrp<0) stop("c must be >0 for priorConstraints@priorDistr=='complexity'")
      parprConstrp <- double(2)
  } else {
    stop('Prior specified in priorConstraints not recognized')
  }

  ans= list(prDelta=prDelta,prDeltap=prDeltap,parprDeltap=parprDeltap,prConstr=prConstr,prConstrp=prConstrp,parprConstrp=parprConstrp)
  return(ans)
}


#Assign a numerical code to the family of likelihoods, to pass on to C
formatFamily= function(family, issurvival) {
    
    if (family=='auto') {
        familyint= 0; familygreedy=1
    } else if (family=='normal') {
        familyint= familygreedy= ifelse(!issurvival,1,11)
    } else if (family=='twopiecenormal') {
        familyint= 2; familygreedy=1
    } else if (family=='laplace') {
        familyint= 3; familygreedy=1
    } else if (family=='twopiecelaplace') {
        familyint= 4; familygreedy=1
    } else if (family %in% c('binomial','binomial logit')) {
        familyint= familygreedy= 21
    } else if (family %in% c('poisson','poisson log')) {
        familyint= familygreedy= 22
    } else stop("family not available")

    ans= list(familyint= as.integer(familyint), familygreedy= as.integer(familygreedy))
    return(ans)
      
}



greedymodelSelectionR <- function(y, x, niter=100, marginalFunction, priorFunction, betaBinPrior, deltaini=rep(FALSE,ncol(x)), verbose=TRUE, ...) {
  #Greedy version of modelSelectionR where variables with prob>0.5 at current iteration are included deterministically (prob<.5 excluded)
  p <- ncol(x)
  if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
  if (!missing(betaBinPrior)) {
    #Initialize probBin
    if ((betaBinPrior['alpha.p']>1) && (betaBinPrior['beta.p']>1)) {
      probBin <- (betaBinPrior['alpha.p']-1)/(betaBinPrior['alpha.p']+betaBinPrior['beta.p']-2)
    } else {
      probBin <- (betaBinPrior['alpha.p'])/(betaBinPrior['alpha.p']+betaBinPrior['beta.p'])
    }
    postOther <- matrix(NA,nrow=niter,ncol=1); colnames(postOther) <- c('probBin')
    priorFunction <- function(sel, logscale=TRUE) dbinom(x=sum(sel),size=length(sel),prob=probBin,log=logscale)
  } else {
    postOther <- matrix(NA,nrow=niter,ncol=0)
  }
  #Greedy iterations
  sel <- deltaini
  mcur <- marginalFunction(y=y,x=x[,sel,drop=FALSE],logscale=TRUE,...) + priorFunction(sel,logscale=TRUE)
  nchanges <- 1; itcur <- 1
  nn <- names(x)
  while (nchanges>0 & itcur<niter) {
    nchanges <- 0; itcur <- itcur+1
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      mnew <- marginalFunction(y=y,x=x[,selnew,drop=FALSE],logscale=TRUE,...) + priorFunction(selnew,logscale=TRUE)
      if (mnew>mcur) { sel[i]=selnew[i]; mcur=mnew; nchanges=nchanges+1; if (verbose) cat(paste(ifelse(sel[i],"Added","Dropped"),nn[i],"\n",collapse=" ")) }
    }
  }
  return(sel)
}

#Gibbs model selection using BIC to approximate marginal likelihood
# - y, x, xadj: response, covariates under selection and adjustment covariates
# - family: glm family, passed on to glm
# - niter: number of Gibbs iteration
# - burnin: burn-in iterations
# - modelPrior: function evaluating model log-prior probability. Takes a logical vector as input
# Returns:
# - postModel, postCoef1, postCoef2, margpp (analogous to pmomPM. postCoef1 & postCoef2 are MLEs under each visited model)
modelselBIC <- function(y, x, xadj, family, niter=1000, burnin= round(.1*niter), modelPrior, verbose=TRUE) {
  pluginJoint <- function(sel) {
    p <- sum(sel)
    ans <- vector("list",2); names(ans) <- c('marginal','coef')
    if (p>0 & p<=length(y)) {
      glm1 <- glm(y ~ x[,sel,drop=FALSE] + xadj -1, family=family)
      ans$marginal <- -.5*glm1$deviance - .5*log(length(y))*(glm1$df.null-glm1$df.residual) + modelPrior(sel)
      ans$coef1 <- coef(glm1)[1:p]; ans$coef2 <- coef(glm1)[-1:-p]
    } else if (p==0) {
      glm1 <- glm(y ~ xadj -1, family=family)
      ans$marginal <- -.5*glm1$deviance - .5*log(length(y))*(glm1$df.null-glm1$df.residual) + modelPrior(sel)
      ans$coef1 <- double(0); ans$coef2 <- coef(glm1)
    } else { ans$marginal <- -Inf; ans$coef1 <- ans$coef2 <- 0}
    return(ans)
  }
  #Greedy iterations
  if (verbose) cat("Initializing...")
  sel <- rep(FALSE,ncol(x))
  mcur <- pluginJoint(sel)$marginal
  nchanges <- 1; it <- 1
  while (nchanges>0 & it<100) {
    nchanges <- 0; it <- it+1
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      mnew <- pluginJoint(selnew)$marginal
      if (mnew>mcur) { sel[i] <- selnew[i]; mcur <- mnew; nchanges <- nchanges+1 }
    }
  }
  if (verbose) { cat(" Done\nGibbs sampling") }
  #Gibbs iterations
  niter10 <- ceiling(niter/10)
  postModel <- matrix(NA,nrow=niter,ncol=ncol(x))
  postCoef1 <- matrix(0,nrow=niter,ncol=ncol(x))
  postCoef2 <- matrix(0,nrow=niter,ncol=ncol(xadj))
  curmod <- pluginJoint(sel)
  for (j in 1:niter) {
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      newmod <- pluginJoint(selnew)
      mnew <- newmod$marginal
      if (runif(1) < 1/(1+exp(mcur-mnew))) { sel[i] <- selnew[i]; mcur <- mnew; curmod <- newmod }
    }
    postModel[j,] <- sel
    postCoef1[j,sel] <- curmod$coef1; postCoef2[j,] <- curmod$coef2
    if (verbose & ((j%%niter10)==0)) cat(".")
  }
  if (verbose) cat("Done\n")
  #Return output
  if (burnin>0) { postModel <- postModel[-1:-burnin,,drop=FALSE]; postCoef1 <- postCoef1[-1:-burnin,,drop=FALSE]; postCoef2 <- postCoef2[-1:-burnin,,drop=FALSE] }
  ans <- list(postModel=postModel, postCoef1=postCoef1, postCoef2=postCoef2, margpp=colMeans(postModel))
}


## Common prior distributions on model space

nselConstraints= function(sel, groups, constraints) {
    #Output
    # - ngroups0: number of unconstrained groups that are selected
    # - ngroups1: number of constrained groups that are selected
    # - ngroups: total number of groups (selected or unselected)
    # - ngroupsconstr: total number of groups that have a constraint
    # - violateConstraint: TRUE if sel violates a group constraint, FALSE otherwise
    ngroups= length(unique(groups))
    if (length(constraints) != ngroups) stop("length(constraints) must be equal to length(unique(groups))")
    if (!is.numeric(groups)) stop("Argument groups must be numeric")
    nconstraints= sapply(constraints,length)
    hasconstraint= (nconstraints>0)
    ngroupsconstr= sum(hasconstraint)
    tab= table(factor(sel,levels=c('FALSE','TRUE')),groups)
    selgroup= tab['TRUE',] >0
    violateConstraint= FALSE
    if (ngroupsconstr>0) { for (i in which(hasconstraint)) { if (selgroup[i] & any(!selgroup[constraints[[i]]])) violateConstraint= TRUE } }
    ngroups0= sum(selgroup[!hasconstraint]); ngroups1= sum(selgroup[hasconstraint])
    return(list(ngroups0=ngroups0, ngroups1=ngroups1, ngroups=ngroups, ngroupsconstr=ngroupsconstr, violateConstraint=violateConstraint, hasconstraint=hasconstraint))
}

#binomPrior <- function(sel, prob=.5, logscale=TRUE) {  dbinom(x=sum(sel),size=length(sel),prob=prob,log=logscale) }

binomPrior <- function(sel, prob=.5, logscale=TRUE, probconstr=prob, groups=1:length(sel), constraints=lapply(1:length(unique(groups)), function(z) integer(0))) {
    nsel= nselConstraints(sel=sel, groups=groups, constraints=constraints)
    if (!nsel$violateConstraint) {
      hasconstraint <- nsel$hasconstraint
      ans <- sum((log(prob) * sel)[!hasconstraint]) + sum((log(1-prob)*(1-sel))[!hasconstraint])
      if (nsel$ngroupsconstr>0) ans= ans+ sum((log(probconstr) * sel)[hasconstraint]) + sum((log(1-probconstr)*(1-sel))[hasconstraint])
    } else { ans= -Inf }
    if (!logscale) ans= exp(ans)
    return(ans)
}


#unifPrior <- function(sel, logscale=TRUE) { ifelse(logscale,-length(sel)*log(2),2^(-length(sel)))  }
unifPrior <- function(sel, logscale=TRUE, groups=1:length(sel), constraints=lapply(1:length(unique(groups)), function(z) integer(0))) {
    binomPrior(sel, prob=.5, logscale=logscale, probconstr=.5, groups=groups, constraints=constraints)
}


#bbPrior <- function(sel, alpha=1, beta=1, logscale=TRUE) {
#  ans <- lbeta(sum(sel) + alpha, sum(!sel) + beta) - lbeta(alpha,beta)
#  ifelse(logscale,ans,exp(ans))
#}

bbPrior <- function(sel, alpha=1, beta=1, logscale=TRUE, alphaconstr=alpha, betaconstr=beta, groups=1:length(sel), constraints=lapply(1:length(unique(groups)), function(z) integer(0))) {
  nsel= nselConstraints(sel=sel, groups=groups, constraints=constraints)
  if (!nsel$violateConstraint) {
      ans= lbeta(nsel$ngroups0 + alpha, nsel$ngroups-nsel$ngroupsconstr-nsel$ngroups0 + beta) - lbeta(alpha,beta)
      if (nsel$ngroupsconstr>0) ans= ans + lbeta(nsel$ngroups1 + alphaconstr, nsel$ngroupsconstr - nsel$ngroups1 + betaconstr) - lbeta(alphaconstr,betaconstr)
  } else { ans= -Inf }
  ifelse(logscale,ans,exp(ans))
}




bbPriorTrunc <- function (sel, logscale=TRUE, maxvars=10, ...) {
  #Same as bbPrior with prob=0 when variables > maxvars
  if (sum(sel)<=maxvars) { ans <- bbPrior(sel, logscale=logscale, ...) } else { ans <- ifelse(logscale, -Inf, 0) }
  return(ans)
}

unifPriorTrunc <- function (sel, logscale=TRUE, maxvars=10, ...) {
  #Same as unifPrior with prob=0 when variables > maxvars
  if (sum(sel)<=maxvars) { ans <- unifPrior(sel, logscale=logscale, ...) } else { ans <- ifelse(logscale, -Inf, 0) }
  return(ans)
}
