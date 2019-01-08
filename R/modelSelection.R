###
### modelSelection.R
###

### Methods for msfit objects

setMethod("show", signature(object='msfit'), function(object) {
  cat('msfit object with',object$p,'variables and',object$family,'residual distribution\n')
  ifelse(any(object$postMode!=0), paste('  Posterior mode: covariate',which(object$postMode==1)), '  Posterior mode: null model')
  cat("Use postProb() to get posterior model probabilities\n")
  cat("Use coef() to get BMA parameter estimates, 95% intervals and marginal inclusion probabilities\n")
  cat("Elements $margpp, $postMode, $postSample and $coef contain further information (see help('msfit') and help('modelSelection') for details)\n")
}
)


#setMethod("bmaCoef", signature(y='ANY',x='matrix',msfit='msfit'), function(y, x, msfit, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
#
#}
#)

coef.msfit <- function(object,...) {
    th= rnlp(msfit=object,niter=10^4)
    ct= (object$stdconstants[-1,'sd']==0)
    if (any(ct)) {
        margpp= c(object$margpp,1)
        nn= c(names(object$margpp),'phi')
    } else {
        margpp= c(mean(th[,1]!=0),object$margpp,1)
        nn= c('intercept',names(object$margpp),'phi')
    }
    ans= cbind(colMeans(th),t(apply(th,2,quantile,probs=c(.025,0.975))),margpp=margpp)
    colnames(ans)= c('estimate','2.5%','97.5%','margpp')
    rownames(ans)= nn
    return(ans)
}


setMethod("postProb", signature(object='msfit'), function(object, nmax, method='norm') {
if (!is.null(object$models)) {
    ans= object$models
} else {
  if (method=='norm') {
    modelpp <- unique(data.frame(object$postSample==1, logpp=object$postProb))
    modelpp <- data.frame(modelid= apply(modelpp[,1:(ncol(modelpp)-1)], 1, function(z) paste(which(z),collapse=',')), logpp=modelpp$logpp)
    modelpp$logpp <- modelpp$logpp - modelpp$logpp[1]
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



#### General model selection routines
modelSelection <- function(y, x, data, groups=1:ncol(x), constraints, center=TRUE, scale=TRUE, enumerate= ifelse(ncol(x)<15,TRUE,FALSE), includevars=rep(FALSE,ncol(x)), maxvars, niter=10^4, thinning=1, burnin=round(niter/10), family='normal', priorCoef=momprior(tau=0.348), priorDelta=modelbbprior(alpha.p=1,beta.p=1), priorVar=igprior(alpha=.01,lambda=.01), priorSkew=momprior(tau=0.348), phi, deltaini=rep(FALSE,ncol(x)), initSearch='greedy', method='auto', hess='asymp', optimMethod='CDA', B=10^5, verbose=TRUE) {
# Input
# - y: vector with response variable
# - x: design matrix with all potential predictors
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
# - priorDelta: prior on model indicator space. Must be object of class 'msPriorSpec' with slot priorType set to 'modelIndicator'. Possible values for slot priorDistr are 'uniform' and 'binomial'
# - priorVar: prior on residual variance. Must be object of class 'msPriorSpec' with slot priorType set to 'nuisancePars'. Slot priorDistr must be equal to 'invgamma'.
# - priorSkew: prior on residual skewness parameter. Ignored unless family=='twopiecenormal' or 'twopiecelaplace'
# - phi: residual variance. Typically this is unknown and therefore left missing. If specified argument priorVar is ignored.
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# - initSearch: algorithm to refine deltaini. initSearch=='greedy' uses a greedy Gibbs sampling search. initSearch=='SCAD' sets deltaini to the non-zero elements in a SCAD fit with cross-validated regularization parameter. initSearch=='none' leaves deltaini unmodified.
# - method: method to compute marginal densities. method=='Laplace' for Laplace approx, method=='MC' for Importance Sampling, method=='Hybrid' for Hybrid Laplace-IS (the latter method is only used for piMOM prior with unknown residual variance phi), method='plugin'
# - hess: only used for asymmetric Laplace errors. When hess=='asymp' the asymptotic hessian is used to compute the Laplace approximation to the marginal likelihood, when hess=='asympDiagAdj' a diagonal adjustment to the asymptotic Hessian is used
# - optimMethod: method to maximize objective function when method=='Laplace' or method=='MC'. Only used for family=='twopiecenormal'. optimMethod=='LMA' uses modified Newton-Raphson algorithm, 'CDA' coordinate descent algorithm
# - B: number of samples to use in Importance Sampling scheme. Ignored if method=='Laplace'.
# - verbose: set verbose==TRUE to print iteration progress
# Output: list
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)

  #Check input
  if (class(y)=="formula") {
      des= createDesign(y, data=data)
      y= des$y; x= des$x; groups= des$groups; constraints= des$constraints
  }
  p= ncol(x); n= length(y)
  if (nrow(x)!=length(y)) stop('nrow(x) must be equal to length(y)')
  if (any(is.na(y))) stop('y contains NAs, this is currently not supported, please remove the NAs')
  if (length(includevars)!=ncol(x) | (!is.logical(includevars))) stop("includevars must be a logical vector of length ncol(x)")
  if (missing(maxvars)) maxvars= ifelse(family=='auto', p+2, p)
  if (maxvars <= sum(includevars)) stop("maxvars must be >= sum(includevars)")
    
  #If there are variable groups, count variables in each group, indicate 1st variable in each group, convert group and constraint labels to integers 0,1,...
  tmp= codeGroupsAndConstraints(p=p,groups=groups,constraints=constraints)
  ngroups= tmp$ngroups; constraints= tmp$constraints; nvaringroup=tmp$nvaringroup; groups=tmp$groups
    
  #Standardize (y,x) to mean 0 and variance 1
  if (!is.vector(y)) { y <- as.double(as.vector(y)) } else { y <- as.double(y) }
  if (!is.matrix(x)) x <- as.matrix(x)
  mx= colMeans(x); sx= sqrt(colMeans(x^2) - mx^2) * sqrt(n/(n-1))
  ct= (sx==0)
  if (any(is.na(ct))) stop('x contains NAs, this is currently not supported, please remove the NAs')
  if (sum(ct)>1) stop('There are >1 constant columns in x (e.g. two intercepts)')
  if (!center) { my=0; mx= rep(0,p) } else { my= mean(y) }
  if (!scale) { sy=1; sx= rep(1,p) } else { sy= sd(y) }
  ystd= (y-my)/sy; xstd= x; xstd[,!ct]= t((t(x[,!ct]) - mx[!ct])/sx[!ct])
  if (missing(phi)) { knownphi <- as.integer(0); phi <- double(0) } else { knownphi <- as.integer(1); phi <- as.double(phi) }
  stdconstants= rbind(c(my,sy),cbind(mx,sx)); colnames(stdconstants)= c('mean','sd')

  #Format arguments for .Call
  if (missing(deltaini)) {
    deltaini= as.integer(which(includevars)); ndeltaini= as.integer(length(deltaini))
  } else {
    if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
    if (!is.logical(deltaini)) { stop('deltaini must be of type logical') } else { ndeltaini <- as.integer(sum(deltaini | includevars)); deltaini <- as.integer(which(deltaini | includevars)-1) }
  }

  method= formatmsMethod(method=method, priorCoef=priorCoef, knownphi=knownphi)
  hess <- as.integer(ifelse(hess=='asympDiagAdj',2,1))
  optimMethod <- as.integer(ifelse(optimMethod=='CDA',2,1))
    
  niter <- as.integer(niter); burnin <- as.integer(burnin); thinning <- as.integer(thinning); B <- as.integer(B)
  sumy2 <- as.double(sum(ystd^2)); XtX <- t(xstd) %*% xstd; ytX <- as.vector(matrix(ystd,nrow=1) %*% xstd)

  tmp= formatmsPriors(priorCoef=priorCoef, priorVar=priorVar, priorSkew=priorSkew, priorDelta=priorDelta)
  r= tmp$r; prior= tmp$prior; tau=tmp$tau; alpha=tmp$alpha; lambda=tmp$lambda; taualpha=tmp$taualpha; fixatanhalpha=tmp$fixatanhalpha;
  prDelta=tmp$prDelta; prDeltap=tmp$prDeltap; parprDeltap=tmp$parprDeltap
    
  if (family=='auto') { familyint <- 0 } else if (family=='normal') { familyint <- 1 } else if (family=='twopiecenormal') { familyint <- 2 } else if (family=='laplace') { familyint <- 3 } else if (family=='twopiecelaplace') { familyint <- 4 } else stop("family not available")
  familyint <- as.integer(familyint)
  if (!is.null(colnames(xstd))) { nn <- colnames(xstd) } else { nn <- paste('x',1:ncol(xstd),sep='') }

  #Run model selection
  if (!enumerate) {
    #Initialize
    includevars <- as.integer(includevars)
    if (familyint==0) { postMode <- rep(as.integer(0),p+2) } else { postMode <- rep(as.integer(0),p) }
    postModeProb <- double(1)
    if (initSearch=='greedy') {
      niterGreed <- as.integer(100)
      ans= .Call("greedyVarSelCI",knownphi,prior,niterGreed,ndeltaini,deltaini,includevars,n,p,ystd,sumy2,xstd,XtX,ytX,method,hess,optimMethod,B,alpha,lambda,phi,tau,taualpha,fixatanhalpha,r,prDelta,prDeltap,parprDeltap,ngroups,nvaringroup,constraints,as.integer(verbose))
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
    ans <- .Call("modelSelectionGibbsCI", postMode,postModeProb,knownphi,familyint,prior,niter,thinning,burnin,ndeltaini,deltaini,includevars,n,p,ystd,sumy2,as.double(xstd),XtX,ytX,method,hess,optimMethod,B,alpha,lambda,phi,tau,taualpha,fixatanhalpha,r,prDelta,prDeltap,parprDeltap,ngroups,nvaringroup,constraints,as.integer(verbose))
    postSample <- matrix(ans[[1]],ncol=ifelse(familyint!=0,p,p+2))
    margpp <- ans[[2]]; postMode <- ans[[3]]; postModeProb <- ans[[4]]; postProb <- ans[[5]]

  } else {

    #Model enumeration
    if (verbose) cat("Enumerating models...\n")
    nincludevars= sum(includevars)
    nvars= ifelse(familyint==0,ncol(xstd)+2-nincludevars,ncol(xstd)-nincludevars)
    if (familyint==0) { includeenum= c(includevars[groups+1],FALSE,FALSE) } else { includeenum= includevars[groups+1] }
    models= listmodels(vars2list=1:ngroups, includevars=includeenum, constraints=sapply(constraints,function(z) z+1), nvaringroup=nvaringroup, maxvars=maxvars) #listmodels expects group indexes 1,2,...
    if (familyint==0) models= rbind(cbind(models,FALSE,FALSE),cbind(models,FALSE,TRUE),cbind(models,TRUE,FALSE),cbind(models,TRUE,TRUE))
    nmodels= as.integer(nrow(models))
    models= as.integer(models)
    includevars= as.integer(includevars)
    ans= .Call("modelSelectionEnumCI", nmodels,models,knownphi,familyint,prior,n,p,ystd,sumy2,as.double(xstd),XtX,ytX,method,hess,optimMethod,B,alpha,lambda,phi,tau,taualpha,fixatanhalpha,r,prDelta,prDeltap,parprDeltap,ngroups,nvaringroup,as.integer(verbose))
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
    } else {
        models= data.frame(modelid=modelid,family=family,pp=pp)
    }
    models= models[order(models$pp,decreasing=TRUE),]
  }

  #Post-process output
  if (familyint!=0) {
    colnames(postSample) <- names(postMode) <- names(margpp) <- nn
  } else {
    colnames(postSample) <- names(postMode)<- c(nn,'asymmetry','laplace')
    names(margpp) <- c(nn,'family.normal','family.tpnormal','family.laplace','family.tplaplace')
  }

  priors= list(priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar, priorSkew=priorSkew)
  ans <- list(postSample=postSample,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb,family=family,p=ncol(xstd),enumerate=enumerate,priors=priors,ystd=ystd,xstd=xstd,stdconstants=stdconstants)
  if (enumerate) { ans$models= models }
  new("msfit",ans)
}




#Create a design matrix for the given formula. Return also variable groups (e.g. from factors) and hierarchical constraints (e.g. from interaction terms), these are the parameters "groups" and "constraints" in modelSelection
createDesign <- function(formula, data, subset, na.action) {
    call <- match.call()
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
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
    x <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0)
    groups <- attr(x,"assign") #group that each variable belongs to, e.g. for factors
    intercept <- ifelse(min(groups)==0,1,0)
    groups2vars <- attr(mt,"factors")[-1,] #for each variable group, hierarchical dependence on other groups
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
    return(list(y=y,x=x,groups=groups,constraints=constraints))
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
        groups= c(0,as.numeric(cumsum(nvaringroup[-length(nvaringroup)]))) #1st variable in each group (0-indexed)
    }
    ans= list(ngroups=ngroups,constraints=constraints,nvaringroup=nvaringroup,groups=groups)
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
formatmsMethod= function(method, priorCoef, knownphi) {
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
  } else if (method=='auto') {
    if (priorCoef@priorDistr=='pMOM') { method <- as.integer(-1) } else { method <- as.integer(0) }
  } else if (method=='plugin') {
    method <- as.integer(2)
  } else {
    stop("Invalid 'method'")
  }
  return(method)
}

#Routine to format modelSelection prior distribution parameters
#Input: priorCoef, priorVar, priorSkew, priorDelta
#Output: parameters for prior on coefficients (r, prior, tau), prior on variance parameter (alpha, lambda), skewness parameter (taualpha, fixatanhalpha), model space prior (prDelta, prDeltap, parprDeltap)
formatmsPriors= function(priorCoef, priorVar, priorSkew, priorDelta) {
  if (priorCoef@priorDistr=='pMOM') {
    r <- as.integer(priorCoef@priorPars['r']); prior <- as.integer(0)
  } else if (priorCoef@priorDistr=='piMOM') {
    r <- as.integer(1); prior <- as.integer(1)
  } else if (priorCoef@priorDistr=='peMOM') {
    r <- as.integer(1); prior <- as.integer(2)
  } else if (priorCoef@priorDistr=='zellner') {
    r <- as.integer(1); prior <- as.integer(3)
  } else {
    stop('Prior specified in priorDistr not recognized')
  }
  tau <- as.double(priorCoef@priorPars['tau'])
  alpha <- as.double(priorVar@priorPars['alpha']); lambda <- as.double(priorVar@priorPars['lambda'])
  #
  if (class(priorSkew)=='msPriorSpec') {
      taualpha <- as.double(priorSkew@priorPars['tau'])
      fixatanhalpha <- as.double(-10000)
  } else {
      taualpha <- 0.358
      fixatanhalpha <- as.double(priorSkew)
  }
  #
  if (priorDelta@priorDistr=='uniform') {
    prDelta <- as.integer(0)
    prDeltap <- as.double(0)
    parprDeltap <- double(2)
  } else if (priorDelta@priorDistr=='binomial') {
    if ('p' %in% names(priorDelta@priorPars)) {
      prDelta <- as.integer(1)
      prDeltap <- as.double(priorDelta@priorPars['p'])
      if ((prDeltap<=0) | (prDeltap>=1)) stop("p must be between 0 and 1 for priorDelta@priorDistr=='binomial'")
      parprDeltap <- double(2)
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
  ans= list(r=r,prior=prior,tau=tau,alpha=alpha,lambda=lambda,taualpha=taualpha,fixatanhalpha=fixatanhalpha,prDelta=prDelta,prDeltap=prDeltap,parprDeltap=parprDeltap)
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
  #browser()
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

modelSelectionR <- function(y, x, niter=10^4, marginalFunction, priorFunction, betaBinPrior, deltaini=rep(FALSE,ncol(x)), verbose=TRUE, ...) {
# Input
# - y: vector with response variable
# - x: design matrix with all potential predictors
# - niter: number of Gibbs sampling iterations
# - marginalFunction: function to compute the marginal density of the data under each model
# - priorFunction: function to compute the model prior probability
# - betaBinPrior: if specified, priorFunction argument is ignored and set to a binomial prior with Beta-hyperprior for the success prob. betaBinPrior should be a vector with named elements 'alpha.p' and 'beta.p', e.g. betaBinPrior= c(alpha.p=1,beta.p=1)
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# ...: other arguments to be passed on to marginalFunction
# Output: list
# - postSample: posterior samples for model indicator
# - postOther: posterior samples for other parameters (probBin: success probability for Binomial prior on number of coefficients in the model)
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)
if (class(y)=='Surv') {
  if ((length(y)/2) != nrow(x)) stop("Dimensions of y and x do not match")
} else {
  if (length(y) != nrow(x)) stop("Dimensions of y and x do not match")
}
if (any(is.na(x))) stop("x cannot have missing values")
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
if (ncol(x)>12) {
  sel <- postMode <- deltaini
  currentJ <- postModeProb <- marginalFunction(y=y,x=x[,sel,drop=FALSE],logscale=TRUE,...) + priorFunction(sel,logscale=TRUE)
  postSample <- matrix(NA,nrow=niter,ncol=p)
  margpp <- double(p)
  postProb <- double(niter)
  k <- 1; postProb[k] <- postModeProb
  names(postProb)[k] <- paste("V",which(sel),collapse=',',sep='')
  niter10 <- ceiling(niter/10)
  for (i in 1:niter) {
    for (j in 1:p) {
      selnew <- sel; selnew[j] <- !sel[j]
      namenew <- paste(which(selnew),collapse=',')
      newJ <- postProb[namenew]
      if (is.na(newJ)) {
        newJ <- marginalFunction(y=y,x=x[,selnew,drop=FALSE],logscale=TRUE,...) + priorFunction(selnew,logscale=TRUE)
        k <- k+1; postProb[k] <- newJ
        names(postProb)[k] <- paste("V",which(selnew),collapse=',',sep='')
      }
      if (newJ>postModeProb) {
        postModeProb <- newJ
        postMode <- selnew
      }
      pp <- 1/(1+exp(-currentJ+newJ))
      if (sel[j]) {  #if variable in the model
        sel[j] <- runif(1)<pp
        margpp[j] <- margpp[j]+pp
      } else {       #if variable not in the model
        sel[j] <- runif(1)>pp
        margpp[j] <- margpp[j]+1-pp
      }
      if (sel[j]==selnew[j]) {  #if value was updated, update marginal and prior densities
        currentJ <- newJ
      }
    }
    if (!missing(betaBinPrior)) {
      probBin <- rbeta(1,betaBinPrior['alpha.p']+sum(sel), betaBinPrior['beta.p']+sum(!sel))
      postOther[i,'probBin'] <- probBin
    }
    postSample[i,] <- sel
    postProb[i] <- currentJ
    if (verbose && ((i%%niter10)==0)) cat('.')
  }
  margpp <- margpp/niter
  if (verbose) cat('\n')
  #Format postProb
  modelid <- sapply(apply(postSample,1,which),function(z) paste("V",z,collapse=',',sep=''))
  postProb <- postProb[modelid]
} else {
  if (verbose) cat(paste("Computing posterior probabilities for all",2^ncol(x),"models..."))
  models <- expand.grid(lapply(1:ncol(x),function(z) c(FALSE,TRUE)))
  postProb <- apply(models,1, function(z) marginalFunction(y=y,x=x[,z,drop=FALSE],logscale=TRUE,...) + priorFunction(z,logscale=TRUE))
  if (verbose) cat('\n')
  postMode <- models[which.max(postProb),]
  postModeProb <- max(postProb)
  pp <- postProb - postModeProb; pp <- exp(pp)/sum(exp(pp))
  sampledmodels <- rep(1:nrow(models), rmultinom(1,size=niter,prob=pp)[,1])
  postSample <- models[sampledmodels,]
  postProb <- postProb[sampledmodels]
  margpp <- as.vector(t(models) %*% matrix(pp,ncol=1))
}
ans <- list(postSample=postSample,postOther=postOther,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb)
ans <- new("msfit",ans)
return(ans)
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

binomPrior <- function(sel, prob=.5, logscale=TRUE) {  dbinom(x=sum(sel),size=length(sel),prob=prob,log=logscale) }
unifPrior <- function(sel, logscale=TRUE) { ifelse(logscale,-length(sel)*log(2),2^(-length(sel)))  }
bbPrior <- function(sel, alpha=1, beta=1, logscale=TRUE) {
  ans <- lbeta(sum(sel) + alpha, sum(!sel) + beta) - lbeta(alpha,beta)
  ifelse(logscale,ans,exp(ans))
}


bbPriorTrunc <- function (sel, logscale=TRUE, maxvars=10) {
  #Same as bbPrior with prob=0 when variables > maxvars
  if (sum(sel)<=maxvars) { ans <- bbPrior(sel, logscale=logscale) } else { ans <- ifelse(logscale, -Inf, 0) }
  return(ans)
}

unifPriorTrunc <- function (sel, logscale=TRUE, maxvars=10) {
  #Same as unifPrior with prob=0 when variables > maxvars
  if (sum(sel)<=maxvars) { ans <- unifPrior(sel, logscale=logscale) } else { ans <- ifelse(logscale, -Inf, 0) }
  return(ans)
}
