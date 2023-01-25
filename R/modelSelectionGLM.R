modelSelectionGLM= function(y, x, data, smoothterms, nknots=9, groups=1:ncol(x), constraints, center=TRUE, scale=TRUE, includevars=rep(FALSE,ncol(x)), maxvars, familyglm= gaussian(), priorCoef, priorGroup, priorVar, priorDelta= modelbbprior(1,1), models, verbose= TRUE) {
    #Enumerate all or specified models and return posterior prob as approximated by BIC. Enumerated models satisfy group and hierarchical restrictions from factors and interaction terms
    #Input
    # - data: data.frame with a column named 'y' containing the response
    # - family: type of GLM, passed onto function 'glm'
    # - maxvars: enumeration is limited to models with <=maxvars
    # - verbose: set to TRUE to display progress information
    # - models: logical matrix with ncol(x) columns selecting variables for each included model
    # Output: object of class msfit, as returned by modelSelection in mombf

    if (missing(priorCoef)) priorCoef= bic()
    if (missing(priorGroup)) priorGroup= bic()
    if (missing(priorVar)) priorVar= igprior(0,0)
    if ((priorCoef@priorDistr != 'bic') | (priorGroup@priorDistr != 'bic')) stop("Only priorCoef=bic() and priorGroup=bic() currently implemented")

    enumerate= TRUE #currently only full enumeration is supported
    family= paste(familyglm$family,familyglm$link)
    familyint= -1 #argument currently ignored, set to value != 0 to prevent inference on the family

    ################# Start of code copied from modelSelection #################

    #Check input
    tmp <- formatInputdata(y=y,x=x,data=data,smoothterms=smoothterms,nknots=nknots,family=family)
    x <- tmp$x; y <- tmp$y; formula <- tmp$formula;
    splineDegree <- tmp$splineDegree
    if (!is.null(tmp$groups)) groups <- tmp$groups
    isgroup <- tmp$isgroup
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
    tmp= codeGroupsAndConstraints(p=p,groups=groups,constraints=constraints)
    ngroups= tmp$ngroups; constraints= tmp$constraints; invconstraints= tmp$invconstraints; nvaringroup=tmp$nvaringroup; groups=tmp$groups

    #Standardize (y,x) to mean 0 and variance 1
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
    #################    End of code copied from modelSelection #################

    #if (missing(phi)) { knownphi <- as.integer(0); phi <- double(0) } else { knownphi <- as.integer(1); phi <- as.double(phi) }
    stdconstants= rbind(c(my,sy),cbind(mx,sx)); colnames(stdconstants)= c('shift','scale')

    #Enumerate models
    if (missing(models)) {
      if (familyint==0) { includeenum= c(includevars[groups+1],FALSE,FALSE) } else { includeenum= includevars[groups+1] }
      models= listmodels(vars2list=1:ngroups, includevars=includeenum, constraints=sapply(constraints,function(z) z+1), nvaringroup=nvaringroup, maxvars=maxvars) #listmodels expects group indexes 1,2,...
      if (familyint==0) models= rbind(cbind(models,FALSE,FALSE),cbind(models,FALSE,TRUE),cbind(models,TRUE,FALSE),cbind(models,TRUE,TRUE))
    } else {
      if (!is.matrix(models) | !is.logical(models)) stop("models must be a matrix of logicals")
      if (dim(models)[2] != ncol(x)) stop("models must have ncol(x) columns")
    }
    nmodels= nrow(models)
    postmean= postvar= matrix(0, nrow=nmodels, ncol=ncol(xstd))
    bicmodel= double(nmodels)
    if (verbose) cat("Enumerating",nmodels,"models")
    iterprogress= round(nmodels/10)
    for (i in 1:nmodels) {
      sel= models[i,]
      if (any(sel)>0) {
        fit1= glm(y ~ -1 + xstd[,sel], data=data, family=familyglm)
      } else { fit1= glm(y ~ -1, data=data, family=familyglm) }
      pm= postmomentsGLM(fit1, priorCoef=priorCoef, priorGroup=priorGroup, isgroup=isgroup)
      postmean[i,sel]= pm$m
      postvar[i,sel]= diag(pm$S)
      bicmodel[i]= BIC(fit1)
      if (verbose && ((i %% iterprogress)==0)) cat('.')
    }
    postProb= -bicmodel/2

    #Add model prior probabilities
    if (priorDelta@priorDistr == "uniform") {
        priorp= apply(models, 1, function(z) { unifPrior(sel=z, logscale=TRUE, groups=groups, constraints=constraints) })
    } else if (priorDelta@priorDistr == "binomial") {
        if ('p' %in% names(priorDelta@priorPars)) {
            priorp= apply(models, 1, function(z) { binomPrior(sel=z, prob=priorDelta@priorPars[['p']], logscale=TRUE, groups=groups, constraints=constraints) })
        } else if (all(c('alpha.p','beta.p') %in% names(priorDelta@priorPars))) {
            priorp= apply(models, 1, function(z) { bbPrior(sel=z, alpha=priorDelta@priorPars['alpha.p'], beta=priorDelta@priorPars['beta.p'], logscale=TRUE, groups=groups, constraints=constraints) })
        } else { stop("priorDelta not recognized") }
    }
    postProb= postProb + priorp
    postMode= models[which.max(postProb),]
    postModeProb= max(postProb)

    #Format output
    postSample= matrix(nrow=0,ncol=p)
    pp= exp(postProb-postModeProb); pp= matrix(pp/sum(pp),ncol=1)
    margpp= as.vector(t(models) %*% pp)
    modelid= apply(models[,1:ncol(xstd),drop=FALSE]==1, 1, function(z) paste(which(z),collapse=','))
    models= data.frame(modelid=modelid,family=family,pp=pp)
    rownames(postmean)= rownames(postvar)= modelid
    colnames(postmean)= colnames(postvar)= colnames(xstd)
    o= order(models$pp,decreasing=TRUE)
    models= models[o,]; postmean= postmean[o,]; postvar= postvar[o,]

    priors= list(priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, priorConstraints=priorDelta, priorVar=priorVar, priorSkew=NULL)
    names(constraints)= paste('group',0:(length(constraints)-1))
    ans= list(postSample=postSample,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb,postmean=postmean,postvar=postvar,family=family,p=ncol(xstd),enumerate=enumerate,priors=priors,ystd=ystd,xstd=xstd,groups=groups,constraints=constraints,stdconstants=stdconstants,outcometype=outcometype,call=call)
    ans$models= models
    new("msfit",ans)
}
