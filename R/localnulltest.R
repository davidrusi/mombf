### Methods for localtest objects

setMethod("show", signature(object='localtest'), function(object) {
  cat("Estimated local covariate effects\n\n")
  print(object$covareffects[1:15,,drop=FALSE])
  if (nrow(object$covareffects)>15) cat("...")
  cat("\n\n")
  cat("Use coef() to extract point estimates, intervals and posterior probabilities for local effects \n")
  #cat("\n\n")
}
)

coef.localtest <- function(object,...) {
    object$covareffects
}

predict.localtest <- function(object, newdata, data, level=0.95, ...) {
    nres= length(object$ms) #number of resolution levels
    if (!missing(newdata)) {
        if (!all(c("x","z") %in% names(newdata))) stop("newdata must be a list with two entries named 'x' and 'z'")
        if (!is.matrix(newdata$x)) stop("newdata$x must be a matrix")
        if (!is.matrix(newdata$z)) newdata$z= as.matrix(newdata$z)
        w= predictDesign(object=object, newdata=newdata)
        ypred= matrix(NA, nrow=nres, ncol=nrow(newdata$x))
        for (i in 1:nres) ypred[i,]= predict(object$ms[[i]], newdata=w[[i]])[,1]
     } else {
        ypred= do.call(rbind,lapply(object$ms, function(z) predict(z)[,1]))
    }
    ypred= colSums(ypred * object$pp_localknots)
    return(ypred)
}


#Create design matrix for newdata given fitted object of type localtest
predictDesign <- function(object, newdata) {
    if (object$Sigma != 'identity') stop("predict method is currently only implemented for independent errors")
    nres= length(object$ms) #number of resolution levels
    ans= vector("list",nres)
    for (j in 1:nres) {
        zdiscrete= zdiscretef= matrix(NA,nrow=nrow(newdata$z),ncol=ncol(newdata$z))
        for (i in 1:ncol(newdata$z)) {
            zdiscretef= cut(newdata$z[,i], breaks=object$regionbounds[[j]][[i]])
            zdiscrete[,i]= as.numeric(zdiscretef)
        }
        region= apply(zdiscrete,1,paste,collapse='.')       
        ans[[j]]= createDesignLocaltest(x=newdata$x, z=newdata$z, region=region, regionbounds=object$regionbounds[[j]], basedegree=object$basedegree, cutdegree=object$cutdegree, knots=object$knots, usecutbasis=object$usecutbasis)$w
    }
    return(ans)
}


setMethod("postProb", signature(object='localtest'), function(object, nmax, method='norm') {
    object$pplocalgrid
}
)


# Local null test for the effect of x on y, at each various regions given by the coordinates in z, using cut spline basis
# Input
# - y: outcome
# - x: matrix with covariate values with nrow(x)==length(y)
# - z: matrix with coordinates with nrow(z)==length(y)
# - x.adjust: optionally, further adjustment covariates to be included in the model with no testing being performed
# - function_id: function identifier. It is assumed that one observes multiple functions over z, this is the identifier of each individual function
# - Sigma: error covariance. By default 'identity', other options are 'MA', 'AR' or 'AR/MA' (meaning that BIC is used to choose between MA and AR). Alternatively the user can supply a function such that Sigma(z[i,],z[j,]) returns the within-function cov(y[i,], y[j,])
# - localgridsize: local test probabilities will be returned for a grid of z values of size localgridsize for each dimension
# - localgrid: regions at which tests will be performed. Defaults to dividing each [min(z[,i]),  max(z[,i])] into 10 equal intervals. If provided, localgrid must be a list with one entry for each z[,i], containing a vector with the desired grid for that z[,i]
# - nbaseknots: number of knots for the spline approximation to the baseline effect of x on y
# - nlocalknots: number of knots for the basis capturing the local effects
# - basedegree: degree of the spline approximation to the baseline
# - cutdegree: degree of the cut spline basis used for testing
# - usecutbasis: if FALSE, then the basis is not cut and a standard spline basis is returned
# - verbose: if set to TRUE some progress information is printed
# - ...: other arguments to be passed on to modelSelection, e.g. family='binomial' for logistic regression
# Output
# - covareffects: estimated local covariate effects at different z's, along with 0.95 posterior intervals and posterior probability for the existence of an effect
# - pplocalgrid: posterior probabilities for the existence of an effect for regions of z values
# - covareffects.mcmc: MCMC output used to build covareffects. Only returned if return.mcmc=TRUE
# - ms: objects of class 'msfit' returned by modelSelection
# - pp_localknots: posterior probability for each resolution level (value of nlocalknots)
# - nlocalknots, basedegree, cutdegree, knots: input parameters
# - regionbounds: list with region bounds defined by the local testing knots at each resolution level
# - Sigma: input parameter

localnulltest= function(y, x, z, x.adjust, localgridsize=100, localgrid, nbaseknots=20, nlocalknots=c(5,10,15), basedegree=3, cutdegree=0, usecutbasis=TRUE, priorCoef=normalidprior(taustd=1), priorGroup=normalidprior(taustd=1), priorDelta=modelbbprior(), mc.cores=min(4,length(nlocalknots)), return.mcmc=FALSE, verbose=FALSE, ...) {
    #Check & format input arguments
    check= checkargs_localnulltest(y=y,x=x,x.adjust=x.adjust,z=z)
    y= check$y; x= check$x; z= check$z; x.adjust= check$x.adjust
    #Define local grid
    if (missing(localgrid)) localgrid= define_localgrid(localgrid=localgrid, localgridsize=localgridsize, z=z)
    #Perform local null tests
    localnulltest_core(y=y, x=x, z=z, x.adjust=x.adjust, localgridsize=localgridsize, localgrid=localgrid, nbaseknots=nbaseknots, nlocalknots=nlocalknots, basedegree=basedegree, cutdegree=cutdegree, usecutbasis=usecutbasis, priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, mc.cores=mc.cores, return.mcmc=return.mcmc, verbose=verbose, ...)
}


#Same as localnulltest, but for data observed on a regular grid (e.g. time series, functional data analysis) where errors may be correlated
localnulltest_fda= function(y, x, z, x.adjust, function_id, Sigma='AR/MA', localgridsize=100, localgrid, nbaseknots=20, nlocalknots=c(5,10,15), basedegree=3, cutdegree=0, usecutbasis=TRUE, priorCoef=momprior(), priorGroup=groupmomprior(), priorDelta=modelbbprior(), mc.cores=min(4,length(nlocalknots)), return.mcmc=FALSE, verbose=FALSE, ...) {
    #Check & format input requirements
    check= checkargs_localnulltest(y=y,x=x,z=z,x.adjust=x.adjust,function_id=function_id)
    y= check$y; x= check$x; z= check$z; x.adjust= check$x.adjust
    #Define local grid
    if (missing(localgrid)) localgrid= define_localgrid(localgrid=localgrid, localgridsize=localgridsize, z=z) #define localgrid
    #Perform local null tests
    localnulltest_core(y=y, x=x, z=z, x.adjust=x.adjust, function_id=function_id, Sigma=Sigma, localgridsize=localgridsize, localgrid=localgrid, nbaseknots=nbaseknots, nlocalknots=nlocalknots, basedegree=basedegree, cutdegree=cutdegree, usecutbasis=usecutbasis, priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, mc.cores=mc.cores, return.mcmc=return.mcmc, verbose=verbose, ...)
}

#Check input arguments to localnulltest
checkargs_localnulltest= function(y, x, z, x.adjust, function_id) {
    if (is.matrix(y)) y= as.vector(y)
    if (!is.matrix(x)) x= as.matrix(x)
    if (!is.matrix(z)) z= as.matrix(z)
    n= length(y)
    if (nrow(x) != n) stop("nrow(x) must be equal to length(y)")
    if (nrow(z) != n) stop("nrow(z) must be equal to length(y)")
    if (!missing(function_id)) if (length(function_id) != n) stop("length(function_id) must be equal to length(y)")
    if (!is.numeric(x)) stop("x must be numeric")
    if (!is.numeric(z)) stop("z must be numeric")
    if (missing(x.adjust)) {
        x.adjust= NULL
    } else {
        if (!is.null(x.adjust)) {
            x.adjust= as.matrix(x.adjust)
            x.adjust= x.adjust[, apply(x.adjust, 2, sd) > 0, drop=FALSE]  #drop the intercept (and constant columns)
            if (!is.numeric(x.adjust)) stop("x.adjust must be numeric")
        }
    }
    return(list(y=y, x=x, z=z, x.adjust=x.adjust))
}

#Core function to run local null tests at multiple resolutions.
# It calls localnulltest_givenknots (for iid errors) or localnulltest_fda_givenknots (for dependent errors observed on regular grids)
localnulltest_core= function(y, x, z, x.adjust, function_id, Sigma, localgridsize=100, localgrid, nbaseknots=20, nlocalknots=c(5,10,15), basedegree=3, cutdegree=0, usecutbasis=TRUE, priorCoef=normalidprior(taustd=1), priorGroup=normalidprior(taustd=1), priorDelta=modelbbprior(), mc.cores=min(4,length(nlocalknots)), return.mcmc=FALSE, verbose=FALSE, ...) {
    #Define function to be run at each resolution level
    if (missing(Sigma)) {
        foo= function(i) { localnulltest_givenknots(y=y, x=x, z=z, x.adjust=x.adjust, localgrid=localgrid, nbaseknots=nbaseknots, nlocalknots=nlocalknots[i], basedegree=basedegree, cutdegree=cutdegree, usecutbasis=usecutbasis, priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, verbose=verbose, ...) }
        Sigma= "identity"
    } else {
        foo= function(i) { localnulltest_fda_givenknots(y=y, x=x, z=z, x.adjust=x.adjust, function_id=function_id, Sigma=Sigma, localgrid=localgrid, nbaseknots=nbaseknots, nlocalknots=nlocalknots[i], basedegree=basedegree, cutdegree=cutdegree, usecutbasis=usecutbasis, priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, verbose=verbose, ...) }
    }
    #Run analysis for each resolution level
    if (("parallel" %in% loadedNamespaces()) & mc.cores>1)  {
      alltests= parallel::mclapply(1:length(nlocalknots), foo, mc.cores=mc.cores)
    } else {
      alltests= lapply(1:length(nlocalknots), foo)
    }
    #Posterior probability of each resolution level (value of nlocalknots)
    logpp= vector("list", length(nlocalknots))
    maxlogpp= double(length(nlocalknots))
    for (i in 1:length(logpp)) {
        logpp[[i]]= logjoint(alltests[[i]]$ms[[1]], return_models=FALSE) - 0.5 * alltests[[i]]$logdetSinv
        maxlogpp[i]= max(logpp[[i]])
    }
    m= max(maxlogpp)
    pp_localknots= sapply(logpp, function(zz) sum(exp(zz-m)))
    pp_localknots= pp_localknots / sum(pp_localknots)
    #Average across resolution levels       
    if(length(nlocalknots)==1) {
        ans= alltests[[1]]
    } else {
        #Posterior probabilities for the local tests
        pplocalgrid= alltests[[1]]$pplocalgrid
        pplocalgrid[,-1]= pplocalgrid[,-1] * pp_localknots[1]
        for (i in 2:length(nlocalknots)) pplocalgrid[,-1]= pplocalgrid[,-1] + alltests[[i]]$pplocalgrid[,-1] * pp_localknots[i]
        #BMA estimate and 95% intervals
        covareffects= alltests[[1]]$covareffects
        p= ncol(alltests[[1]]$covareffects.mcmc)
        B= sapply(alltests, function(zz) nrow(zz$covareffects.mcmc))
        w= rep(pp_localknots, B)
        for (i in 1:p) {
            samples= sapply(alltests, function(zz) zz$covareffects.mcmc[,i])
            covareffects[i,'estimate']= weighted.mean(samples, w)
            covareffects[i,4:5]= quantile_weighted(samples, weights=w, probs = c(0.025,0.975))
        }
        if (return.mcmc) {
            covareffects.mcmc= cbind(do.call(rbind, lapply(alltests, "[[", "covareffects.mcmc")), weight=w/sum(w))
        } else {
            covareffects.mcmc= NULL
        }
        ms= lapply(alltests, function(zz) zz[["ms"]][[1]])
        regionbounds= lapply(alltests, function(zz) zz[["regionbounds"]][[1]])
        knots= lapply(alltests, function(zz) zz[["knots"]][[1]])
        ans= list(pplocalgrid=pplocalgrid, covareffects=covareffects, covareffects.mcmc=covareffects.mcmc, ms=ms, pp_localknots=pp_localknots, nlocalknots=nlocalknots, regionbounds=regionbounds, basedegree=basedegree, cutdegree=cutdegree, usecutbasis=usecutbasis, knots=knots, Sigma=Sigma)
    }
    new("localtest",ans)
}




localnulltest_givenknots= function(y, x, z, x.adjust, localgridsize=100, localgrid, nbaseknots=20, nlocalknots=10, basedegree=3, cutdegree=0, usecutbasis=TRUE, priorCoef=normalidprior(taustd=1), priorGroup=normalidprior(taustd=1), priorDelta=modelbbprior(), verbose=FALSE, ...) {
    #Check & format input requirements
    check= checkargs_localnulltest(y=y, x=x, z=z, x.adjust=x.adjust)
    y= check$y; x= check$x; z= check$z; x.adjust= check$x.adjust
    # Define local tests
    if (missing(localgrid)) localgrid= define_localgrid(localgrid=localgrid, localgridsize=localgridsize, z=z)
    # Define knots and corresponding knot-based testing regions
    kk= define_knots_localnulltest(z=z, localgrid=localgrid, nbaseknots=nbaseknots, basedegree=basedegree, nlocalknots=nlocalknots)
    knots= kk$knots; regionbounds= kk$regionbounds; region=kk$region; regioncoord= kk$regioncoord; testov= kk$testov; testxregion= kk$testxregion; testIntervals= kk$testIntervals
    #Create design matrix & run Bayesian model selection
    desnew= estimationPoints(x=x, regioncoord=regioncoord, regionbounds=regionbounds, testov=testov) #Points at which local effects will be estimated
    des= createDesignLocaltest(x=rbind(x,desnew$x), z=rbind(z,desnew$z), y=y, region=c(region,desnew$region), regionbounds=regionbounds, basedegree=basedegree, cutdegree=cutdegree, knots=knots, usecutbasis=usecutbasis, useSigma=FALSE)
    w= des$w[1:nrow(x),]; wnew= des$w[-1:-nrow(x),]
    if (!is.null(x.adjust)) {
        w= cbind(w, x.adjust)
        groups= c(des$vargroupsn, rep(max(des$vargroupsn)+1, ncol(x.adjust)))
    } else {
        groups= des$vargroupsn
    }
    ms= modelSelection(y, x=w, groups=groups, priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, verbose=verbose, ...)        
    pp= postProb(ms)
    #Obtain posterior probability for each region defined by the knots. Store in regionid, ppregion
    modelid= modelid2logical(pp$modelid, nvars=ncol(ms$xstd))
    regionnames= unique(des$vargroups[-1:-des$ncolw0])
    firstvaringroup= match(regionnames, des$vargroups)
    regionid= modelid[,firstvaringroup,drop=FALSE]
    colnames(regionid)= des$vargroups[firstvaringroup]
    regionidtext= apply(regionid, 1, function(z) paste(which(z), collapse=','))
    ppregion= aggregate(pp$pp ~ regionidtext, FUN=sum)
    names(ppregion)= c('regionid','pp')
    regionid= modelid2logical(ppregion$regionid, nvars=ncol(regionid))
    colnames(regionid)= des$vargroups[firstvaringroup]
    #Find posterior probability for each local test
    varids= unique(des$w1varname)
    pplocaltest= matrix(NA, nrow=length(testxregion), ncol=ncol(x))
    colnames(pplocaltest)= varids
    for (i in 1:length(varids)) {
        sel= grep(paste("^",varids[i],sep=""), colnames(regionid)) #inclusion indicators for variable varids[i]
        nn= colnames(regionid)[sel]
        nncut= sub("\\.","", sub(varids[i],'',nn))
        colnames(regionid)[sel]= nncut
        pplocaltest[,i]= postProblocaltest(regionid[,sel,drop=FALSE], ppregion$pp, testxregion)
        colnames(regionid)[sel]= nn
    }
    #Estimated effect at the center of each local test region
    if (!is.null(x.adjust)) { m.adjust= colMeans(x.adjust) } else { m.adjust= NULL }
    est= estimateLocaleffect(ms, wnew, desnew, x.adjust=m.adjust)
    #Return output
    pplocalgrid= data.frame(localtest=1:nrow(testIntervals), testIntervals, pplocaltest)
    ans= list(pplocalgrid=pplocalgrid, covareffects=est$covareffects, covareffects.mcmc=est$covareffects.mcmc, ms=list(ms), pp_localknots=1, logdetSinv=0, nlocalknots=nlocalknots, regionbounds=list(regionbounds), basedegree=basedegree, cutdegree=cutdegree, usecutbasis=usecutbasis, knots=knots, Sigma='identity')
    new("localtest",ans)
}



localnulltest_fda_givenknots= function(y, x, z, x.adjust, function_id, Sigma='AR/MA', localgridsize=100, localgrid, nbaseknots=20, nlocalknots=10, basedegree=3, cutdegree=0, usecutbasis=TRUE, priorCoef=normalidprior(taustd=1), priorGroup=normalidprior(taustd=1), priorDelta=modelbbprior(), verbose=FALSE, ...) {
    #Check & format input requirements
    check= checkargs_localnulltest(y=y, x=x, z=z, x.adjust=x.adjust, function_id=function_id)
    y= check$y; x= check$x; z= check$z; x.adjust= check$x.adjust
    if (inherits(Sigma, "character")) {
        if (!(Sigma %in% c('MA','AR','AR/MA'))) stop("This Sigma is not implemented. Use 'MA', 'AR', 'AR/MA', or a user-supplied function")
    } else {
        if (!inherits(Sigma, "function")) stop("Sigma must either be character or a function such that Sigma(i,j)=cov(epsilon_i,epsilon_j)")
    }
    # Sort data according to function_id and z
    o= do.call(order, data.frame(function_id, z))
    y= y[o]; x= x[o,,drop=FALSE]; z= z[o,,drop=FALSE]; function_id= function_id[o]
    if (!is.null(x.adjust)) x.adjust= x.adjust[o,,drop=FALSE]
    # Define local tests
    if (missing(localgrid)) localgrid= define_localgrid(localgrid=localgrid, localgridsize=localgridsize, z=z)
    # Define knots and corresponding knot-based testing regions
    kk= define_knots_localnulltest(z=z, localgrid=localgrid, nbaseknots=nbaseknots, basedegree=basedegree, nlocalknots=nlocalknots)
    knots= kk$knots; regionbounds= kk$regionbounds; region= kk$region; regioncoord= kk$regioncoord; testov= kk$testov; testxregion= kk$testxregion; testIntervals= kk$testIntervals
    #Create design matrix & run Bayesian model selection
    xc= scale(x, center=TRUE, scale=FALSE)
    desnew= estimationPoints(x=xc, regioncoord=regioncoord, regionbounds=regionbounds, testov=testov) #Points at which local effects will be estimated
    des= createDesignLocaltest(x=xc, z=z, y=y, x.adjust=x.adjust, xnew=desnew$x, znew=desnew$z, function_id=function_id, region=region, regionnew=desnew$region, regionbounds=regionbounds, basedegree=basedegree, cutdegree=cutdegree, knots=knots, usecutbasis=usecutbasis, useSigma=TRUE, Sigma=Sigma)
    ytilde= des$ytilde; wtilde= des$wtilde; logdetSinv= des$logdetSinv
    wnew= des$wnew
    if (!is.null(x.adjust)) {
        wtilde= cbind(wtilde, des$wtilde.adjust)
        groups= c(des$vargroupsn, rep(max(des$vargroupsn)+1, ncol(x.adjust)))
    } else {
        groups= des$vargroupsn
    }
    ms= modelSelection(y=ytilde, x=wtilde, groups=groups, priorCoef=priorCoef, priorGroup=priorGroup, priorDelta=priorDelta, verbose=verbose, ...)
    pp= postProb(ms)
    #Obtain posterior probability for each region defined by the knots. Store in regionid, ppregion
    modelid= modelid2logical(pp$modelid, nvars=ncol(ms$xstd))
    regionnames= unique(des$vargroups[-1:-des$ncolw0])
    firstvaringroup= match(regionnames, des$vargroups)
    regionid= modelid[,firstvaringroup,drop=FALSE]
    colnames(regionid)= des$vargroups[firstvaringroup]
    regionidtext= apply(regionid, 1, function(z) paste(which(z), collapse=','))
    ppregion= aggregate(pp$pp ~ regionidtext, FUN=sum)
    names(ppregion)= c('regionid','pp')
    regionid= modelid2logical(ppregion$regionid, nvars=ncol(regionid))
    colnames(regionid)= des$vargroups[firstvaringroup]
    #Find posterior probability for each local test
    varids= unique(des$w1varname)
    pplocaltest= matrix(NA, nrow=length(testxregion), ncol=ncol(x))
    colnames(pplocaltest)= varids
    for (i in 1:length(varids)) {
        sel= grep(paste("^",varids[i],sep=""), colnames(regionid)) #inclusion indicators for variable varids[i]
        nn= colnames(regionid)[sel]
        nncut= sub("\\.","", sub(varids[i],'',nn))
        colnames(regionid)[sel]= nncut
        pplocaltest[,i]= postProblocaltest(regionid[,sel,drop=FALSE], ppregion$pp, testxregion)
        colnames(regionid)[sel]= nn
    }
    #Estimated effect at the center of each local test region
    if (!is.null(x.adjust)) { m.adjust= colMeans(x.adjust) } else { m.adjust= NULL }
    est= estimateLocaleffect(ms, wnew, desnew, x.adjust=m.adjust)
    #Return output
    pplocalgrid= data.frame(localtest=1:nrow(testIntervals), testIntervals, pplocaltest)
    ans= list(pplocalgrid=pplocalgrid, covareffects=est$covareffects, covareffects.mcmc=est$covareffects.mcmc, ms=list(ms), pp_localknots=1, logdetSinv=logdetSinv, nlocalknots=nlocalknots, regionbounds=list(regionbounds), basedegree=basedegree, cutdegree=cutdegree, usecutbasis=usecutbasis, knots=knots)
    new("localtest",ans)
}


# Define knots and corresponding knot-based testing regions
# Input
# - z: matrix with n rows and d colums indicating the d-dimensional coordinates for each observation
# - localgrid: regions at which tests will be performed. Defaults to dividing each [min(z[,i]),  max(z[,i])] into 10 equal intervals. If provided, localgrid must be a list with one entry for each z[,i], containing a vector with the desired grid for that z[,i]
# - nbaseknots: number of knots for the baseline basis
# - basedegree: degree of the baseline basis
# - nlocalknots: number of knots for the local testing basis
# Output
# - knots: a list with d elements containing the knots for each dimension
# - regionbounds: a list with d elemants containing the local testing region bounds
# - region: character vector indicating the region of each observation (row in z)
# - regioncoord: region coordinates
# - testov: overlaps between localgrid and regioncoord
# - testxregion: list with an entry for each test, indicating the regions overlapping that test
# - testIntervals: local test grid, formatted as intervals
define_knots_localnulltest= function(z, localgrid, nbaseknots, basedegree, nlocalknots) {
    knots= vector("list", ncol(z))
    regionbounds= vector("list",ncol(z))
    for (i in 1:ncol(z)) {
        zlim= range(z[,i])
        baseknotwidth= (zlim[2] - zlim[1]) / nbaseknots
        knots[[i]]= seq(zlim[1] - baseknotwidth*(0.5+basedegree), zlim[2] + baseknotwidth*(0.5+basedegree), by=baseknotwidth)
        if (zlim[2] - zlim[1] > 0.01) knots[[i]]= round(knots[[i]], 3)
        for (i in 1:ncol(z)) {
            zseq= seq(zlim[1], zlim[2], length=nlocalknots)
            regionbounds[[i]]= c(zseq[1] - (zseq[2]-zseq[1]), zseq, zseq[length(zseq)] + (zseq[2]-zseq[1]))
            if (zlim[2] - zlim[1] > 0.01) regionbounds[[i]]= round(regionbounds[[i]], 3)
        }
    }
    zdiscrete= zdiscretef= matrix(NA,nrow=nrow(z),ncol=ncol(z))
    regionlabels= vector("list", length(ncol(z)))
    for (i in 1:ncol(z)) {
        zdiscretef= cut(z[,i], breaks=regionbounds[[i]])
        regionlabels[[i]]= data.frame(regionid=1:length(levels(zdiscretef)), start=regionbounds[[i]][-length(regionbounds[[i]])], end=regionbounds[[i]][-1])
        zdiscrete[,i]= as.numeric(zdiscretef)
    }
    region= apply(zdiscrete,1,paste,collapse='.')
    #Store coordinates for each region
    regioncoord= regionCoordinates(unique(zdiscrete), regionlabels)
    #Find regions overlapping each local test
    testov= testOverlaps(localgrid, regioncoord)
    testxregion= testov$testxregion                           #regions overlapping each local test
    testIntervals= do.call(data.frame, testov$testIntervals)  #local tests formatted as intervals
    label= rep(paste('z',1:ncol(z),sep=''), each=2)
    label= paste(label, rep(c('.low','.high'),ncol(z)), sep='')
    names(testIntervals)= label
    #Return output
    ans= list(knots=knots, regionbounds=regionbounds, region=region, regioncoord=regioncoord, testov=testov, testxregion=testxregion, testIntervals=testIntervals)
    return(ans)
}


#Estimate error correlation matrix
# Input
# - y: outcome
# - w: basis to be fed into modelSelection
# - z: matrix with coordinates with nrow(z)==length(y)
# - function_id: function identifier. Independence is assumed across functions
# - region: factor indicating what region (based on a partition of z) each observation belongs to
# - method: method to estimate the covariance. 'MA', 'AR', or 'AR/MA' (meaning trying both AR and MA and choose the one with best BIC)
# Output: a function Sigma such that Sigma(z[i,], z[j,])= cov(z[i,], z[j,])
#   For method=='MA', based on estimating the MA1 parameter y_t= theta * epsilon_{t-1} + epsilon_t
#   For method=='AR', based on estimating the AR1 parameter y_t= theta * y_{t-1} + epsilon_t

estimateSigma_localtest= function(y, w, z, function_id, region, method='AR/MA') {
    e= residuals(lm(y ~ w))
    if (ncol(z)>1) stop("If z is univariate, Sigma must be provided")
    if (!(method %in% c('MA','AR','AR/MA'))) stop("Method to estimate Sigma must be 'MA', 'AR', or 'AR/MA'")
    #
    ids= unique(function_id)
    theta= bics= matrix(NA, nrow=length(ids), ncol=2)
    colnames(theta)= colnames(bics)= c('MA','AR')
    for (i in 1:length(ids)) {
        sel= (function_id == ids[i])
        if (method %in% c('MA','AR/MA')) {
            fit.ma= try(arima(e[sel], order=c(0,0,1), method='ML'), silent=TRUE)
            if (inherits(fit.ma, "try-error")) fit.ma= try(arima(e[sel], order=c(0,0,1)), silent=TRUE)
            if (!inherits(fit.ma, "try-error")) {
                theta[i,'MA']= coef(fit.ma)['ma1']
                bics[i,'MA']= BIC(fit.ma)
            }
        }
        if (method %in% c('AR','AR/MA')) {
            fit.ar= try(arima(e[sel], order=c(1,0,0), method='ML'), silent=TRUE)
            if (inherits(fit.ar, "try-error")) fit.ar= try(arima(e[sel], order=c(1,0,0)), silent=TRUE)
            if (!inherits(fit.ar, "try-error")) {
                theta[i,'AR']= coef(fit.ar)['ar1']
                bics[i,'AR']= BIC(fit.ar)
            }
        }
    }
    if (method=='AR/MA') method= ifelse(which.min(colMeans(bics, na.rm=TRUE)) == 1, 'MA', 'AR')
    #Store the correlation matrix
    zregion= unique(cbind(as.numeric(region), z))
    if (method=='MA') {
        theta= mean(theta[,'MA'], na.rm=TRUE)
        ans= evalSigma_MA1_localtest(theta=theta, region=zregion[,1])   
        #fun= function(z1, z2) { if (z1==z2) return(1) else if (abs(z1-z2)==zdif) return(theta/(1+theta^2)) else return(0) }
    } else {
        theta= mean(theta[,'AR'], na.rm=TRUE)
        ans= evalSigma_AR1_localtest(theta=theta, region=zregion[,1])   
    }
    return(ans)
}



#Compute correlation matrix for MA1 model. The covariance is cov(i,i)= 1+theta^2, cov(i,i+1)= theta, cov(i,j)=0 otherwise
# Input
# - theta: MA1 parameter, i.e. y[t]= epsilon[t-1] * theta + epsilon[t]
# - region: factor indicating what region each observation belongs to
# - function_id: function identifier. Independence is assumed across functions
# Output: a list with one element per region, containing Sigma= corr(e) for that region. The list names match the unique values in region
evalSigma_MA1_localtest= function(theta, region) {
    regionids= unique(region)
    ans= vector("list", length(regionids))
    names(ans)= regionids
    for (i in 1:length(regionids)) {
        sel= which(region == regionids[i])
        ans[[i]]= diag(length(sel))
        ans[[i]][abs(col(ans[[i]]) - row(ans[[i]])) == 1]= theta / (1+theta^2)
        rownames(ans[[i]])= colnames(ans[[i]])= sel #set row/column names, for later use
    }
    return(ans)
}


#Compute correlation matrix for AR1 model, i.e. is cor(i,j)= theta^{|i-j|}
# Input
# - theta: AR1 parameter, i.e. y[t]= epsilon[t-1] * theta + epsilon[t]
# - region: factor indicating what region each observation belongs to
# - function_id: function identifier. Independence is assumed across functions
# Output: a list with one element per region, containing Sigma= corr(e) for that region. The list names match the unique values in region
evalSigma_AR1_localtest= function(theta, region) {
    regionids= unique(region)
    ans= vector("list", length(regionids))
    names(ans)= regionids
    for (i in 1:length(regionids)) {
        sel= which(region == regionids[i])
        ans[[i]]= matrix(NA, nrow=length(sel), ncol=length(sel))
        ans[[i]]= theta^abs(col(ans[[i]]) - row(ans[[i]]))
        rownames(ans[[i]])= colnames(ans[[i]])= sel #set row/column names, for later use
    }
    return(ans)
}


#Compute covariance matrix from a user-given function
# Input
# - Sigma: function such that Sigma(i,j)= cov(y[i],y[j]), which is assumed to only depend on z[i,] and z[j,]
# - region: factor indicating what region (based on a partition of z) each observation belongs to
# - z: coordinates
# Output: a list with one element per region, containing Sigma= cov(e) for that region. The list names match the unique values in region
evalSigma_localtest= function(Sigma, region, z) {
    regionids= unique(region)
    ans= vector("list", length(regionids))
    names(ans)= regionids
    for (i in 1:length(regionids)) {
        sel= which(region == regionids[i])
        ans[[i]]= matrix(NA, nrow=length(sel), ncol=length(sel))
        ans[[i]][1,1]= Sigma(z1=z[sel[1],], z2=z[sel[1],])
        if (length(sel)>1) {
            for (j in 2:length(sel)) {
                for (k in 1:j) {
                    ans[[i]][j,k]= ans[[i]][k,j]= Sigma(z1=z[sel[j],], z2=z[sel[k],])
                }
            }
        }
        rownames(ans[[i]])= colnames(ans[[i]])= sel #set row/column names, for later use
    }
    return(ans)
}



#Compute uncorrelated outcome and design matrix given inverse of error covariance, considering only within-region terms
# Input
# - y: outcome
# - w: design matrix
# - x.adjust: optional argument. Adjustment covariates (no testing performed on these)
# - region: factor indicating the regions to which each observation in y belongs to
# - S: error covariance. A list with an entry for each region with Sigma for that region
# - function_id: function identifier. Independence across functions is assumed
# Output: a list with elements
# - ytilde: block-diag(S)^{-1/2} y, where the blocks are defined by region
# - wtilde: block-diag(S)^{-1/2} w, where the blocks are defined by region
# - wtilde.adjust: block-diag(S}^{-1/2} x.adjust
# - logdetSinv: log-determinant of Sinv
uncorrelate_outcome= function(y, w, x.adjust, region, S, function_id) {
    n= length(y)
    ytilde= double(n)
    wtilde= matrix(NA, nrow=nrow(w), ncol=ncol(w))
    if (!is.null(x.adjust)) {
        wtilde.adjust= matrix(NA, nrow=nrow(x.adjust), ncol=ncol(x.adjust))
    } else {
        wtilde.adjust= NULL
    }
    logdetSinv= 0
    regionids= unique(region)
    functionids= unique(function_id)
    nfunctions= length(functionids)
    for (i in 1:length(regionids)) {
        e= eigen(S[[regionids[i]]])
        sqSinv= e$vectors %*% diag(1/sqrt(e$values), nrow=nrow(S[[regionids[i]]])) %*% t(e$vectors)
        logdetSinv= logdetSinv - nfunctions * sum(log(e$values))
        seli= (region == regionids[i])
        for (j in functionids) {
            sel= seli & (function_id == j)
            ytilde[sel]= sqSinv %*% matrix(y[sel],ncol=1)
            wtilde[sel,]= sqSinv %*% w[sel,,drop=FALSE]
            if (!is.null(x.adjust)) wtilde.adjust[sel,]= sqSinv %*% x.adjust[sel,,drop=FALSE]
        }
    }
    ans= list(ytilde=ytilde, wtilde=wtilde, wtilde.adjust=wtilde.adjust, logdetSinv=logdetSinv)
    return(ans)
}



#Compute quantiles from weighted sample. Adapted from function Quantile in package DescTools
quantile_weighted= function (x, weights = NULL, probs = c(0.025,0.975), na.rm = FALSE, names = TRUE, type = 7, digits = 7) {
    format_perc <- function(x, digits = max(2L, getOption("digits")), 
        probability = TRUE, use.fC = length(x) < 100, ...) {
        if (length(x)) {
            if (probability) x <- 100 * x
            ans <- paste0(if (use.fC) 
                formatC(x, format = "fg", width = 1, digits = digits)
            else format(x, trim = TRUE, digits = digits, ...), 
                "%")
            ans[is.na(x)] <- ""
            ans
        }
        else character(0)
    }
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    n <- length(x)
    if (n == 0 || (!isTRUE(na.rm) && any(is.na(x)))) {
        return(rep.int(NA, length(probs)))
    }
    if (!is.null(weights)) {
        if (!is.numeric(weights)) stop("'weights' must be a numeric vector")
        else if (length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        }
        else if (!all(is.finite(weights))) 
            stop("missing or infinite weights")
        if (any(weights < 0)) stop("negative weights")
        if (!is.numeric(probs) || all(is.na(probs)) || isTRUE(any(probs < 0 | probs > 1))) {
            stop("'probs' must be a numeric vector with values in [0,1]")
        }
        if (all(weights == 0)) {
            warning("all weights equal to zero")
            return(rep.int(0, length(probs)))
        }
    }
    if (isTRUE(na.rm)) {
        indices <- !is.na(x)
        x <- x[indices]
        n <- length(x)
        if (!is.null(weights)) weights <- weights[indices]
    }
    order <- order(x)
    x <- x[order]
    weights <- weights[order]
    if (is.null(weights)) rw <- (1:n)/n
    else rw <- cumsum(weights)/sum(weights)
    if (type == 5) {
        qs <- sapply(probs, function(p) {
            if (p == 0) 
                return(x[1])
            else if (p == 1) 
                return(x[n])
            select <- min(which(rw >= p))
            if (rw[select] == p) 
                mean(x[select:(select + 1)])
            else x[select]
        })
    }
    else if (type == 7) {
        if (is.null(weights)) {
            index <- 1 + max(n - 1, 0) * probs
            lo <- pmax(floor(index), 1)
            hi <- ceiling(index)
            x <- sort(x, partial = if (n == 0) numeric() else unique(c(lo, hi)))
            qs <- x[lo]
            i <- which((index > lo & x[hi] != qs))
            h <- (index - lo)[i]
            qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
        }
        else {
            n <- sum(weights)
            ord <- 1 + (n - 1) * probs
            low <- pmax(floor(ord), 1)
            high <- pmin(low + 1, n)
            ord <- ord%%1
            allq <- approx(cumsum(weights), x, xout = c(low, high), method = "constant", f = 1, rule = 2, ties="mean")$y
            k <- length(probs)
            qs <- (1 - ord) * allq[1:k] + ord * allq[-(1:k)]
        }
    }
    else {
        qs <- NA
        warning(gettextf("type %s is not implemented", type))
    }
    if (names && length(probs) > 0L) {
        stopifnot(is.numeric(digits), digits >= 1)
        names(qs) <- format_perc(probs, digits = digits)
    }
    return(qs)
}


#Define the local testing regions. In ncol(z)>=2 dimensions these are rectangular-shaped
# Input
# - localgrid: if not missing, this is returned as the function's output after checking it has the right format
# - localgridsize: number of local tests that one wishes to perform. The number of grid points in each dimension is localgridsize^(1/dim)
# - z: coordinates
# Output: list of length ndim defining a grid for each coordinate
define_localgrid= function(localgrid, localgridsize, z) {
    ndim= ncol(z)
    if (!missing(localgrid)) {
        if (!list(localgrid)) stop(paste("localgrid should be a list of length",ndim,"defining a grid for each coordinate in z (e.g. just 1 if z is univariate"))
    } else {
        if (localgridsize <= 1) stop("You defined too few localgridsize")
        splitrange= function(zz) {
            zlim= range(zz)
            ans= seq(zlim[1], zlim[2], length=localgridsize)
            if (zlim[2] - zlim[1] > 0.01) ans= round(ans,3)
            return(ans)
        }
        localgrid= apply(z, 2, splitrange, simplify=FALSE)
        #localgrid= apply(z, 2, function(zz) seq(min(zz), max(zz), length=localgridsize), simplify=FALSE)
    }
    names(localgrid)= paste('dim', 1:length(localgrid), sep='')
    return(localgrid)
}


#Estimate via MCMC the local effects at wnew
# Input
# - ms: output from modelSelection
# - wnew: design points at which to estimate the local effects
# - desnew: output of estimationPoints containing further info about wnew
# - x.adjust: optional argument. The mean value of adjustment covariates included in ms, for which no testing is performed
# Output
# - covareffects: data.frame indicating the effect of each covariate at each region (posterior mean and 95% intervals)
# - covareffects.mcmc: mcmc draws for the effects summarized in covareffects
estimateLocaleffect= function(ms, wnew, desnew, x.adjust) {
    th= rnlp(msfit=ms, niter=10^4)
    newdata= wnew[desnew$covareffect$value==1,] - wnew[desnew$covareffect$value==0]
    if (!is.null(x.adjust)) {
        newdata= cbind(newdata, matrix(rep(x.adjust, nrow(newdata)), nrow=nrow(newdata), byrow=TRUE))
    }
    sel= (desnew$covareffect$value==1)
    covareffects= data.frame(covariate=desnew$covareffect$covariate[sel], desnew$z[sel,,drop=FALSE])
    thsel= !(colnames(th) %in% c('intercept','phi'))
    covareffects.mcmc= th[,thsel] %*% t(newdata)
    covareffects$estimate= colMeans(covareffects.mcmc)
    qq= t(apply(covareffects.mcmc, 2, quantile, probs=c(.025,.975)))
    colnames(qq)= c("2.5%","97.5%")
    margpp= colMeans(covareffects.mcmc != 0)
    covareffects= cbind(covareffects, qq, margpp=margpp)
    ans= list(covareffects= covareffects, covareffects.mcmc=covareffects.mcmc)
    return(ans)
}


#Obtain points at which local effects will be estimated
# - If z is not missing, local effects are estimated at these z coordinates
# - If z is missing, these are the center of each local test's region with each x set to 0 or 1
# Input
# - x: x values
# - z: z values at which to estimate the local effects
# - regioncoord: region coordinates
# - regionbounds: boundaries of each region
# - testov: overlap between local tests and regions, as returned by testOverlaps
# Output
# - x: x values at which to estimate the effect. These will be (0,1) for each test region
# - z: z values at which to estimate the effect. If z is missing these are the local test region centers, else the input z's
# - region: region corresponding to each row in z
# - covareffect: data.frame indicating the covariate corresponding to each row in x, and its value (0 or 1)
estimationPoints= function(x, z, regioncoord, regionbounds, testov) {
    if (missing(z)) {
        #Obtain center point of each test interval
        z= lapply(testov$testIntervals, rowMeans)
        z= as.matrix(expand.grid(z))
        z= rbind(z, z) #duplicate matrix for x at mean(x) and at mean(x)+1
    }
    colnames(z)= paste('z',1:ncol(z),sep='')
    #Figure out region for each row in z
    zregion= matrix(NA, nrow=nrow(z), ncol=ncol(z))
    for (i in 1:ncol(z)) {
        regionid= intervals::interval_overlap(z[,i], intervals::Intervals(regioncoord[[i]]))
        regionid= sapply(regionid, '[[', 1) #if there are multiple hits, then return the 1st one
        regionnames= rownames(regioncoord[[i]])[regionid]
        zregion[,i]= as.numeric(sub("R","",regionnames))
    }
    zregion= apply(zregion,1,paste,collapse='.')
    #Build design matrix
    m= colMeans(x)
    xmid= matrix(rep(m,each=nrow(z)), nrow=nrow(z))
    xdes= zdes= vector("list",ncol(x))
    for (j in 1:ncol(x)) {
        xmid[,j]= rep(0:1, each=nrow(xmid)/2)
        xdes[[j]]= xmid
        zdes[[j]]= z
        xmid[,j]= m[j]
    }
    covariate= rep(1:ncol(x), each=nrow(z))
    value= rep(rep(0:1,each=nrow(z)/2), ncol(x))
    ans= list(x=do.call(rbind,xdes), z=do.call(rbind,zdes), region=rep(zregion,ncol(x)), covareffect=data.frame(covariate,value))
    return(ans)
}


#Return coordinates of each region
# Input:
# - zdiscrete: discretized z coordinates, with each dimension in a separate column
# - regionlabels: a list with an entry for each dimension, indicating the start and end of the interval for each discretized z coordinate
# Output: a list with an entry per dimension, indicating the coordinates of each region in that dimension
regionCoordinates= function(zdiscrete, regionlabels) {
    #zdiscrete= unique(zdiscrete)
    nn= paste('R', apply(zdiscrete,1,paste,collapse='.'), sep='')
    regioncoord= lapply(1:ncol(zdiscrete), function(i) { ans= matrix(NA, nrow=nrow(zdiscrete), ncol=2); rownames(ans)= nn; return(ans) })
    names(regioncoord)= paste('dim',1:length(regioncoord),sep='')
    for (j in 1:ncol(zdiscrete)) {
        for (i in 1:nrow(regioncoord[[j]])) {
            sel= zdiscrete[i,j]
            regioncoord[[j]][i,]= as.double(regionlabels[[j]][sel,c('start','end')])
        }
        regioncoord[[j]]= regioncoord[[j]][order(nn),]
    }
    return(regioncoord)
}


#Create design matrix with baseline effects of x and orthogonalized cut basis effects of x interacted with z
# Input
# - x: covariate values
# - z: coordinates
# - y: outcome (only used if useSigma is TRUE)
# - x.adjust: optional argument. Adjustment covariate values (no testing performed on these)
# - xnew: covariate values for new observations
# - znew: covariate values for new observations
# - function_id: function identifier (only used if useSigma is TRUE)
# - region: identifier of the region that each row in z belongs to
# - regionnew: identifier of the region that each row in znew belongs to
# - regionbounds: boundaries of the regions
# - basedegree: baseline spline basis degree
# - cutdegree: cut spline basis degree
# - knots: knots for baseline spline
# - dropzeroes: if TRUE, columns with <5 non-zero entries are dropped
# - usecutbasis: if FALSE, then the basis is not cut and a standard spline basis is returned
# - useSigma: if TRUE, ytilde= Sigma^{-1/2} y and Sigma^{-1/2} w are returned
# - Sigma: covariance matrix, as given to localnulltest
# Output
# - w: design matrix. Its first columns corresponds to baseline design matrix w0, the rest to the cut spline basis w1. The output satisfies cor(w0,w1)=0
# - wnew: design matrix for (xnew, znew). Note that cor(w0new, w1new) need not be zero
# - ncolw0, ncolw1: number of columns in w0 and w1
# - vargroups: variable groups to be used in modelSelection, indicating variables which should be added/dropped as a group. For example, if there are multiple columns in w1 coding for a local effect, these are assigned to the same group
# - vargroupsn: same as vargroups, using numerical instead of text values
# - w1varname: name of the covariate associated to each column in w1
# - ytilde: block-diag(Sigma)^{-1/2} y
# - wtilde: block-diag(Sigma)^{-1/2} w
# - wtilde.adjust: block-diag(Sigma)^{-1/2} x.adjust
# - logdetSinv: log-determinant of block-diag(Sigma)^{-1}
createDesignLocaltest= function(x, z, y, x.adjust, xnew, znew, function_id, region, regionnew, regionbounds, basedegree, cutdegree, knots, dropzeroes=TRUE, usecutbasis=TRUE, useSigma=FALSE, Sigma) {
    if (!missing(xnew)) {
        xall= rbind(x, xnew)
        zall= rbind(z, znew)
        regionall= c(region, regionnew)
    } else {
        xall= x
        zall= z
        regionall= region
    }
    if (ncol(z)==1) {
       w0= bspline(zall, degree=basedegree, knots=knots[[1]])
    } else {
       w0= tensorbspline(zall, degree=basedegree, knots=knots, maineffects=TRUE)
       w0= cbind(w0$main, w0$tensor)
    }
    w0= w0[,apply(w0,2,'sd')>0]  #remove columns with no observations
    wcut= cutbasis(zall, degree=cutdegree, region=regionall, knots=regionbounds, dropzeroes=dropzeroes, usecutbasis=usecutbasis)
    w1= vargroups= vector("list",ncol(x))
    for (j in 1:ncol(x)) {
        w1[[j]]= xall[,j] * wcut$basis
        vargroups[[j]]= paste('x',j,'.R',wcut$regionid,sep='')
    }
    w1varname= rep(paste('x',1:ncol(x),sep=''), sapply(w1, ncol))
    w1= do.call(cbind, w1)
    vargroups= do.call(c, vargroups)
    vargroupsn= c(1:ncol(w0), ncol(w0) + as.numeric(factor(vargroups)))
    vargroups= c(paste("baseline",1:ncol(w0),sep=''), vargroups)
    # Separate w from wnew
    if (!missing(xnew)) {
        sel= 1:nrow(x)
        w0new= w0[-sel,]; w1new= w1[-sel,]
        w0= w0[sel,]; w1= w1[sel,]
    } else {
        w0new= w1new= NULL
    }
    # Orthogonalize design matrix coding for local covariate effects
    w_decorrelated= decorrelate(w0=w0, w1=w1, region=region, w0new=w0new, w1new=w1new, regionnew=regionnew)
    # Remove columns that became constant after removing xnew
    w1o= w_decorrelated$w1o
    colsel0= (apply(w0, 2, 'sd') > 0); colsel1= (apply(w1o, 2, 'sd') > 0)
    w0= w0[,colsel0,drop=FALSE]; w1o= w1o[,colsel1,drop=FALSE]
    vargroups= vargroups[c(colsel0,colsel1)]
    vargroupsn= vargroupsn[c(colsel0,colsel1)]
    w1varname= w1varname[colsel1]
    w= cbind(w0, w1o)
    #w= cbind(w0, w_decorrelated$w1o)
    if (!missing(xnew)) {
        wnew= cbind(w0new[,colsel0,drop=FALSE], w_decorrelated$w1newo[,colsel1,drop=FALSE])
    } else {
        wnew= NULL
    }
    #If errors are dependent, transform the outcome and the design matrix
    if (useSigma) {
        if (!is.function(Sigma)) {
            S= estimateSigma_localtest(y=y, w=w, z=z, function_id=function_id, region=region, method=Sigma)
        } else {
            zregion= unique(cbind(as.numeric(region), z))
            S= evalSigma_localtest(Sigma, region=zregion[,1], z=zregion[,-1,drop=FALSE])
        }
        yt= uncorrelate_outcome(y=y, w=w, x.adjust=x.adjust, region=region, S=S, function_id=function_id) #return Sigma^{-1/2} y, Sigma^{-1/2} w, log |Sinv|
        ytilde= yt$ytilde; wtilde= yt$wtilde; logdetSinv= yt$logdetSinv
        wtilde.adjust= w_decorrelated$wtilde.adjust
    } else {
        ytilde= wtilde= NULL
        logdetSinv= 0
    }
    ans= list(w=w, wnew=wnew, ncolw0= ncol(w0), ncolw1= ncol(w_decorrelated$w1o), vargroups=vargroups, vargroupsn=vargroupsn, w1varname=w1varname, ytilde=ytilde, wtilde=wtilde, logdetSinv=logdetSinv)
    return(ans)
}




cutbasis= function(z, degree, region, knots, dropzeroes=TRUE, usecutbasis=TRUE) {
    #For each region, create a cut B-spline basis for z that is zero when z is outside the region. Tensor B-splines are used when ncol(z)>1
    #Input
    # - z: matrix for which cut spline basis is to be created
    # - degree: spline degree
    # - region: vector of length == nrow(z), indicating the region of each row in z
    # - knots: list of length == ncol(z) indicating the knots for each dimension of z
    # - dropzeroes: if TRUE, columns that have <4 non-zero elements are dropped
    # - usecutbasis: if FALSE, then the basis is not cut and a standard spline basis is returned
    #Output
    # - basis: the cut basis
    # - regionid: indicates the region that each column in basis corresponds to
    if (length(region) != nrow(z)) stop("length(region) must be equal to nrow(z)")
    #Obtain standard B-splines
    if (ncol(z)==1) {
       baseline= bspline(z, degree=degree, knots=knots[[1]])
    } else {
       baseline= tensorbspline(z, degree=degree, knots=knots, maineffects=FALSE)$tensor
    }
    #Obtain cut splines by placing 0's in the standard B-splines
    regionids= unique(region)
    regionids= regionids[order(as.numeric(regionids))]
    nregions= length(regionids)
    basis= vector("list",nregions)
    for (i in 1:nregions) {
        rowsel= (region == regionids[i])
        basissel= (colSums(baseline[rowsel,,drop=FALSE] != 0) >= ifelse(dropzeroes,5,1)) #remove columns that are essentially always 0
        if (usecutbasis) {
            basis[[i]]= matrix(0, nrow=nrow(z), ncol=sum(basissel))
            basis[[i]][rowsel,]= baseline[rowsel,basissel]
        } else {
            basis[[i]]= baseline[,basissel,drop=FALSE]
        }
        if (sum(basissel)==1) {
            colnames(basis[[i]])= paste("R",regionids[i],sep='')
        } else if (sum(basissel)>1) {
            colnames(basis[[i]])= paste("R",regionids[i],".",1:ncol(basis[[i]]),sep='')
        }
    }
    regionid= rep(regionids, sapply(basis,ncol))
    basis= do.call(cbind, basis)
    ans= list(basis=basis, regionid=regionid)
    return(ans)    
}



cutbasisuniv= function(z, degree, regionbounds) {
    #Same as cutbasis, when z is a vector (i.e. one-dimensional coordinates)
    if (!is.vector(z)) stop("z must be a vector")
    regioninterval= cbind(regionbounds[-length(regionbounds)], regionbounds[-1])
    regioninterval= regioninterval[(regioninterval[,2] > min(z)) & (regioninterval[,1] < max(z)),] #restrict regions to range of observed z's
    nregions= nrow(regioninterval)
    #Range of z values spanned by each standard B-spline
    basisbounds= matrix(NA, nrow=length(regionbounds)-degree-1, ncol=2)
    for (i in 1:nrow(basisbounds)) basisbounds[i,]= c(regionbounds[i],regionbounds[i+degree+1])
    #Obtain standard B-splines
    baseline= bspline(z, degree=degree, knots=regionbounds)
    #Obtain cut splines by placing 0's in the standard B-splines
    basis= matrix(0, nrow=length(z), ncol= nregions * (degree+1))
    subsetid= integer(ncol(basis))
    for (i in 1:nregions) {
        rowsel= ((z >= regioninterval[i,1]) & (z < regioninterval[i,2]))
        colsel= (1 + (i-1)*(degree+1)): (i*(degree+1))
        basissel= (basisbounds[,1] < regioninterval[i,2]) & (basisbounds[,2] > regioninterval[i,1])
        if (length(colsel) > sum(basissel)) colsel= colsel[1:sum(basissel)]
        basis[rowsel,colsel]= baseline[rowsel, basissel]
        subsetid[colsel]= i
    }
    nonzero= (colSums(basis != 0) >= 5) #remove columns that are essentially always 0
    basis= basis[,nonzero]
    subsetid= subsetid[nonzero]
    regionwithid= data.frame(1:nrow(regioninterval), regioninterval)
    names(regionwithid)= c('subsetid','regionstart','regionend')
    subsetid= merge(data.frame(subsetid=subsetid), regionwithid, by='subsetid')
    ans= list(basis=basis, subsetid=subsetid, regions=regioninterval)
    return(ans)    
}



decorrelate= function(w0, w1, region, w0new, w1new, regionnew) {
    # Decorrelate basis w1 coding for local covariate effects in each region from basis w0 coding for the baseline mean
    # Input
    # - w0: design matrix for baseline
    # - w1: design matrix for effect of covariates
    # - region: vector indicating the region that each row in (w0,w1) corresponds to
    # - w0new: design matrix for baseline for new observations
    # - w1new: design matrix for effect of covariates for new observations
    # - regionnew: region that each row in (w0new, w1new) corresponds to
    # Ouput:
    # - w1o: residuals from regressing w1 onto w0
    # - w1newo: residuals from predicting w1new from w0new, using the least-squares coefficient regressing w1 onto w0
    # Note: entries in cor(w0, w1o) are guaranteed to all be zero, but not so for (w0new, w1newo)
    if (nrow(w0) != nrow(w1)) stop("w0 and w1 must have the same number of rows")
    if (nrow(w0) != length(region)) stop("region must be a vector of length equal to nrow(w0)")
    w1o= matrix(0, nrow=nrow(w1), ncol=ncol(w1))
    if (!is.null(w0new)) {
        if (nrow(w0new) != nrow(w1new)) stop("w0new and w1new must have the same number of rows")
        if (nrow(w0new) != length(regionnew)) stop("regionnew must be a vector of length equal to nrow(w0new)")
        w1newo= matrix(0, nrow=nrow(w1new), ncol=ncol(w1new))
    } else {
        w1newo= NULL
    }
    regionids= unique(region)
    for (j in 1:length(regionids)) {
        rowsel= (region == regionids[j])     #rows in w0 corresponding to observations in region j
        colsel1= (colSums(w1[rowsel,,drop=FALSE]!=0) > 0) #columns in w1 coding for region j
        colsel0= (colSums(w0[rowsel,,drop=FALSE]!=0) > 0) #columns in w0 coding for region j
        if (any(colsel1)) {
            fit= lm(w1[rowsel,colsel1] ~ w0[rowsel,colsel0])
            w1o[rowsel,colsel1]= residuals(fit)
            if (!is.null(w0new)) {
                rowselnew= (regionnew == regionids[j])   #rows in w0new corresponding to observations in region j
                if (any(rowselnew)) {
                    b= coef(fit)
                    b[is.na(b)]= 0
                    w1newo[rowselnew,colsel1] = w1new[rowselnew,colsel1,drop=FALSE] - cbind(rep(1,sum(rowselnew)), w0new[rowselnew,colsel0,drop=FALSE]) %*% b
                }
            }
        }
    }
    ans= list(w1o=w1o, w1newo=w1newo)
    return(ans)
}


#For each local test, find the regions that have an overlap with the test
# Input
# - localgrid: list where each element corresponds to a coordinate, and contains the grid defining the local test boundaries in each coordinate
# - regioncoord: list where each element corresponds to a coordinate, and contains the region start and end formatted as a data.frame
# Output
# - testIntervals: localgrid formatted as intervals
# - testxregion: list with an entry for each test, indicating the regions overlapping that test
testOverlaps= function(localgrid, regioncoord) {
    ndim= length(localgrid)
    testIntervals= overlaps= vector("list", length(localgrid))
    #For each coordinate, find overlaps between local tests & regions
    for (i in 1:ndim) {
        testIntervals[[i]]= Intervals(cbind(localgrid[[i]][-length(localgrid[[i]])], localgrid[[i]][-1])) #from package intervals
        regionIntervals= Intervals(regioncoord[[i]]) 
        overlaps[[i]]= interval_overlap(testIntervals[[i]], regionIntervals) #from package intervals
    }
    #For each local test, select regions that overlap in all coordinates
    testxregion= overlaps[[1]]
    if (ndim > 1) {
      tests= expand.grid(lapply(localgrid, function(l) 1:(length(l)-1)))
      names(tests)= paste('dim',1:ncol(tests),sep='')
      for (j in 1:nrow(tests)) {
          testjoverlap= vector("list",ndim)
          for (i in 2:ndim) {
              testdimi= tests[j,i]
              testjoverlap[[i]]= overlaps[[i]][[testdimi]]  #regions overlapping j^th test in i^th coordinate
              testxregion[[j]]= intersect(testxregion[[j]], testjoverlap[[i]])
          }
      }
    }
    #Convert region ids to names
    regionnames= rownames(regioncoord[[1]])
    testxregion= lapply(testxregion, function(zz) regionnames[zz])
    #Return output
    ans= list(testIntervals=testIntervals, testxregion=testxregion)
    return(ans)
}


#For each local test, compute its posterior probability as the prob of any region overlapping the test being active
# Input
# - regionids: binary matrix with 1 column per region
# - ppregions: posterior probability of each row in regionids
# - testxregion: list with an entry for each local test, indicating the regions overlapping that test. These region labels must match column names in regionids
postProblocaltest= function(regionids, ppregions, testxregion) {
    ntests= length(testxregion)
    ans= double(ntests)
    for (i in 1:ntests) {
        sel= colnames(regionids) %in% testxregion[[i]]         #regions overlapping with test i
        ans[i]= sum(ppregions[rowSums(regionids[,sel,drop=FALSE]) > 0])   #sum pp whenever any region is active
    }
    return(ans)
}




