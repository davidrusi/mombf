## FIND THE MODEL ATTAINING THE BEST VALUE OF AN INFORMATION CRITERION

checkargs_IC= function(...) {
    params= eval(quote(list(...)))
    forbiddenpars= c('priorCoef','priorDelta','center','scale')
    if (any(forbiddenpars %in% names(params))) stop(paste("Arguments",paste(forbiddenpars,collapse=", "),"have set values so they cannot be passed on to modelSelection"))
}


topmodelnames= function(ans) {
    topvarids= strsplit(ans$models$modelid[1], ',')[[1]]
    ans$varnames[as.numeric(topvarids)]
}

bestBIC= function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=bic(), priorDelta=modelunifprior(), center=FALSE, scale=FALSE)
    ans= vector("list",3)
    names(ans)= c('topmodel','models','varnames')
    ans$models= getBIC(ms)
    ans$varnames= colnames(ms$xstd)
    ans$topmodel= topmodelnames(ans)
    return(ans)
}

bestAIC= function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=aic(), priorDelta=modelunifprior(), center=FALSE, scale=FALSE)
    ans= vector("list",3)
    names(ans)= c('topmodel','models','varnames')
    ans$models= getAIC(ms)
    ans$varnames= colnames(ms$xstd)
    ans$topmodel= topmodelnames(ans)
    return(ans)
}

bestEBIC= function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=bic(), priorDelta=modelbbprior(), center=FALSE, scale=FALSE)
    ans= vector("list",3)
    names(ans)= c('topmodel','models','varnames')
    ans$models= getEBIC(ms)
    ans$varnames= colnames(ms$xstd)
    ans$topmodel= topmodelnames(ans)
    return(ans)
}

bestIC= function(..., penalty) {
    if (missing(penalty)) stop("penalty must be specified. Alternatively consider using bestBIC(), bestEBIC() or bestAIC()")
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=ic(penalty), priorDelta=modelunifprior(), center=FALSE, scale=FALSE)
    ans= vector("list",3)
    names(ans)= c('topmodel','models','varnames')
    ans$models= getIC(ms)
    ans$varnames= colnames(ms$xstd)
    ans$topmodel= topmodelnames(ans)
    return(ans)
}



## EXTRACT INFORMATION CRITERIA FROM AN msfit object

setMethod("getAIC", signature(object='msfit'), function(object) {
    ans= getBIC(object)
    names(ans)[names(ans)=='bic']= 'aic'
    return(ans)
}
)


setMethod("getBIC", signature(object='msfit'), function(object) {
    pc= object$priors$priorCoef
    pm= object$priors$priorDelta
    if ((pc@priorDistr != 'bic') || (pm@priorDistr != 'uniform')) stop("To obtain BIC you should set priorCoef=bic() and priorDelta=modelunifprior() when calling modelSelection")
    ans= getIC(object)
    names(ans)[names(ans)=='ic']= 'bic'
    return(ans)
}
)

setMethod("getIC", signature(object='msfit'), function(object) {
    ic = -2 * ( object$postProb + object$p * log(2) )
    ans= tibble(modelid=object$modelid, ic=ic)
    ans= ans[order(ans$ic),]
    return(ans)
}
)


setMethod("getEBIC", signature(object='msfit'), function(object) {
    pc= object$priors$priorCoef
    pm= object$priors$priorDelta
    isbbprior= (pm@priorDistr == 'binomial') && all(pm@priorPars == c(1,1))
    if ((pc@priorDistr != 'bic') || !isbbprior) stop("To obtain BIC you should set priorCoef=bic() and priorDelta=modelbbprior() when calling modelSelection")
    ebic = -2  * ( object$postProb + log(object$p+1) )
    ans= tibble(modelid=object$modelid, ebic=ebic)
    #ans= cbind(getmodelid(object), ebic)
    ans= ans[order(ans$ebic),]
    return(ans)
}
)



