### Methods for icfit objects

setMethod("show", signature(object='icfit'), function(object) {
  cat("icfit object\n\n")
  cat("Model with best",object$criterion,":",object$topmodel,"\n\n")
  cat("Use summary(), coef() and predict() to get inference for the top model\n")
  cat("Use coef(object$msfit) and predict(object$msfit) to get BMA estimates and predictions\n")
}
)

confint.icfit <- function(object, ...) {
    confint(object$topmodel.fit)
}

coef.icfit <- function(object,...) {
    coef(object$topmodel.fit)
}


predict.icfit <- function(object, ...) {
    predict(object$topmodel.fit, ...)
}

summary.icfit <- function(object, ...) {
    summary(object$topmodel.fit, ...)
}


## FIND THE MODEL ATTAINING THE BEST VALUE OF AN INFORMATION CRITERION

checkargs_IC= function(...) {
    params= eval(quote(list(...)))
    forbiddenpars= c('priorCoef','priorDelta','center','scale')
    if (any(forbiddenpars %in% names(params))) stop(paste("Arguments",paste(forbiddenpars,collapse=", "),"have set values so they cannot be passed on to modelSelection"))
}


topmodelnames= function(ans) {
    topvarids= as.numeric(strsplit(ans$models$modelid[1], ',')[[1]])
    ans= list(topvarids= topvarids, varnames= ans$varnames[topvarids])
    return(ans)
}

family2glm= function(family) {
    if (family=="normal") {
        ans= gaussian()
    } else if (family=="binomial") {
        ans= binomial()
    } else if (family=="poisson") {
        ans= poisson()
    } else stop("Only the Normal, Binomial and Poisson families are currently implemented")
    return(ans)
}


#Extract info criteria from an msfit and return icfit object
extractmsIC= function(ms, getICfun) {
    ans= vector("list",5)
    names(ans)= c('topmodel','topmodel.fit','models','varnames','msfit')
    ans$models= getICfun(ms)
    if (!is.null(colnames(ms$xstd))) {
        ans$varnames= colnames(ms$xstd)
    } else {
        ans$varnames= paste('x[,',1:ncol(ms$xstd),']',sep='')
    }
    tm= topmodelnames(ans)
    ans$topmodel= tm$varnames
    ans$msfit= ms
    data= data.frame(y=ms$ystd, ms$xstd[,tm$topvarids])
    names(data)[-1]= ans$varnames[tm$topvarids]
    f= formula(y ~ -1 + .)
    ans$topmodel.fit= glm(f , data=data, family=family2glm(ms$family))
    new("icfit",ans)
}


bestBIC= function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=bic(), priorDelta=modelunifprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getBIC)
    ans$criterion= 'BIC'
    return(ans)
}

bestAIC= function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=aic(), priorDelta=modelunifprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getAIC)
    ans$criterion= 'AIC'
    return(ans)
}

bestEBIC= function(...) {
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=bic(), priorDelta=modelbbprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getEBIC)
    ans$criterion= 'EBIC'
    return(ans)
}

bestIC= function(..., penalty) {
    if (missing(penalty)) stop("penalty must be specified. Alternatively consider using bestBIC(), bestEBIC() or bestAIC()")
    checkargs_IC(...)
    ms= modelSelection(..., priorCoef=ic(penalty), priorDelta=modelunifprior(), center=FALSE, scale=FALSE)
    ans= extractmsIC(ms, getIC)
    ans$criterion= paste('GIC (penalty=',penalty,')',collapse='')
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



