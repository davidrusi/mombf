

setMethod("getBIC", signature(object='msfit'), function(object) {
    pc= object$priors$priorCoef
    pm= object$priors$priorDelta
    if ((pc@priorDistr != 'bic') || (pm@priorDistr != 'uniform')) stop("To obtain BIC you should set priorCoef=bicprior() and priorDelta=modelunifprior() when calling modelSelection")
    bic = -2 * ( object$postProb + object$p * log(2) )
    ans= cbind(getmodelid(object), bic)
    return(ans)
}
)



setMethod("getEBIC", signature(object='msfit'), function(object) {
    pc= object$priors$priorCoef
    pm= object$priors$priorDelta
    isbbprior= (pm@priorDistr == 'binomial') && all(pm@priorPars == c(1,1))
    if ((pc@priorDistr != 'bic') || !isbbprior) stop("To obtain BIC you should set priorCoef=bicprior() and priorDelta=modelbbprior() when calling modelSelection")
    ebic = -2  * ( object$postProb + log(object$p+1) )
    ans= cbind(getmodelid(object), ebic)
    return(ans)
}
)

