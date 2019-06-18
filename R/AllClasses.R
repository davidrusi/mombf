###
### AllClasses.R
###

#require(methods)


##=============================================================================
setClass("mixturebf",
         representation(postprob="data.frame",
                        p='numeric',
                        n='numeric',
                        priorpars="list",
                        postpars="list",
                        mcmc="list"),
         prototype(priorpars=list(), postpars=list(), mcmc=list()))


##=============================================================================
setClass("msPriorSpec",
         representation(priorType="character",
                        priorDistr="character",
                        priorPars="vector"),
         prototype(priorPars=NA))


##=============================================================================
setClass("msfit", representation("list"), prototype = prototype(elementType = "list"), contains="list")
