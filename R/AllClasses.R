###
### AllClasses.R
###

require(methods)

##=============================================================================
setClass("msPriorSpec",
         representation(priorType="character",
                        priorDistr="character",
                        priorPars="vector"),
         prototype(priorPars=NA))


setValidity("msPriorSpec", function(object){
  msg <- NULL
  if (!any(object@priorType %in% c('coefficients','modelIndicator','nuisancePars'))) {
    msg <- "priorType must be 'coefficients', 'modelIndicator' or 'nuisancePars'"
  } else {
    if (object@priorType=='coefficients') {
      
      if (!any(object@priorDistr %in% c('pMOM','piMOM','peMOM'))) {
        msg <- "priorDistr must be 'pMOM','piMOM' or 'peMOM'"
      } else {
        if (object@priorDistr=='pMOM') { if (!('r' %in% names(object@priorPars))) msg <- "priorPars must contain an element named 'r'" } 
        if (!('tau' %in% names(object@priorPars)) & !all(c('a.tau','b.tau') %in% names(object@priorPars))) msg <- "priorPars must either specify 'tau' or 'a.tau' and 'b.tau'"        }
      
    } else if (object@priorType=='modelIndicator') {

      if (!any(object@priorDistr %in% c('uniform','binomial'))) {
        msg <- "priorDistr must be 'uniform' or 'binomial'"
      } else {
        if (object@priorDistr=='binomial') {
          n <- c('p','alpha.p','beta.p') %in% names(object@priorPars)
          if ((!n[1]) & (!all(n[-1]))) msg <- "For priorDistr=='binomial' either 'p' or 'alpha.p' and 'beta.p' must be specified in priorPars"
        }
      }
      
    } else {
      if (object@priorDistr!='invgamma') {
        msg <- "priorDistr must be invgamma"
      } else {
        if (!all(c('alpha','lambda') %in% names(object@priorPars))) msg <- "priorPars must contain elements named 'alpha', 'lambda'"
      }
    }
  }

  ifelse(is.null(msg),TRUE,msg)
}
)



            
