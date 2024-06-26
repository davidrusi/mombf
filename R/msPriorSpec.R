###
### msPriorSpec.R
###


##-----------------------------------------------------------------------------
## Invoked by validObject() method.
valid_msPriorSpec <- function(object) {

    msg <- NULL

    ## Ensure first two slots are scalar
    if (length(object@priorType) != 1) {
        msg <- c(msg,
                 sprintf("slot %s must be of length 1",
                         sQuote("priorType")))
        object@priorType <- object@priorType[1]         # process first only
    }
    if (length(object@priorDistr) != 1) {
        msg <- c(msg,
                 sprintf("slot %s must be of length 1",
                         sQuote("priorDistr")))
        object@priorDistr <- object@priorDistr[1]       # process first only
    }

    ## Ensure valid type
    valid_prior_types <- c("coefficients",
                           "modelIndicator",
                           "nuisancePars")
    found <- object@priorType %in% valid_prior_types
    if (!found) {
        msg <- c(msg,
                 sprintf("slot %s must be one of: %s",
                         sQuote("priorType"),
                         paste(dQuote(valid_prior_types), collapse=", ")))
    }

    ##-------------------------------------------------------------------------
    has_either_or <- function(this, that, among) {
        this_in <- this %in% among
        that_in <- that %in% among

        (all(this_in) & !any(that_in)) || (!any(this_in) & all(that_in))
    }


    ## Validate rest based on type
    switch(object@priorType,
           coefficients = {
               valid_coef_prior_distrs <- c(
                 "pMOM","groupMOM","piMOM","peMOM","zellner","groupzellner","normalid","bic","icarplus"
               )
               found <- object@priorDistr %in% valid_coef_prior_distrs
               if (!found) {
                   msg <- c(msg,
                            sprintf("slot %s must be one of: %s",
                                    sQuote("priorDistr"),
                                    paste(dQuote(valid_coef_prior_distrs),
                                          collapse=", ")))
               }
               if (object@priorDistr == "pMOM") {
                   reqd_pmom_prior_par <- "r"
                   found <- reqd_pmom_prior_par %in% names(object@priorPars)
                   if (!found) {
                       msg <- c(msg,
                                sprintf("slot %s must be vector with named element %s (when using %s distr)",
                                        sQuote("priorPars"),
                                        dQuote(reqd_pmom_prior_par),
                                        sQuote("pMOM")))
                   }
               }
               prPars_names <- names(object@priorPars)
               found <- (
                  (("a.tau" %in% prPars_names) & ("b.tau" %in% prPars_names)) ||
                  ("tau" %in% prPars_names) ||
                  ("taustd" %in% prPars_names)
               )
               if (!found) {
                   msg <- c(msg,
                            sprintf("slot %s must be vector with named element(s) tau, taustd or, a.tau and b.tau",
                                    sQuote("priorPars")))
               }
           },
           modelIndicator = {
               valid_ind_prior_distrs <- c("uniform",
                                           "binomial",
                                           "complexity")
               found <- object@priorDistr %in% valid_ind_prior_distrs
               if (!found) {
                   msg <- c(msg,
                            sprintf("slot %s must be one of: %s",
                                    sQuote("priorDistr"),
                                    paste(dQuote(valid_ind_prior_distrs), collapse=", ")))
               }
               if (object@priorDistr == "binomial") {
                   reqd_prior_pars <- c("p",
                                        "alpha.p",
                                        "beta.p")
                   found <- has_either_or(reqd_prior_pars[1],
                                          reqd_prior_pars[-1],
                                          names(object@priorPars))
                   if (!found) {
                       msg <- c(msg,
                                sprintf("slot %s must be vector with named element(s) %s, or %s (when using %s distr)",
                                        sQuote("priorPars"),
                                        dQuote(reqd_prior_pars[1]),
                                        paste(dQuote(reqd_prior_pars[-1]),
                                              collapse=" and "),
                                        sQuote("binomial")))
                   }
               } else if (object@priorDistr == "complexity") {
                   reqd_prior_pars <- "c"
                   found <- reqd_prior_pars %in% names(object@priorPars)
                   if (!found) {
                       msg <- c(msg,
                                sprintf("slot %s must be vector with names element %s",
                                        sQuote("priorPars"),
                                        dQuote(reqd_prior_pars[1])))
                   }
               }
           },
           nuisancePars = {
               valid_nuisance_prior_distrs <- c("gamma","invgamma")
               found <- object@priorDistr %in% valid_nuisance_prior_distrs
               if (!found) {
                 msg <- c(msg,
                          sprintf("slot %s must be %s",
                                  sQuote("priorDistr"),
                                  dQuote(valid_nuisance_prior_distrs)))
               }
               reqd_prior_pars <- c("alpha",
                                    "lambda")
               found <- reqd_prior_pars %in% names(object@priorPars)
               if (!all(found)) {
                   msg <- c(msg,
                            sprintf("slot %s must contain named elements %s",
                                    sQuote("priorPars"),
                                    paste(dQuote(reqd_prior_pars),
                                          collapse=", ")))
               }
           })

    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
}

setValidity("msPriorSpec", valid_msPriorSpec)


##-----------------------------------------------------------------------------
is.msPriorSpec <- function(x) {
    extends(class(x), "msPriorSpec")
}


##-----------------------------------------------------------------------------
## Generators for this object.
msPriorSpec <- function(priorType=c("coefficients",
                                    "modelIndicator",
                                    "nuisancePars"),
                        priorDistr,
                        priorPars) {
    ## Check arguments
    stopifnot(is.character(priorType)  && length(priorType) == 1)
    stopifnot(is.character(priorDistr) && length(priorDistr) == 1)
    stopifnot(is.vector(priorPars))

    ## Create new class
    new("msPriorSpec",
        priorType=match.arg(priorType),
        priorDistr=priorDistr,
        priorPars=priorPars)
}


bic <- function() {
    new("msPriorSpec", priorType="coefficients", priorDistr="bic", priorPars=c(tau=-1))
}

bicprior <- function() {
    bic()
}

aic <- function() {
    ic(2)
}

ic <- function(penalty) {
    new("msPriorSpec", priorType="coefficients", priorDistr="bic", priorPars=c(tau=penalty))
}


exponentialprior <- function(lambda=1) {
    gammaprior(alpha=1, lambda=lambda)
}

gammaprior <- function(alpha=.01, lambda=.01) {
    new("msPriorSpec", priorType='nuisancePars', priorDistr='gamma', priorPars=c(alpha=alpha,lambda=lambda))
}


igprior <- function(alpha=.01, lambda=.01) {
    new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=alpha,lambda=lambda))
}

momprior <- function(taustd=1, tau, tau.adj=10^6, r=1) {
  if (missing(tau)) {
    priorPars <- c(taustd=taustd, tau.adj=tau.adj,r=r)
  } else {
    priorPars <- c(tau=tau,tau.adj=tau.adj,r=r)
  }
  new("msPriorSpec", priorType="coefficients", priorDistr="pMOM", priorPars=priorPars)
}

icarplusprior <- function(a=0.5, taustd=1, tau.adj=10^6) {
    if ((a < 0) | (a>1)) stop("a must be in [0,1]")
    new("msPriorSpec", priorType="coefficients", priorDistr="icarplus", priorPars=c(a=a, taustd=taustd, tau.adj=tau.adj))
}


imomprior <- function(tau, tau.adj=10^6) {
    new("msPriorSpec", priorType="coefficients", priorDistr="piMOM", priorPars=c(tau=tau, tau.adj=tau.adj))
}

emomprior <- function(tau, tau.adj=10^6) {
    new("msPriorSpec", priorType="coefficients", priorDistr="peMOM", priorPars=c(tau=tau, tau.adj=tau.adj))
}

zellnerprior <- function(taustd=1, tau, tau.adj=10^6) {
  if (missing(tau)) {
    priorPars <- c(taustd=taustd, tau.adj=tau.adj)
  } else {
    priorPars <- c(tau=tau,tau.adj=tau.adj)
  }
  new("msPriorSpec", priorType="coefficients", priorDistr="zellner", priorPars=priorPars)
}

normalidprior <- function(taustd=1, tau, tau.adj=10^6) {
  if (missing(tau)) {
    priorPars <- c(taustd=taustd, tau.adj=tau.adj)
  } else {
    priorPars <- c(tau=tau,tau.adj=tau.adj)
  }
  new("msPriorSpec", priorType="coefficients", priorDistr="normalid", priorPars=priorPars)
}

groupmomprior <- function(taustd=1, tau, tau.adj=10^6) {
  if (missing(tau)) {
    priorPars <- c(taustd=taustd, tau.adj=tau.adj)
  } else {
    priorPars <- c(tau=tau,tau.adj=tau.adj)
  }
  new("msPriorSpec", priorType="coefficients", priorDistr="groupMOM", priorPars=priorPars)
}

groupimomprior <- function(tau, tau.adj=10^6) {
    new("msPriorSpec", priorType="coefficients", priorDistr="groupiMOM", priorPars=c(tau=tau, tau.adj=tau.adj))
}

groupemomprior <- function(tau, tau.adj=10^6) {
    new("msPriorSpec", priorType="coefficients", priorDistr="groupeMOM", priorPars=c(tau=tau, tau.adj=tau.adj))
}

groupzellnerprior <- function(taustd=1, tau, tau.adj=10^6) {
  if (missing(tau)) {
    priorPars <- c(taustd=taustd, tau.adj=tau.adj)
  } else {
    priorPars <- c(tau=tau,tau.adj=tau.adj)
  }
  new("msPriorSpec", priorType="coefficients", priorDistr="groupzellner", priorPars=priorPars)
}


modelunifprior <- function() {
    new("msPriorSpec",priorType='modelIndicator',priorDistr='uniform',priorPars=double(0))
}

modelbinomprior <- function(p=0.5) {
    new("msPriorSpec",priorType='modelIndicator',priorDistr='binomial',priorPars=list(p=p))
}

modelbbprior <- function(alpha.p=1, beta.p=1) {
    new("msPriorSpec",priorType='modelIndicator',priorDistr='binomial',priorPars=c(alpha.p=alpha.p,beta.p=beta.p))
}

modelcomplexprior <- function(c=1) {
    new("msPriorSpec",priorType='modelIndicator',priorDistr='complexity',priorPars=c(c=c))
}


igprior <- function(alpha=.01, lambda=.01) {
    new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=alpha,lambda=lambda))
}
