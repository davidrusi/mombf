###
### AllGenerics.R
###

if (!isGeneric("coefByModel")) {
  setGeneric("coefByModel", function(object, maxmodels, alpha=0.05, niter=10^3, burnin=round(niter/10)) standardGeneric("coefByModel"))
}

if (!isGeneric("coefOneModel")) {
  setGeneric("coefOneModel", function(y, x, m, V, outcometype, family, priorCoef, priorGroup, priorVar, alpha=0.05, niter=10^3, burnin=round(niter/10)) standardGeneric("coefOneModel"))
}


if (!isGeneric("dpimom")) {
    setGeneric("dpimom",
               function(x,
                        tau=1,
                        phi=1,
                        logscale=FALSE) standardGeneric("dpimom"))

}

if (!isGeneric("dpmom")) {
    setGeneric("dpmom",
               function(x,
                        tau,
                        a.tau,
                        b.tau,
                        phi=1,
                        r=1,
                        baseDensity='normal',
                        logscale=FALSE) standardGeneric("dpmom"))
}

if (!isGeneric("demom")) {
    setGeneric("demom",
               function(x,
                        tau,
                        a.tau,
                        b.tau,
                        phi=1,
                        logscale=FALSE) standardGeneric("demom"))
}


if (!isGeneric("getAIC")) {
  setGeneric("getAIC", function(object) standardGeneric("getAIC"))
}

if (!isGeneric("getBIC")) {
  setGeneric("getBIC", function(object) standardGeneric("getBIC"))
}

if (!isGeneric("getEBIC")) {
  setGeneric("getEBIC", function(object) standardGeneric("getEBIC"))
}

if (!isGeneric("getIC")) {
  setGeneric("getIC", function(object) standardGeneric("getIC"))
}


if (!isGeneric("marginalIW")) {
  setGeneric("marginalNIW", function(x, xbar, samplecov, n, z, g,  mu0=rep(0,ncol(x)), nu0=ncol(x)+4, S0,logscale=TRUE) standardGeneric("marginalNIW"))
}

if (!isGeneric("plotprior")) {
  setGeneric("plotprior", function(object, xlab, ylab, ylim=c(0,1), ...) standardGeneric("plotprior"))
}

if (!isGeneric("logjoint")) {
  setGeneric("logjoint", function(object, return_models=TRUE, models_as_char=FALSE) standardGeneric("logjoint"))
}

if (!isGeneric("postProb")) {
  setGeneric("postProb", function(object, nmax, method='norm') standardGeneric("postProb"))
}

if (!isGeneric("postProbSubset")) {
  setGeneric("postProbSubset", function(object, varsubset, nmax, method='norm') standardGeneric("postProbSubset"))
}

if (!isGeneric("marglhood_acrossmodels")) {
  setGeneric("marglhood_acrossmodels", function(object, logscale=TRUE) standardGeneric("marglhood_acrossmodels"))
}

if (!isGeneric("postSamples")) {
    setGeneric("postSamples", function(object) standardGeneric("postSamples"))
}



if (!isGeneric("rnlp")) {
  setGeneric("rnlp", function(y, x, m, V, msfit, outcometype, family, priorCoef, priorGroup, priorVar, priorprec, isgroup, niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') standardGeneric("rnlp"))
}

