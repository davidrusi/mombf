###
### AllGenerics.R
###

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

