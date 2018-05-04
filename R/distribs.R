#Density of a Dirichlet(q) ditribution evaluated at x
ddir= function (x, q, logscale=TRUE) {
    if (is.vector(x)) x= matrix(x,nrow=1)
    if (length(q)==1) q= rep(q,ncol(x))
    sel= rowSums((x<0) | (x>1))>0
    ans= double(nrow(x))
    ans[sel]= -Inf
    xnorm= t(x[!sel,,drop=FALSE])/rowSums(x[!sel,,drop=FALSE])
    ans[!sel]= colSums((q - 1) * log(xnorm)) + lgamma(sum(q)) - sum(lgamma(q))
    if (!logscale) ans= exp(ans)
    return(ans)
}


#Adapted from LaplacesDemon
dinvgamma= function (x, shape = 1, scale = 1, log = FALSE) {
    x <- as.vector(x)
    shape <- as.vector(shape)
    scale <- as.vector(scale)
    if (any(shape <= 0) | any(scale <= 0)) stop("The shape and scale parameters must be positive.")
    NN <- max(length(x), length(shape), length(scale))
    x <- rep(x, len = NN)
    shape <- rep(shape, len = NN)
    scale <- rep(scale, len = NN)
    alpha <- shape
    beta <- scale
    ans <- alpha * log(beta) - lgamma(alpha) - {alpha + 1} * log(x) - {beta/x}
    if (log == FALSE) ans <- exp(ans)
    return(ans)
}
