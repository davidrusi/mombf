## Routines for Cox proportional hazards model

pmomCoxMarginalR <- function(y, x, tau, r=1, method=ifelse(ncol(x)<=10,'Normal','plugin'), logscale=TRUE) {
  #Integrated Cox partial likelihood (Normal approximation) wrt product moment prior
  # - Partial Likelihood: f(y|th)= exp(-.5 (th-thhat)' solve(V) (th-thhat))
  # - Prior proportional to N(th; 0, tau*I) * prod(th^2/tau)^r
  # Input
  # - y: Surv object
  # - x: covariates
  # - tau: prior dispersion
  # - r: prior power parameter
  # - method: set to 'Normal' to compute E(prod(th^2)^i) exactly, to 'plugin' to use prod(E(th)^2). By default 'Normal' is used for up to 10 dimensions
  if (class(y)!="Surv") stop("y must be of class 'Surv'")
  if (r != 1) stop("Only r=1 currently implemented")
  if (missing(x)) {
    p <- 0
  } else {
    p <- ncol(x)
    if (!is.data.frame(x)) x <- data.frame(x)
  }
  if (p==0) {
    ans <- coxph(y ~ 1)$loglik
  } else {
    fit <- coxph(y ~ ., data=x)
    thhat <- matrix(coef(fit),ncol=1); Vinv <- solve(fit$var)
    Sinv <- Vinv + diag(p)/tau; S <- solve(Sinv)
    m <- S %*% Vinv %*% thhat
    ans <- as.numeric(fit$loglik[2] - (p+p/2)*log(tau) - .5*(t(thhat) %*% Vinv %*% thhat - t(m) %*% Sinv %*% m) + .5*as.numeric(determinant(S,logarithm=TRUE)$modulus))
    if (method=='Normal') {
      ans <- ans + log(eprod(m, S, power=2*r, dof=-1))
    } else {
      ans <- ans + prod(m^2+diag(S))
    }
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}



logPL <- function(beta, time, event, x, mc.cores=1) {
  #Cox log-partial likelihood for several values of beta
  # - beta: matrix with ncol(beta)==ncol(x) and rows containing the various beta values. When ncol(x)==1 beta can also be a vector
  # - time: survival times (censored or uncensored)
  # - event: event indicator
  # - x: covariate values
  # - mc.cores: number of cores, passed on to mclapply
  if (mc.cores > 1) {
    if (!require("parallel", character.only=TRUE)) {
      warning("cannot load parallel package - continuing without it...")
    }
  }
  event <- as.numeric(event)
  if (!all(event %in% c(0,1))) stop("event must contain only 0's and 1's")
  if (!is.matrix(x) & !is.data.frame(x)) stop("x must be a matrix or data.frame")
  o <- order(time)
  time <- time[o]; event <- event[o]; x <- x[o,,drop=FALSE]

  if (is.vector(beta)) {
    if (ncol(x)!=1) stop('beta should be a matrix for ncol(x)>1')
    if ("parallel" %in% loadedNamespaces())  {
      ans <- mclapply(1:length(beta), function(i) { logPLSingle(length(time), event=event, x=x, beta[i]) }, mc.cores=mc.cores)
    } else {
      ans <- lapply(1:length(beta), function(i) { logPLSingle(length(time), event=event, x=x, beta[i]) })
    }
  } else {
    if ("parallel" %in% loadedNamespaces())  {
      ans <- mclapply(1:length(beta), function(i) { logPLSingle(length(time), event=event, x=x, beta[i,]) }, mc.cores=mc.cores)
    } else {
      ans <- lapply(1:length(beta), function(i) { logPLSingle(length(time), event=event, x=x, beta[i,]) })
    }
  }
  unlist(ans)
}

logPLSingle <- function(n, event, x, beta) {
  #Cox log-partial likelihood for a single value of beta
  #Note: event, x should be ordered according to the ranks of time
  #Downloaded from: http://www4.stat.ncsu.edu/~lu/research/example-pbc-cox.R
  L = 0.0
  for(i in 1:n){
     temp = 0.0
     for(j in i:n) temp = temp + exp(sum(beta*x[j,])) 
     L = L - event[i]*(sum(beta*x[i,])-log(temp))
  }         
  return(-L)
} 
