plot_taylor_approximation <- function(
    f,          # Original function (takes vector input)
    grad_f,     # Gradient function (returns vector)
    hessian_f,  # Hessian function (returns matrix)
    x0,         # Current point (vector)
    coord_idx,  # Index of coordinate to vary (e.g., 1 for x1)
    x_lim = c(x0[coord_idx] - 1, x0[coord_idx] + 1),  # Range for varying coordinate
    n_points = 100,  # Number of points for plotting
    ...         # Other arguments passed to f, grad_f and hessian_f
) {
  # Check inputs
  if (length(x0) < coord_idx)  stop("coord_idx exceeds the dimension of x0")
  
  # Generate sequence for varying coordinate
  x_var <- seq(x_lim[1], x_lim[2], length.out = n_points)
  
  # Evaluate original function along the varying coordinate
  f_values <- sapply(x_var, function(xi) {
    x <- x0
    x[coord_idx] <- xi
    f(x, ...)
  })
  
  # Compute Taylor approximations at x0
  grad <- grad_f(x0, ...)
  hess <- hessian_f(x0, ...)
  delta <- x_var - x0[coord_idx]
  
  # First-order Taylor approximation (linear)
  taylor1 <- f(x0, ...) + grad[coord_idx] * delta
  
  # Second-order Taylor approximation (quadratic)
  taylor2 <- f(x0, ...) + grad[coord_idx] * delta + 0.5 * hess[coord_idx, coord_idx] * delta^2
  
  # Plot
  plot(x_var, f_values, type = "l", lwd = 2, col = "black",
       xlab = paste("Coordinate", coord_idx, "(x", coord_idx, ")"),
       ylab = "Function value",
       main = paste("Function and Taylor Approximations at x[", coord_idx, "] =", x0[coord_idx]))
  abline(v = x0[coord_idx], lty = 2, col = "gray")  # Mark the current point
  lines(x_var, taylor1, col = "red", lwd = 2, lty = 2)    # 1st-order approx
  lines(x_var, taylor2, col = "blue", lwd = 2, lty = 3)   # 2nd-order approx
  legend("topright", legend = c("Original", "1st-order Taylor", "2nd-order Taylor"),
         col = c("black", "red", "blue"), lwd = 2, lty = c(1, 2, 3))
}



# - Z: p x q matrix containing the q meta-covariates for the p covariates. It should not include an intercept, as it's already automatically added
# - wini: Optional. Initial value for the q hyper-parameters
# - niter.mcmc= number of iterations in the final MCMC, run once after hyper-parameter estimates have been obtained
# - niter.mstep= number of MCMC iterations in each M-step required to update hyper-parameter estimates
# - niter.eBayes: max number of iterations in the empirical Bayes optimization algorithm. The algorithm also stop when the objective function didn't improve in >=2 iterations
# - priorvar.w: hyper-parameters w follow a prior w ~ N(0, w_priorvar I) where priorvar.w is the prior variance
# - verbose: if TRUE, progress in hyper-parameter estimation is printed
modelSelection_eBayes= function(Z, wini, niter.mcmc= 5000, niter.mstep= 1000, niter.eBayes= 20, priorvar.w, verbose=TRUE, ...) {
    # Add intercept and orthogonalize Z wrt intercept
    Z= model.matrix(~ ., data= as.data.frame(Z))
    zerosd= which(apply(Z, 2, 'sd')==0)
    if (length(zerosd) > 1) {
        Z= Z[, -zerosd[-1], drop=FALSE] #if there's >1 intercept, keep only the first
    }
    sel= (colnames(Z) != "(Intercept)")
    if (sum(sel) > 0) Z[,sel]= round(residuals(lm(Z[,sel] ~ 1)), 3)
    # Default hyper-parameters prior variance
    if (missing(priorvar.w)) priorvar.w= eBayes_priorelicit(Z=Z, targetprob=0.95, min_priorprob=0.001, prior="zellner")
    # Get initial parameter estimates
    if (verbose) cat("Initializing hyper-parameter estimates... ")
    if (missing(wini)) {
        # Obtain posterior inclusion prob under Beta-Binomial(1,1) model prior
        ms= modelSelection(priorDelta=modelbbprior(), niter=niter.mstep, verbose=FALSE, ...)
    } else {
        # Obtain posterior inclusion prob under Bernoulli prior for given wini
        priorprob= 1 / (1 + exp(- Z %*% matrix(wini,ncol=1)))
        priorDelta= modelbinomprior(priorprob)
        ms= modelSelection(priorDelta=priorDelta, niter=niter.mstep, verbose=FALSE, ...)
    }
    # Set initial estimate using solution from full-rank/clustering case 
    wini= eBayes_mstep_init(Z=Z, postprob=ms$margpp, priorvar.w=priorvar.w)
    # Iterate
    i = 1
    w= wini$w
    found= FALSE
    Vinv= (t(Z) %*% Z) / nrow(Z)
    priorhess= - Vinv / priorvar.w
    fval= double(niter.eBayes + 1)
    fval[1]= fbest= NA
    noimprove= 0
    while (!found) {
        # E-step: find posterior probabilities for current w
        priorprob= 1 / (1 + exp(- Z %*% w))
        priorDelta= modelbinomprior(priorprob)
        ms= modelSelection(priorDelta=priorDelta, niter=niter.mstep, verbose=FALSE, ...)
        # M-step
        if (wini$fullrankcase) {
            wnew= eBayes_mstep_fullrank(Z=Z, Uinv=wini$Uinv, postprob=ms$margpp, wcur=w, priorvar.w=priorvar.w)
        } else {
            wnew= eBayes_mstep(wcurrent=w, Z=Z, postprob=ms$margpp, maxiter=10, priorvar.w=priorvar.w, Vinv=Vinv, priorhess=priorhess)$w
        }
        maxstep= max(abs(wnew - w))
        # Objective function
        fval[i+1]= eBayes_logit_objective(wnew, msfit=ms, Z=Z, priorvar.w=priorvar.w, Vinv=Vinv)
        #fval[i+1]= eBayes_em_logit_objective(wnew, msfit=ms, Z=Z, priorvar.w=priorvar.w, Vinv=Vinv)
        if (is.na(fbest)) {
            fbest= fval[i]= eBayes_logit_objective(w, msfit=ms, Z=Z, priorvar.w=priorvar.w, Vinv=Vinv)
            #fbest= fval[i]= eBayes_em_logit_objective(w, msfit=ms, Z=Z, priorvar.w=priorvar.w, Vinv=Vinv)
            wbest= w
            if (verbose) cat("Done\n\n","Iteration", "Objective function", "Estimate", "\n", i, w, fval[i], "\n")
        }
        if (fval[i+1] >= fval[i]) { noimprove= 0 } else { noimprove= noimprove + 1 }
        if (fval[i+1] > fbest) { fbest= fval[i+1]; wbest= wnew }
        w= wnew
        found = (i >= niter.eBayes) || (maxstep < 0.01) || (noimprove > 1)
        i= i + 1
        if (verbose) cat(" ", i, w, fval[i], "\n")
    }
    # Run modelSelection with empirical Bayes prior probabilities
    if (verbose) cat("Done. \n","Final MCMC with best hyper-parameter value... ")
    w= wbest
    priorprob= 1 / (1 + exp(- Z %*% w))
    priorDelta= modelbinomprior(priorprob)
    ms= modelSelection(priorDelta=priorDelta, niter=niter.mcmc, verbose=FALSE, ...)
    ms$eBayes_hyperpar= w
    ms$Z= Z
    if (verbose) cat("Done\n")
    return(ms)
}


# Find prior variance of hyper-parameters such that all prior inclusion probabilities are in (min_priorprob, 1 - min_priorprob) with >= targetprob probability
# Hyper-parameters w follow a prior w ~ N(0, priorvar.w V), where V=I when prior == "normalshrinkage" and V= (Z^T Z / nrow(Z))^{-1} when prior == "zellner"
# The prior inclusion probability for a covariate that has meta-covariates z is m(w)= 1 / (1 + e^{-z^T w})
#
# Input
# - Z: p x q matrix containing the q meta-covariates for the p covariates. It should not include an intercept, as it's already automatically added
# - targetprob
# - min_priorprob
# Output: value of priorvar.w such that the prior prob that m(z) >= targetprob for all z values obtained by taking a row from Z
eBayes_priorelicit= function(Z, targetprob=0.95, min_priorprob=0.001, prior="zellner") {
    if (min_priorprob >= 0.5) stop("min_priorprob must be <= 0.5")
    #logit= function(x) log(x / (1-x))
    priorvar.w= (log(min_priorprob/(1-min_priorprob)) / qnorm((1-targetprob)/2))^2
    if (prior == "zellner") {
        V= solve((t(Z) %*% Z) / nrow(Z))
        v= mahalanobis(Z, center=rep(0,ncol(Z)), cov=V, inverted=TRUE)
        vmax= max(v)
    } else if (prior == "normalshrinkage") {
        vmax= max(rowSums(z^2))
    }
    priorvar.w= priorvar.w / vmax
    #unique(apply(Z, 1, function(z) 1 - 2 * pnorm(logit(0.05) / sqrt(priorvar.w * mahalanobis(z, center=rep(0,ncol(Z)), cov=V, inverted=TRUE))))) #check: should all be >= targetprob
    return(priorvar.w)    
}



# Exact M-step solution when the number of unique rows in Z is ncol(Z), and U which is the ncol(Z) x ncol(Z) matrix containing these unique rows has full rank
# The M-step seeks w such that Z^T m(w) = Z^T postprob, where m(z) = 1 / (1 + exp{- Z w})
eBayes_mstep_fullrank= function(Z, Uinv, postprob, wcur, priorvar.w) {
    if (missing(Uinv)) {
        U= unique(Z)
        if (nrow(U) != ncol(U)) stop("The number of unique rows in Z is not ncol(Z)")
        Uinv= solve(U)
    } else {
        if (ncol(Uinv) != ncol(Z)) stop("Uinv must have the same number of columns as Z")
        if (nrow(Uinv) != ncol(Uinv)) stop("Uinv must be a square matrix")
    }
    Ztilde= Z %*% Uinv
    meanpp= (t(Ztilde) %*% postprob) / colSums(Ztilde)
    if (missing(wcur)) wcur= log(meanpp / (1-meanpp))
    wtilde= double(ncol(Z))
    for (j in 1:length(wtilde)) {
        wtilde[j]= eBayes_fullrank_solve_modifiedlogistic(a=meanpp[j], c= priorvar.w * nrow(Z), w_initial=wcur[j])
    }
    what= Uinv %*% wtilde
    #grad= eBayes_em_logit_grad(what, postprob=postprob, Z=Z)  #check: gradient should be 0
    return(what)
}


# Newton-Raphson solver for f(w) = a, where f(w) = 1 / (1 + exp(-w)) + w/c
eBayes_fullrank_solve_modifiedlogistic <- function(a, c, w_initial, tol = 1e-3, max_iter = 100) {
  if (missing(w_initial)) w_initial= log(a / (1-a))  #solution for the c=infinity case
  # Define f(w) and its derivative f'(w)
  f <- function(w) { 1 / (1 + exp(-w)) + w / c - a }
  f_prime <- function(w) { exp(-w) / (1 + exp(-w))^2 + 1 / c }
  # Initialize
  w <- w_initial
  iter <- 0
  converged <- FALSE
  # Newton-Raphson iterations
  while (iter < max_iter && !converged) {
    f_val <- f(w)
    f_deriv <- f_prime(w)
    if (abs(f_val) < tol) {
      converged <- TRUE
    } else {
      w <- w - f_val / f_deriv  # Newton update
      iter <- iter + 1
    }
  }
  return(w)
}



# Find an initial value for w solving the M-step equation Z^T m(w) = Z^T postprob, where m(z) = 1 / (1 + exp{- Z w})
#
# The exact solution is returned when the matrix U containing the unique rows in Z has ncol(Z) rows and is full-rank
# Otherwise, an approximated w is returned by clustering Z into ncol(Z) clusters and replacing Z by the cluster means (which, if they're full-rank, admit a closed-form solution
# Input
# - Z: p x q matrix recording q meta-covariates for p covariates
# - postprob: vector of length p with posterior probabilities for each covariate at current w
# - priorvar.w: the prior on w is N(0, priorvar.w * V) where V= (Z^T Z/nrow(Z))^{-1}
# Output: a list with the following components
# - w: value of w
# - fullrankcase: TRUE if Z satisfies the full-rank case conditions
# - Uinv: inverse of U if fullrankcase==TRUE, and NULL otherwise
eBayes_mstep_init= function(Z, postprob, priorvar.w=priorvar.w) {
    if (is.vector(postprob)) postprob= matrix(postprob, ncol=1)
    U= unique(Z)
    wini= NULL
    fullrankcase= FALSE
    if (nrow(U) == ncol(U)) {
        Uinv= try(solve(U), silent=TRUE)
        if (!inherits(Uinv, "try-error")) {
            wini= eBayes_mstep_fullrank(Z, Uinv=Uinv, postprob=postprob, priorvar.w=priorvar.w)
            fullrankcase= TRUE
        }
    }
    #if (is.null(wini)) {
    #    clus= kmeans(Z, centers=ncol(Z))$cluster
    #    clusmeans= aggregate(Z, by = list(clus = clus), FUN = mean)
    #    Zc= merge(data.frame(rowid=1:nrow(Z), clus), clusmeans, by="clus")
    #    Zc= Zc[order(Zc$rowid,decreasing=FALSE),]
    #    Zc= as.matrix(Zc[,-1:-2]) #remove cluster and row id
    #    Uinv= try(solve(clusmeans[,-1]), silent=TRUE)
    #    if (!inherits(Uinv, "try-error")) wini= eBayes_mstep_fullrank(Zc, Uinv=Uinv, postprob=postprob, priorvar.w=priorvar.w)
    #}
    if (is.null(wini)) {
        logitpp= log((postprob + .001) / (1 - postprob + .001))
        wini= coef(lm(logitpp ~ -1 + Z))
    }
    # Return output
    if (!fullrankcase) Uinv= NULL
    ans= list(w= matrix(wini, ncol=1), fullrankcase= fullrankcase, Uinv= Uinv)
    return(ans)
}


# Use Newton-Raphson to solve the M-step equation Z^T m(w) = Z^T postprob, where m(z) = 1 / (1 + exp{- Z w})
# Input
# - wcurrent: current value of w
# - Z: p x q matrix containing the q meta-covariates for the p covariates. It should not include an intercept, as it's already automatically added
# - postprob: vector of length p with posterior probabilities for each covariate at w= wcurrent
# - priorvar.w: prior on w is N(0, priorvar.w V) where V= (Z^T Z/nrow(Z))
# - maxiter: maximum number of iterations
# Output
# - w: solution of the M-step equation
# - grad: gradient of the objective function Z^T ( m(w) - postprob ) at w. Upon convergence, it should be zero
eBayes_mstep= function(wcurrent, Z, postprob, priorvar.w=priorvar.w, Vinv=Vinv, priorhess=priorhess, maxiter=10) {
    if (missing(priorhess)) priorhess= - Vinv / priorvar.w
    # Peform a first Newton-Raphson step
    w = wcurrent
    priorprob= 1 / (1 + exp(- Z %*% w))
    priorgrad= - (Vinv %*% w) / priorvar.w
    grad= priorgrad + eBayes_em_logit_grad(w, postprob=postprob, Z=Z, priorprob=priorprob)
    hess= priorhess + eBayes_em_logit_hess(w, postprob=postprob, Z=Z, priorprob=priorprob)
    delta= - base::solve(hess, grad)
    #Id= diag(ncol(Z))
    #delta= - base::solve(hess - Id, grad)
    w= w + delta
    maxstep= max(abs(delta))
    if (maxstep < 0.01) found= TRUE else found= FALSE
    iter= 0
    while (!found) {
        #eBayes_em_logit_objective(w, msfit=ms, Z=Z, priorvar.w=priorvar.w, Vinv=Vinv)
        priorprob= 1 / (1 + exp(- Z %*% w))
        priorgrad= - (Vinv %*% w) / priorvar.w
        grad= priorgrad + eBayes_em_logit_grad(w, postprob=postprob, Z=Z, priorprob=priorprob)                
        hess= priorhess + eBayes_em_logit_hess(w, postprob=postprob, Z=Z, priorprob=priorprob)
        #numDeriv::grad(eBayes_em_logit_objective, w, msfit=ms, Z=Z, priorvar.w=priorvar.w, Vinv=Vinv) #matches
        #numDeriv::hessian(eBayes_em_logit_objective, w, msfit=ms, Z=Z, priorvar.w=priorvar.w, Vinv=Vinv)  #matches
        #plot_taylor_approximation(f=eBayes_em_logit_objective, grad_f= eBayes_em_logit_grad, hessian_f= eBayes_em_logit_hess, x0= w, coord_idx=1, x_lim = c(-7,-4), n_points = 50, Z=Z, msfit=ms, postprob=pp0)
        #plot_taylor_approximation(f=eBayes_em_logit_objective, grad_f= eBayes_em_logit_grad, hessian_f= eBayes_em_logit_hess, x0= w, coord_idx=2, x_lim = c(2,5), n_points = 50, Z=Z, msfit=ms, postprob=pp0)
        #delta_damp= - base::solve(hess - Id, grad)
        delta= - base::solve(hess, grad)
        w= w + delta
        iter= iter + 1
        maxstep= max(abs(delta))
        if ((iter > maxiter) || maxstep < 0.001) found= TRUE
    }
    # Return output
    ans= list(w=w, grad=grad)
}


# Objective function in the empirical Bayes problem
# The objective function is sum_gamma pi(y | gamma) pi(gamma | w)= sum_gamma h(gamma) pi(gamma | y, w)
# where h(gamma)= pi(y | gamma) pi(gamma | w) / pi(gamma | y, w)
eBayes_logit_objective= function(w, msfit, Z, priorvar.w, Vinv) {
    ans= log(sum(exp(msfit$postProb)))
    ans= ans - (0.5 / priorvar.w) * (matrix(w,nrow=1) %*% Vinv %*% matrix(w,ncol=1))  #subtract log-prior
    return(ans)
}


# Objective function in the M-step under a logit link
# The objective function is f(w)= E[ log pi(gamma | w) | y, w=wcur ]
# - w: q-dimensional vector with hyper-parameter values
# - msfit: output of modelSelection when w=wcur
# - postprob: posterior inclusion probabilities pi(gamma_j = 1 | wcur), where wcur is the current hyper-parameter value
# - Z: p x q matrix containing the q meta-covariates for the p covariates. It should not include an intercept, as it's already automatically added
# - postprob, priorprob: ignored. Kept for compability with eBayes_em_logit_grad
# - priorvar.w, Vinv: prior on w is N(0, priorvar.w V), and Vinv= V^{-1}
# Output: f(w) at w
eBayes_em_logit_objective= function(w, msfit, Z, postprob, priorprob, priorvar.w, Vinv) {
    if (missing(priorprob)) priorprob= 1 / (1 + exp(- Z %*% w))
    ppmodel= postProb(msfit)
    logprior= double(nrow(ppmodel))
    for (i in 1:nrow(ppmodel)) {
        gamma= as.numeric(strsplit(ppmodel$modelid[i], split=",")[[1]])
        logprior[i]= sum(log(priorprob[gamma])) + sum(log((1-priorprob)[-gamma]))
    }
    ans= weighted.mean(logprior, w=ppmodel$pp)  #log marginal likelihood
    ans= ans - (0.5 / priorvar.w) * (matrix(w,nrow=1) %*% Vinv %*% matrix(w,ncol=1))  #subtract log-prior
    return(ans)
}


# Gradient of the objective function in the M-step under a logit link
# The objective function is f(w)= E[ log pi(gamma | w) | y, w=wcur ]
# - w: q-dimensional vector with hyper-parameter values
# - postprob: posterior inclusion probabilities pi(gamma_j = 1 | wcur), where wcur is the current hyper-parameter value
# - Z: p x q matrix containing the q meta-covariates for the p covariates. It should not include an intercept, as it's already automatically added
# - priorprob: prior inclusion probabilities 1/(1+exp(Z w)). If not specified, they're computed from (Z,w)
# - msfit: ignored (kept for compability with eBayes_em_logit_objective)
# Output: gradient of f(w) at w
eBayes_em_logit_grad= function(w, postprob, Z, priorprob, msfit) {
    if (missing(priorprob)) priorprob= 1 / (1 + exp(- Z %*% w))
    grad= t(Z) %*% (postprob - priorprob)
    return(grad)
}

# Same as eBayes_em_logit_grad, but returns the hessian of f(w)
eBayes_em_logit_hess= function(w, postprob, Z, priorprob, msfit) {
    if (missing(priorprob)) priorprob= 1 / (1 + exp(- Z %*% w))
    Zstd= sqrt(priorprob[,1]) * Z
    hess= - t(Zstd) %*% Zstd
    return(hess)
}


