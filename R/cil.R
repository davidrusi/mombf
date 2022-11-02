################################################################################
# List of functions
################################################################################
# * inside
# * lasso.bic
# * binary.id
# * compact.id
# * pinc.mt
# * grid.opt.mt
# * Of.EB.cil
# * Gf.EB.cil
# * Of.EP.cil
# * Gf.EP.cil
# * check.parameter.format
# * exposureCoef
# * model.pprobs.cil
# * bma.cil.teff
# * check.input.format
# * cil.teff
################################################################################
# List of methods
################################################################################
# * postProb
# * coef
################################################################################

# ##############################################################################
# # Dependencies
# library(mombf)
# library(pracma)
# library(glmnet)
# library(hdm)
# ##############################################################################

################################################################################
# Test if value "x" is inside interval ab = (a, b)
inside <- function(x, ab, include.ab = TRUE) {
################################################################################
  ifelse(include.ab == TRUE,
    (x >= ab[1] && x <= ab[2]), (x > ab[1] && x < ab[2]))
}



biclm <- function(y, pred, x, beta) {
  if (missing(pred)) pred <- x %*% beta
  n <- ifelse(is.vector(y), length(y), nrow(y))
  npar <- sum(beta != 0)
  ans <- n * log(colSums((y - pred)^2) / length(y)) + n * (log(2 * pi) + 1) + log(n) * npar
  return(ans)
}


biclogreg <- function(y, pred, x, beta) {
  if (missing(pred)) pred <- x %*% beta
  n <- ifelse(is.vector(y), length(y), nrow(y))
  npar <- sum(beta != 0)
  ans <- - 2 * sum(dbinom(y, size=1, prob=expit(-pred), log=TRUE)) + n * log(npar)
  return(ans)
}

################################################################################
# LASSO estimation with BIC
lasso.bic <- function(y, x, intercept = TRUE, standardize = TRUE, family = 'gaussian') {
################################################################################
  if (!is.matrix(x)) x <- as.matrix(x)
  fit <- glmnet(x = x, y = y, family = family, alpha = 1, intercept = intercept, standardize = standardize)
  if (intercept == TRUE) {
    ct <-  which(apply(x, 2, 'sd') == 0)
    if (length(ct)==0) x <- cbind(1,x)
    beta <- as.matrix(rbind(fit[['a0']], fit[['beta']]))
  } else {
    beta <- as.matrix(fit[['beta']])
  }
  
  n <- length(y)
  p <- colSums(beta != 0)

  if (family == 'gaussian') {
    bic <-  sapply(1:ncol(beta), function(i) biclm(y, x=x, beta=beta[,i,drop=FALSE]))
  } else if (family == 'binomial') {
    bic <- sapply(1:ncol(beta), function(i) biclogreg(y=y, x=x, beta=beta[,i,drop=FALSE]))
  } else {
     stop(paste("family=",family,"not implemented"))
  }
  sel <- which.min(bic)
  ypred <- x %*% beta[, sel, drop=FALSE]
  if (intercept) {
    beta <- beta[-1,sel] #if intercept was added by glmnet, exclude it from the estimated coef
  } else {
    beta <- beta[,sel]
  }
  ans <- list(coef = beta, ypred = ypred, lambda.opt = fit[['lambda']][sel], lambda = data.frame(lambda = fit[['lambda']], bic = bic, nvars = p))
  return(ans)
}

################################################################################
# Converting model indexes to binary strings
binary.id <- function(x, p, nt = 1, conc = FALSE) {
################################################################################
  out <- rep(0, p + nt)
  out[as.numeric(unlist(strsplit(x, ',')))] <- 1
  if (conc == TRUE) { out <- paste(out, collapse = '') }
  return(out)
}

################################################################################
# Converting binary strings to model indexes
compact.id <- function(x) {
################################################################################
  paste(which(unlist(strsplit(as.character(x), '')) == 1), collapse = ',')
}

################################################################################
# Prior inclusion probability function (conditional on coefficient size)
pinc.mt <- function(x, th, rho.min = 0, rho.max = 1, squared = FALSE, includeX) {
################################################################################
  #  If there are interactions (CAREFUL: main effects + ORDERED interactions)
  th <- as.numeric(th)
  if (ncol(x) / (length(th) - 1) > 1) {
    th <- c(th, rep(th[-1], each = ncol(x) / (length(th) - 1) - 1))
  }

  # If using abs(beta) or beta^2
  if (squared == TRUE) {
    probs <- (1 + exp(-apply(x, 1, function(z) {
        sum(c(1, z^2) * th)
      })))^(-1)
  } else {
    probs <- (1 + exp(-apply(x, 1, function(z) {
        sum(c(1, abs(z)) * th)
      })))^(-1)
  }

  # Limit upper and lower bounds
  probs <- rho.min + (rho.max - rho.min) * probs

  # Force inclusion of variables indicated by includeX
  probs[includeX] <- 1
    
  # End
  return(probs)
}

################################################################################
grid.opt.mt <- function(G0, betad, pj1, method = 'EB', ws = NA, th.grid,
  th.prior, rho.min = 0, rho.max = 1, ret.plotly = FALSE, includeX, Dvars) {
################################################################################
  opt.th <- th.grid
  if (th.prior == 'tunif') {
    opt.th <- opt.th[which(opt.th[, 1] >= -log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 2] >= -log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 1] <= log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 2] <= log(length(betad))), ]
  }

  #opt.th <- expand.grid(th0s, th1s)
  if (method == 'EB') {
    # Evaluate -log of objective function
    fval <- exp(-apply(opt.th, 1, Of.EB.cil, betad = betad, G0 = G0, ws = ws, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, includeX = includeX, Dvars = Dvars))
  } else if (method == 'EP') {
    fval <- exp(-apply(opt.th, 1, Of.EP.cil, betad = betad, pj1 = pj1, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, includeX = includeX))
  } else {
    stop('supported theta search methods are "EB" and "EP".')
  }

  # Optimal values. If there are multiple modes, choose the one with smallest L2 norm
  th.opt <- opt.th[fval == max(fval),,drop=FALSE]
  th.opt <- as.numeric(th.opt[which.min(th.opt[,1]^2 + th.opt[,2]^2),])
  #th.opt <- as.numeric(opt.th[which.max(fval),])
  names(th.opt) <- paste('th', 1:length(th.opt) - 1, sep = '')

  # Output
  return(list(th.opt = th.opt))
}#; grid.opt.mt <- compiler::cmpfun(grid.opt.mt)

################################################################################
# OBJECTIVE FUNCTIONS FOR EB METHOD
# - x: theta hyper-parameter value
# - betad: measures of association between treatment d and controls
# - G0: model identifiers. Each element is a character string where 1's indicate inclusion and 0's exclusion
# - ws: marginal likelihood for each model in G0
# - th.prior: prior on theta, currently not used
# - rho.min, rho.max: the prior inclusion prob given by any theta is restricted to [rho.min, rho.max], to avoid prior inclusion prob=0 or 1
# - includeX: variables whose inclusion is forced into the model
# - Dvars: index of variables corresponding to treatment. These have prior inclusion prob = 1/2, regardless of theta
Of.EB.cil <- function(x, betad, G0, ws, th.prior, rho.min = 0, rho.max = 1, includeX, Dvars) {
################################################################################
  # log p(y | theta)
  p1 <- rep(NA, nchar(G0)[1])
  p1[-Dvars] <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max, includeX = includeX)
  p1[Dvars] <- 1/2
  foo <- function(m) {
    gm <- as.numeric(unlist(strsplit(m, '')))
    pg <- prod((p1 ** gm) * (1 - p1) ** (1 - gm))  #prior model probability
    wg <- ws[which(G0 == m)] #marginal likelihood of the model
    return(wg * pg)
  }
  sc <- sum(sapply(G0, foo))
  pt1 <- log(sum(sc))  # log p(y | theta)

  # log p(theta)
  pt2 <- 0

  return(-(pt1 + pt2))
}

################################################################################
Gf.EB.cil <- function(x, betad, G0, ws, th.prior, rho.min = 0, rho.max = 1, includeX, Dvars) {
################################################################################
  # Compute (unnormalised) model inclusion probabilities
  p1 <- rep(NA, nchar(G0)[1])
  p1[-Dvars] <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max, includeX = includeX)
  p1[Dvars] <- 1/2
  #sD <- nchar(G0[1]) - length(p1)
  #p1 <- c(rep(1/2, sD), p1)
  sc <- sapply(G0, function(m) {
    gm <- as.numeric(unlist(strsplit(m, '')))
    wg <- ws[which(G0 == m)]
    pg <- (p1 ** gm) * (1 - p1) ** (1 - gm)
    return(wg * prod(pg))
  })

  # Marginal inclusion probabilities
  pprobs <- data.frame('model_id' = G0, 'pprob' = sc / sum(sc))
  pprobs[, 1] <- sapply(pprobs[, 1], compact.id)
  mpprobs <- apply(t(sapply(pprobs[, 1], binary.id,
    p = nchar(G0[1]) - ncol(betad), nt = ncol(betad))), 2,
    function(x) { sum(x * pprobs[, ncol(pprobs)]) })

  # grad log p(y | theta)
  gs <- rep(NA, length(x))
  gs[1] <- sum((mpprobs - p1)[-Dvars])
  #gs[1] <- sum((mpprobs - p1)[-(1:sD)])
  for (k in 2:length(gs)) {
    aux0 <- (ncol(betad) / (length(x) - 1) - 1)
    idxs <- c(k - 1)
    if (aux0 > 0) { idxs <- c(idxs, (length(x) - 1) + aux0 * (k - 2) + 1:aux0) }
    aux1 <- mpprobs - p1
    aux2 <- rowSums(abs(as.matrix(betad[, idxs])))
    gs[k] <- sum(aux1[-Dvars] * aux2)
    #gs[k] <- sum(aux1[-(1:sD)] * aux2)
  }

  # grad log p(theta)
  t0 <- t1 <- 0

  # - grad log p(theta | y) = -(grad log p(y | theta) + grad log p(theta))
  return(-(gs + t0))
}

################################################################################
# OBJECTIVE FUNCTION FOR EP METHOD
# - x: theta hyper-parameter value
# - betad: measures of association between treatment d and controls
# - pj1: marginal posterior inclusion probabilities in outcome model under uniform model prior
# - th.prior: prior on theta, currently not used
# - rho.min, rho.max: the prior inclusion prob given by any theta is restricted to [rho.min, rho.max], to avoid prior inclusion prob=0 or 1
# - includeX: variables whose inclusion is forced into the model
Of.EP.cil <- function(x, betad, pj1, th.prior, rho.min = 0, rho.max = 1, includeX) {
################################################################################
  # log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max, includeX = includeX)
  #David: the line below seems wrong: the EP objective fun is a sum involving control prob in p1, excluding the treatment
  #p1 <- c(rep(1/2, length(pj1) - length(p1)), p1)
  #David: the corrected line is below (make pj1 shorter, rather than make p1 longer)
  sD <- length(pj1) - length(p1); pj1 <- pj1[-(1:sD)]
  fj <- pj1 * p1 + (1 - pj1) * (1 - p1)
  pt1 <- sum(log(fj))

    # log p(theta)
  pt2 <- 0

  # - log p(theta | y) \propto -(log p(y | theta) + log p(theta))
  return(-(pt1 + pt2))
}

################################################################################
Gf.EP.cil <- function(x, betad, nt, pj1, th.prior, rho.min = 0, rho.max = 1, includeX) {
################################################################################
  # grad log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max, includeX= includeX)
  sD <- length(pj1) - length(p1)
  p1 <- c(rep(1/2, sD), p1)
  fj <- pj1 * p1 + (1 - pj1) * (1 - p1)

  gs <- rep(NA, length(x))
  gs[1] <- sum(((((1 - p1) * p1) / fj) * (2 * pj1 - 1))[-(1:sD)])
  for (k in 2:length(gs)) {
    aux0 <- (ncol(betad) / (length(x) - 1) - 1)
    idxs <- c(k - 1)
    if (aux0 > 0) { idxs <- c(idxs, (length(x) - 1) + aux0 * (k - 2) + 1:aux0) }
    aux1 <- (((1 - p1) * p1) / fj) * (2 * pj1 - 1)
    aux2 <- rowSums(abs(as.matrix(betad[, idxs])))
    gs[k] <- sum(aux1[-(1:sD)] * aux2)
  }

  # grad log p(theta)
  t0 <- 0

  # - grad log p(theta | y) = -(grad log p(y | theta) + grad log p(theta))
  return(-(gs + t0))
}

################################################################################
check.parameter.format <- function(rho.min, th.range, max.mod, lpen, eps, bvs.fit0, th.EP, D) {
################################################################################
  # Error messages accounted for
  err1 <- 'if informed, argument "rho.min" must be a scalar in the interval (0, 1/2).'
  #err2 <- 'if informed, argument "rho.max" must be a scalar in the interval (1/2, 1).'
  err3 <- 'if informed, argument "th.range" must be numeric and have at least two distinct values.'
  #err4 <- 'if informed, argument "tau" must be a positive scalar.'
  err5 <- 'parameter "max.mod" must be a positive integer scalar.'
  err6 <- 'parameter "lpen" must be set to "lambda.min" or "lambda.1se" (see documentation for "glmnet").'
  err7 <- 'argument "eps" must be a scalar in the interval (0, 1/2).'
  err8 <- 'if informed, parameter "bvs.fit0" must be of class "msfit".'
  err9 <- 'if informed, parameter "th.EP" must be a numeric vector of length equal to the number of columns in "D" plus 1.'

  # Parameter: "rho.min"
  if (! is.null(rho.min) | length(rho.min) > 1) {
    if (! is.numeric(rho.min) | length(rho.min) > 1) {
      stop(err1)
    } else if (! inside(rho.min, c(0, 1/2), include.ab = FALSE)) {
      stop(err1)
    }
  }

  # Parameter: "th.range"
  if (! is.null(th.range)) {
    if (! is.numeric(th.range)) {
      stop(err3)
    } else if (! length(unique(th.range)) == length(th.range)) {
      stop(err3)
    }
  }

  # Parameter: "tau"
  #if (! is.null(tau) | length(tau) > 1) {
  #  if (! is.numeric(tau) | length(tau) > 1) {
  #    stop(err4)
  #  } else if (sign(tau) != 1) {
  #    stop(err4)
  #  }
  #}

  # Parameter: "max.mod"
  if (! is.null(max.mod) | length(max.mod) > 1) {
    if (! is.numeric(max.mod) | length(max.mod) > 1) {
      stop(err5)
    } else if (sign(max.mod) != 1) {
      stop(err5)
    }
  }

  # Parameter: "lpen"
  if (length(lpen) != 1) {
    stop(err6)
  } else if (! lpen %in% c('lambda.min', 'lambda.1se')) {
    stop(err6)
  }

  # Parameter: "eps"
  if (! is.null(eps) | length(eps) > 1) {
    if (! is.numeric(eps) | length(eps) > 1) {
      stop(err7)
    } else if (! inside(eps, c(0, 1/2), include.ab = FALSE)) {
      stop(err7)
    }
  }

  # Object: "bvs.fit0"
  if (! is.null(bvs.fit0)) {
    if (!inherits(bvs.fit0, 'msfit')) {
      stop(err8)
    }
  }

  # Object: "th.EP"
  if (! is.null(th.EP)) {
    if (! is.numeric(th.EP) | length(th.EP) != ncol(D) + 1) {
      stop(err9)
    }
  }
}



################################################################################
# Estimate coefficients regressing exposure (or treatments) on controls
# - Di: exposure (or treatment) variables
# - A: controls
# - familyD: type of outcome for each column in Di, e.g. c('normal','binomial')
# - typeofvar: does each column in A correspond to a 'numeric' or 'factor'. The latter are not standardized to zero mean, unit variance
# - addintcpt: should an intercept be added to A? (TRUE/FALSE)
# - mod1: method to estimate the regression model
# - priorCoef, Rinit: only used if mod1=='bvs', 'bma' or 'bms'. Then it's the prior on the coefficients and number of MCMC iterations used by modelSelection to estimate the coefficients
# - lpen: only used if mod1=='lasso'. It's an option when choosing lambda in LASSO with cross-validation, 'lambda.min' or 'lambda.1se'
exposureCoef <- function(Di, A, familyD, typeofvar = rep('numeric',ncol(A)), addintcpt, mod1, priorCoef, Rinit, lpen) {
  n <- nrow(A); p <- ncol(A)
  if (length(familyD)==1) familyD <- rep(familyD, ncol(Di))
  familyD.glmnet <- familyD
  familyD.glmnet[familyD.glmnet == 'normal'] <- 'gaussian'

  #Scale design matrix
  mx <- colMeans(A); sx <- sqrt(colMeans(A^2) - mx^2) * sqrt(n/(n-1))
  isct <- (sx==0)
  mx[typeofvar=='factor'] <- 0; sx[typeofvar=='factor'] <- 1
  A[,!isct]= t((t(A[,!isct]) - mx[!isct])/sx[!isct])

  #Obtain coefficients for each exposure (treatment) variable
  betad <- matrix(NA, nrow = ncol(A), ncol = ncol(Di))
  for (i in 1:ncol(Di)) {
    # Define matrices
    b <- Di[, i]
    exc <- 2 * p  # phantom column: not exclude anything
    if (length(unique(Di[, i])) <= 3 &  # This controls sum-to-zero constraint
        all(unique(Di[, i]) %in% c(-1, 0, 1))) {
      b[which(b == -1)] <- 0
    }

    # Exposure family
    ntreatvals <- length(unique(b))
    if (ntreatvals > 1) {
    #if (ntreatvals > 1 & min(table(b)) > 1) {  #David: this line seemed incorrect, for continuous b we had min(table(d))==1 and hence betad wasn't initialized
      if (mod1 == 'ginv') {
        if (familyD[i] != 'normal') stop("The generalized inverse method cannot be applied when familyD != 'normal'")
        betad[, i] <- as.matrix(pinv(t(A) %*% A) %*% t(A) %*% b)[-exc] #pracma::pinv
      } else if (mod1 %in% c('bvs', 'bma', 'bms')) {
        m1.fit <- modelSelection(b, A, family=familyD[i], priorCoef = priorCoef, priorVar = igprior(0.01,0.01), priorDelta = modelbbprior(1,1),
          niter = Rinit, verbose = FALSE)
        betad[,i] <- coef(m1.fit)[1:nrow(betad),1]
        #post.th <- colMeans(rnlp(msfit = m1.fit, priorCoef = priorCoef, niter = 1e4))
        #betad[, i] <- unname(post.th[2:(length(post.th) - 1)])[-exc]
      } else if (mod1 == 'lasso') {
        m1.fit <- cv.glmnet(x=A, y=b, intercept = addintcpt, family = familyD.glmnet[i])
        m1.coef <- as.numeric(coef(m1.fit, s = m1.fit[[lpen]])[-exc])
        if (nrow(betad) == length(m1.coef)) betad[,i] <- m1.coef else betad[,i] <- m1.coef[-1]
        #betad[, i] <- unname(coef(m1.fit, s = m1.fit[[lpen]])[-1, 1])[-exc]
      } else if (mod1 == 'lasso_bic') {
        m1.fit <- lasso.bic(y=b, x=A, intercept = addintcpt, family = familyD.glmnet[i])
        betad[, i] <- unname(m1.fit[['coef']][-exc])
        #if (addintcpt == FALSE) { m1.fit[['coef']][2] <- m1.fit[['coef']][1] }
        #betad[, i] <- unname(m1.fit[['coef']][-1][-exc])
      } else if (mod1 == 'ridge') {
        m1.fit <- cv.glmnet(A, b, alpha = 0, intercept = addintcpt, family = familyD.glmnet[i])
        betad[, i] <- unname(coef(m1.fit, s = m1.fit[['lambda.min']]))[-exc]
        #betad[, i] <- unname(coef(m1.fit, s = m1.fit[['lambda.min']])[-1, ])[-exc]
      } else {
        stop('currently unsupported "mod1" method.')
      }
    } else {
      betad[, i] <- rep(0, nrow(betad))  # If treatment has no variability
    }
  }
  colnames(betad) <- colnames(Di)
  rownames(betad) <- colnames(A)
  return(betad)
}
################################################################################


################################################################################
#Set range of values to consider for hyper-parameter theta in the EB/EP optimization
setThetaRange <- function(ncolD) {
  if (ncolD == 1) {
    th.range <- seq(-40, 40, 1)
  } else if (ncolD == 2) {
    th.range <- seq(-40, 40, 2)
  } else if (ncolD <= 4) {
    th.range <- seq(-10, 10, 2)
  } else if (ncolD <= 7) {
    th.range <- seq(-7.5, 7.5, 2.5)
  } else if (ncolD <= 14) {
    th.range <- seq(-2, 2, 2)
    #th.range <- seq(-5, 5, 2.5)
  } else {
    th.range <- c(-2, 2)
  }
  return(th.range)
}
################################################################################


################################################################################
model.pprobs.cil <- function(y, D, X, I = NULL, family = 'normal', familyD = 'normal',
  mod1, th.search = 'EB', th.prior = 'unif', priorCoef, rho.min = NULL,
  th.range = NULL, max.mod = 2^20, lpen = 'lambda.1se', eps = 1e-10,
  R, Rinit= 500, bvs.fit0 = NULL, th.EP = NULL, center = center, scale = scale,
  includevars = includevars, verbose = TRUE) {
################################################################################
  # Make sure parameter inputs are in the correct format
  check.parameter.format(rho.min = rho.min, th.range = th.range, max.mod = max.mod, lpen = lpen, eps = eps, bvs.fit0 = bvs.fit0, th.EP = th.EP, D = D)

  # Renaming (DEPENDENCIES: mombf, pracma, glmnet, hdm)
  ms <- modelSelection#mombf::modelSelection
  igp <- igprior(0.01, 0.01)#mombf::igprior(0.01, 0.01)
  bbp <- modelbbprior(1, 1)#mombf::modelbbprior(1, 1)
  ufp <- modelunifprior()#mombf::modelunifprior()
  mlik <- nlpMarginal#mombf::nlpMarginal

  # These parameters must be integers 
  max.mod <- round(max.mod, 0)
  R <- round(R, 0)

  # Create design matrices
  if (!is.null(I)) {
    Z <- cbind(D, I, X)
    Di <- cbind(D, I)
    ncolI <- ncol(I)
  } else {
    Z <- cbind(D, X)
    Di <- D
    ncolI <- 0
  }

  if (is.data.frame(X)) {
    f <- as.formula(paste('~', paste(names(X), collapse=' + ')))
    A <- createDesign(formula=f, data=X)
    typeofvar <- A$typeofvar
    A <- A$x
  } else if (is.matrix(X)) {
      A <- X
      typeofvar <- rep('numeric',ncol(A))
  } else {
    stop("X must be a matrix or a data.frame")
  }

  if (missing(includevars)) includevars <- rep(FALSE, ncol(D)+ncolI+ncol(A))
  includeX <-  includevars[(ncol(D)+ncolI+1):length(includevars)]
  if (length(includevars) != ncol(D) + ncolI + ncol(A)) stop(paste("includevars has length",length(includevars),"but there are",ncol(D)+ncolI+ncol(A),"columns in (D,I,X). Note: if X is a data.frame, an intercept was automatically added"))

  isct <- (apply(A, 2, 'sd') == 0)
  if (sum(isct) > 1) stop("There are >1 constant columns (e.g. intercepts) in X. Try removing the intercept")
  addintcpt <- sum(isct)==0
  includeX[isct] <- TRUE #always add the intercept (so it doesn't affect EB/EP estimates)
   
  # Other fixed parameters
  if (is.null(rho.min)) rho.min <- 1 / (ncol(Z) + 1)
  #if (is.null(rho.min)) rho.min <- 1 / (ncol(Z)^2 + 1)
  rho.max <- 1 - rho.min

  # Estimate coefficients on exposure model
  if (verbose) cat("ESTIMATING ASSOCIATION BETWEEN TREATMENTS AND CONTROLS\n\n")
  betad <- exposureCoef(Di=Di, A=A, familyD=familyD, typeofvar=typeofvar, addintcpt=addintcpt, mod1=mod1, priorCoef=priorCoef, Rinit=Rinit, lpen=lpen)
    
  # THETA SEARCH ###############################################################
  # Initial fit: th0 = 0; th1 = 0 (UNIFORM MODEL PRIOR)
  if (is.null(bvs.fit0)) {
    if (verbose) cat("INITIAL FIT UNDER UNIFORM MODEL PRIOR \n\n")
    if (is.data.frame(Z)) {
      f <- as.formula(paste('y ~', paste(names(Z), collapse=' + ')))
      data <- cbind(y, Z)
      Zdes <- formatInputdata(f, data=data, family=family)$x
      p <- ncol(Zdes) - ncol(Di)
      isct <- apply(Zdes, 2, var) == 0
      includevars[isct] <- TRUE #always include the intercept (so it doesn't affect EB/EP estimates)
      Dvars <- 2:(ncol(D)+1) #columns in Z corresponding to D (1st column is the intercept)
      bvs.fit0 <- ms(f, data=data, priorCoef = priorCoef, priorVar = igp, priorDelta = ufp, niter = Rinit, verbose = verbose, center = center, scale = scale, includevars = includevars)
    } else {
        p <- ncol(Z) - ncol(Di)
        isct <- apply(Z, 2, var) == 0
        includevars[isct] <- TRUE #always include the intercept (so it doesn't affect EB/EP estimates)
        Dvars <- 1:ncol(D)  #columns in Z corresponding to D
        bvs.fit0 <- ms(y, Z, priorCoef = priorCoef, priorVar = igp, priorDelta = ufp, niter = Rinit, verbose = verbose, center = center, scale = scale, includevars = includevars)
    }
  }
  # Posterior model probabilities
  pprobs <- postProb(bvs.fit0)[, c(1, 3)]

  # Marginal inclusion probabilities
  pj1 <- bvs.fit0[['margpp']]
  pj1[which(is.nan(pj1))] <- 0  # Potential problems in mombf with undefined
  pj1[which(pj1 == 1 & !includevars)] <- 1 - eps  # In case there are extreme values
  pj1[which(pj1 == 0)] <- eps

  # Set of visited models
  max.mod <- min(max.mod, nrow(pprobs))
  G0 <- unname(sapply(as.character(pprobs[1:max.mod, 1]), binary.id, p = p, nt = ncol(Di), conc = TRUE)) # Prevent a dictionary too big

  if (verbose) cat("\n ESTIMATING HYPER-PARAMETER VALUES \n\n")
  # Determine the values of hat{theta} #########################################
  # Set limits for optimisation (in principle: unbounded)
  lws <- rep(-Inf, ncol(D) + 1)
  ups <- rep(+Inf, ncol(D) + 1)

  # Range for grid search
  if (is.null(th.range)) th.range <- setThetaRange(ncol(D))

  if (is.null(th.EP)) {
    # EP approximation (to set starting point)
    ths <- vector(ncol(D) + 1, mode = 'list')
    for (i in 1:length(ths)) { ths[[i]] <- th.range }
    EP.is <- grid.opt.mt(G0, betad, pj1, th.grid = expand.grid(ths), method = 'EP', th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, includeX = includeX, Dvars = Dvars)

    # Load gradient functions and optimise
    st <- unlist(EP.is[[1]])
    opt.EP <- nlminb(st, objective = Of.EP.cil, gradient = Gf.EP.cil, 
      lower = lws, upper = ups, pj1 = pj1, betad = betad, 
      th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, includeX = includeX)

    # Set values if there is convergence
    if (opt.EP[['convergence']] %in% 0:1) {
      th.EP <- opt.EP[['par']]
    } else {
      stop('theta search under "EP" did not converge.')
    }
  }

  # Empirical Bayes search (starting at EP solution)
  if (th.search == 'EB') {
    # Model search at EP optimum
    auxprpr <- pinc.mt(betad, th = th.EP, rho.min = rho.min, rho.max = rho.max, includeX = includeX)
    auxprpr <- c(rep(1/2, ncol(betad)), auxprpr)
    auxprpr[which(auxprpr == 1)] <- 1 - eps
    auxprpr[which(auxprpr == 0)] <- eps
    auxmp <- modelbinomprior(p = auxprpr)#mombf::modelbinomprior(p = auxprpr)

    if (is.data.frame(Z)) {
      update.fit <- ms(f, data=data, priorCoef=priorCoef, priorVar=igp, priorDelta=auxmp, niter=round(R/2), verbose=verbose, center=center, scale=scale, includevars=includevars)
    } else {
      update.fit <- ms(y, Z, priorCoef = priorCoef, priorVar = igp, priorDelta = auxmp, niter = round(R/2), verbose = verbose, center = center, scale = scale, includevars = includevars)
    }

    # Append model set
    upprobs <- postProb(update.fit)[, c(1, 3)]
    G1 <- unname(sapply(
      as.character(upprobs[, 1]), binary.id, p = p, nt = ncol(D), conc = TRUE))
    G0 <- unique(c(G0, G1[1:min(max.mod, length(G1))]))

    # Pre-compute marginal likelihoods
    ws <- unname(unlist(sapply(G0, function(m) {
      js <- which(unlist(strsplit(m, '')) == 1)
      return(mlik(js, y, update.fit$xstd, priorCoef = priorCoef, logscale = TRUE)) #David: line below seemed wrong, ws was exponentiated later
      #return(mlik(js, y, Z, priorCoef = priorCoef, logscale = FALSE))
    })))

    # Avoid numerical problems with MLs that are too small
    ws <- exp(ws - max(ws)) # Multiply by constant exp(-max(ws)) (NO EFFECT ON ARG MAX)

    # Optimise starting at the EP optimum
    opt.EB <- try(nlminb(th.EP, objective = Of.EB.cil, gradient = Gf.EB.cil,
      lower = lws, upper = ups, betad = betad, G0 = G0,
      ws = ws, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, includeX = includeX, Dvars = Dvars), silent = TRUE)
    if (inherits(opt.EB, 'try-error')) {
      opt.EB <- try(nlminb(th.EP, objective = Of.EB.cil,
        gradient = Gf.EB.cil, lower = lws, upper = ups, betad = betad,
        G0 = G0, ws = ws, th.prior = th.prior, rho.min = rho.min,
        rho.max = rho.max, includeX = includeX, Dvars = Dvars))
    }

    # Optimal values
    if (opt.EB[['convergence']] %in% 0:1) {
      th.hat <- opt.EB[['par']]
    } else {
      stop('theta search under th.search=="EB" did not converge.')
    }
  } else if (th.search == 'EP') {
    # Optimal values
    th.hat <- th.EP
  } else {
    stop('th.search method unsupported -- choose amongst "EB" and "EP".')
  }

  # MODEL SEARCH (under fixed theta hat) #######################################
  # Prior probabilities under theta
  pg1cth <- rep(NA, nchar(G0)[1])
  pg1cth[-Dvars] <- pinc.mt(betad, th = th.hat, rho.min = rho.min, rho.max = rho.max, includeX = includeX)
  pg1cth[Dvars] <- 1/2
  #pg1cth[which(pg1cth == 1 & !includevars)] <- 1 - eps

  # Bernoulli model prior with unequal probabilities
  nbp <- modelbinomprior(p = pg1cth)#mombf::modelbinomprior(p = pg1cth)
 
  # New model search
  #new.fit <- ms(y, Z, priorCoef = priorCoef, priorVar = igp, priorDelta = nbp, niter = R, verbose = verbose, center = center, scale = scale, includevars = includevars)
  if (verbose) cat("\n FINAL FIT UNDER ESTIMATED ASYMMETRIC BERNOULLI PRIOR \n\n")
  if (is.data.frame(Z)) {
    new.fit <- ms(f, data=data, priorCoef = priorCoef, priorVar = igp, priorDelta = nbp, niter = R, verbose = verbose, center = center, scale = scale, includevars = includevars)
    treatidx <- Dvars
  } else {
    new.fit <- ms(y, Z, priorCoef = priorCoef, priorVar = igp, priorDelta = nbp, niter = R, verbose = verbose, center = center, scale = scale, includevars = includevars)
    treatidx <- Dvars + ifelse(all(apply(Z,2,'sd') > 0), 1, 0) #if Z has no intercept, coef(new.fit) will add one automatically
  }

  # BMA inference and Posterior model probabilities
  beta <- coef(new.fit)
  teff <- beta[treatidx,]
  pprobs <- postProb(new.fit)#mombf::postProb(new.fit)[, c(1, 3)]
  mpprobs <- new.fit[['margpp']]

  # Finish
  return(list(teff = teff, beta = beta, pprob.y = pprobs, mpprob.y = mpprobs, th.hat = th.hat, priorprobs = pg1cth, 
              msfit = new.fit, init.msfit = bvs.fit0, th.EP = th.EP, mod1.coef = betad, Dvars=Dvars, includeX=includeX, rho.min=rho.min, rho.max=rho.max))
}#; model.pprobs.cil <- compiler::cmpfun(model.pprobs.cil)


################################################################################
check.input.format <- function(y, D, X, I, R) {
################################################################################
  # Format errors that can be encountered
  err1 <- 'argument "y" must be of class "matrix" and contain one column.'
  err2 <- 'argument "D" must have at least one column and the same number of rows as "y".'
  #err2 <- 'argument "D" must be of class "matrix" with at least one column and the same number of rows as "y".'
  err3 <- 'argument "X" must have at least one column and the same number of rows as "y".'
  #err3 <- 'argument "X" must be of class "matrix" with at least one column and the same number of rows as "y".'
  err4 <- 'argument "I" must have at least one column and the same number of rows as "y".'
  #err4 <- 'if argument "I" is informed, it must be of class "matrix" with at least one column and the same number of rows as "y".'
  err5 <- 'parameter "R" must be a single positive integer.'

  # Object "y"
  if (! 'matrix' %in% class(y)) {
    stop(err1)
  } else if (ncol(y) != 1) {
    stop(err1)
  }

  # Object "D"
  if (nrow(D) != nrow(y) | ncol(D) == 0) stop(err2)

  # Object "X"
  if (nrow(X) != nrow(y) | ncol(X) == 0) stop(err3)

  # Object "I"
  if (! is.null(I)) {
    if (nrow(I) != nrow(y) | ncol(I) == 0) stop(err4)
    ncolI= ncol(I)
  } else {
      ncolI= 0
  }
  
  # Parameter "R"
  if (! is.numeric(R) | length(R) != 1) {
    stop(err5)
  } else if (R < 1) {
    stop(err5)
  }

  # NA warnings
  if (any(is.na(y))) { stop('"y" cannot contain NAs.') }
  if (any(is.na(D))) { stop('"D" cannot contain NAs.') }
  if (any(is.na(X))) { stop('"X" cannot contain NAs.') }
  if (! is.null(I) & any(is.na(I))) {
    stop('if informed, "I" cannot contain NAs.')
  }

  # End
  return(format.is.ok = TRUE)
}

################################################################################
cil <- function(y, D, X, I = NULL, family = 'normal', familyD = 'normal', R = 1e4, Rinit = 500, th.search = 'EB',
  mod1 = 'lasso_bic', th.prior = 'unif', priorCoef, rho.min = NULL, th.range = NULL, max.mod = 2^20,
  lpen = 'lambda.1se', eps = 1e-10, bvs.fit0 = NULL, th.EP = NULL, center = TRUE, scale = TRUE, includevars, verbose=TRUE) {
################################################################################
  # Assert inputs are in the correct format
  check.input.format(y, D, X, I, R)

  # Posterior model probabilities
  pprobs <- model.pprobs.cil(y, D, X, I, family = family, familyD = familyD, R = R, Rinit = Rinit, th.search = th.search,
    mod1 = mod1, th.prior = th.prior, priorCoef = priorCoef,
    rho.min = rho.min, th.range = th.range,
    max.mod = max.mod, lpen = lpen, eps = eps, bvs.fit0 = bvs.fit0,
    th.EP = th.EP, center = center, scale = scale, includevars = includevars, verbose = verbose)

  # BMA computation
  #teff <- bma.cil(pprobs[['msfit']], nt = ncol(cbind(D, I)))

  # Output object
  out <- list(cil.teff = pprobs[['teff']],  # BMA estimates, 95\% CIs and marginal posterior prob for treatment variables
              coef = pprobs[['beta']], # BMA estimates, 95\% CI's and marginal posterior probs for all variables
#cil.teff = teff[['bma.e.cil']],  # BMA Point estimates
#              cil.teff.postdist = teff[['bma.d.cil']],  # BMA Distribution
#              cil.bma.mcmc = teff[['mcmc.cil']],  # BMA MCMC runs
              model.postprobs = pprobs[['pprob.y']],  # Model post. probs.
              margpp = pprobs[['mpprob.y']],  # Marginal post. probs under CIL prior
              margprior = pprobs[['priorprobs']], # CIL-estimated marginal prior inclusion probabilities
              margprior.limits = c(pprobs[['rho.min']], pprobs[['rho.max']]), # (lower, upper) limits on margprior
              theta.hat = pprobs[['th.hat']],  # Estimated theta hat
              treat.coefs = pprobs[['mod1.coef']],  # Coefs. treatment models
              treat.varindex = pprobs[['Dvars']], # Indexes of columns in full design matrix corresponding to treatments
              msfit = pprobs[['msfit']],  # Model selection fit of the model
              theta.EP = pprobs[['th.EP']],  # Estimated theta for EP method
              includeX = pprobs[['includeX']], # Covariates forced to be included in the outcome model
              init.msfit = pprobs[['init.msfit']])  # Initial model sel. fit

  # Exit
  return(invisible(new('cilfit', out)))
  #return(invisible(out))
}



################################################################################
# DEPRECATED FUNCTION
#bma.cil <- function(msfit, nt = 1, ret.mcmc = TRUE) {
#################################################################################
#  # Posterior draws
#  draws <- rnlp(msfit = msfit, niter = 1e4)
# 
#  # BMA Treatment Effect Estimates
#  idxs <- 1 + 1:nt
#  if (nt == 1) {
#    teff.est <- mean(draws[, idxs])
#  } else {
#    teff.est <- colMeans(draws[, idxs])
#  }
# 
#  # BMA distribution
#  if (nt == 1) {
#    bma.distr <- quantile(draws[, idxs], seq(0, 1, 0.005))
#  } else {
#    bma.distr <- apply(draws[, idxs], 2, quantile, seq(0, 1, 0.005))
#  }
# 
#  # Output
#  out <- list(bma.e.cil = teff.est, bma.d.cil = bma.distr, mcmc.cil = NA)
#  if (ret.mcmc == TRUE) { out[['mcmc.cil']] <- draws[, 1 + 1:nt] }
# 
#  # End
#  return(invisible(out))
#}



################################################################################
# Adapt relevant methods for objects of class "cilfit"
################################################################################

# plot()
setMethod("plotprior", signature(object='cilfit'), function(object, xlab, ylab, ylim=c(0,1), ...) {
  if (missing(xlab)) xlab <- 'Treatment regression coefficients'
  if (missing(ylab)) ylab <- 'Posterior inclusion probability under uniform model prior'
  ntreat <- ncol(object$treat.coefs)
  treat.varindex <- object$treat.varindex
  bmean <- colMeans(object$treat.coefs)
  b <- t(matrix(bmean, nrow=ntreat, ncol=200))
  for (i in 1:ntreat) {
    devAskNewPage(ask = TRUE)
    pp.unifprior <- object$init.msfit$margpp[-treat.varindex]
    pp.unifprior <- pp.unifprior[!object$includeX]
    treat.coef <- object$treat.coef[!object$includeX,i]
    isct <- (apply(object$init.msfit$xstd[,,drop=FALSE], 2, 'sd') == 0)[-treat.varindex]
    isct <- isct[!object$includeX]
    pp.unifprior <- pp.unifprior[!isct]
    treat.coef <- treat.coef[!isct]
    plot(treat.coef, pp.unifprior, xlab=xlab, ylab=ylab, ylim=ylim, ...)
    b[,i] <- matrix(seq(min(treat.coef), max(treat.coef), length=200), ncol=1)
    priorprobfit <- pinc.mt(b, th=object$theta.hat, rho.min=object$margprior.limits[1], rho.max=object$margprior.limits[2], includeX=rep(FALSE,length(b)))
    lines(b, priorprobfit)
    b[,i] <- bmean[i]
  }
  devAskNewPage(ask = FALSE)
}
)


# show()
setMethod("show", signature(object='cilfit'), function(object) {
  cat('cilfit object with outcome of type',object$msfit$outcometype,',',object$msfit$p,'covariates and',object$msfit$family,'error distribution\n')
  cat(" Use coef() to obtain BMA estimates for treatment variables\n")
  cat(" Method available for 'msfit' objects can be applied to element 'msfit', e.g. coef(object$msfit) returns BMA estimates for all variables\n")
  cat(" Elements $margpp and $margprior contain marginal variable posterior and prior inclusion probabilities\n")
}
)


# postProb()
setMethod('postProb', signature(object = 'cilfit'),
  function(object, nmax, method = 'norm') {
    ans <- object[['model.postprobs']]
    if (missing(nmax)) {
      nmax <- nrow(ans)
    }
    return(ans[1:nmax, ])
  }
)

# coef()
# Returns BMA inference for treatment variables. For inference on all parameters use coef(object$msfit)
coef.cilfit <- function(object, ...) {
  return(object$cil.teff)
  #qls <- c(0.025, 0.975)
  #pes <- as.matrix(object[['cil.teff']], ncol = 1)
  #if (length(object[['cil.teff']]) == 1) {
  #  cis <- quantile(object[['cil.teff.postdist']], probs = qls)
  #} else {
  #  cis <- apply(object[['cil.teff.postdist']], 2, quantile, probs = qls)
  #}
  #mpp <- as.matrix(object[['marg.postprobs']][1:length(pes)], ncol = 1)
  #out <- cbind(pes, t(as.matrix(cis, nrow = 1, ncol = 2)), mpp)
  #colnames(out) <- c('estimate', '2.5%', '97.5%', 'margpp')
  #rownames(out) <- paste('treat', 1:nrow(out), sep = '')
  #return(out)
}

# # coef()
# setMethod('coef', signature(object = 'cilfit'), function(object) {
#   qls <- c(0.025, 0.975)
#   pes <- as.matrix(object[['cil.teff']], ncol = 1)
#   if (length(object[['cil.teff']]) == 1) {
#     cis <- quantile(object[['cil.teff.postdist']], probs = qls)
#   } else {
#     cis <- apply(object[['cil.teff.postdist']], 2, quantile, probs = qls)
#   }
#   mpp <- as.matrix(object[['marg.postprobs']][1:length(pes)], ncol = 1)
#   out <- cbind(pes, t(as.matrix(cis, nrow = 1, ncol = 2)), mpp)
#   colnames(out) <- c('estimate', '2.5%', '97.5%', 'margpp')
#   rownames(out) <- paste('treat', 1:nrow(out), sep = '')
#   return(out)
# })
# END OF SCRIPT
