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

################################################################################
# LASSO estimation with BIC
lasso.bic <- function(y, x, intercept = TRUE, standardize = TRUE,
  family = 'gaussian') {
# (c) David Rossell
################################################################################
  #require(glmnet)
  fit <- glmnet(x = x, y = y, family = family, alpha = 1,
    intercept = intercept, standardize = standardize)
  if (intercept == TRUE) {
    pred <- cbind(1, x) %*% rbind(fit[['a0']], fit[['beta']])
  } else {
    pred <- x %*% fit[['beta']]
  }
  
  transf <- function(x, family) {  
    if (family == 'binomial') { x <- exp(pred) / (1 + exp(pred)) }
    return(x)
  }
  pred <- transf(pred, family = family)

  n <- length(y)
  p <- colSums(fit[['beta']] != 0) + ifelse(intercept == TRUE, 1, 0)
  bic <- n * log(colSums((y - pred)^2) / length(y)) +
         n * (log(2 * pi) + 1) + log(n) * p
  sel <- which.min(bic)
  beta <- c(fit[['a0']][sel], fit[['beta']][, sel])
  ypred <- pred[, sel]
  ans <- list(coef = beta, ypred = ypred, lambda.opt = fit[['lambda']][sel],
    lambda = data.frame(lambda = fit[['lambda']], bic = bic, nvars = p))
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
pinc.mt <- function(x, th, rho.min = 0, rho.max = 1, squared = FALSE) {
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

  # End
  return(probs)
}

################################################################################
grid.opt.mt <- function(G0, betad, pj1, method = 'EB', ws = NA, th.grid,
  th.prior = 'unif', rho.min = 0, rho.max = 1, V = diag(2),
  ret.plotly = FALSE) {
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
    opt.th[, ncol(opt.th) + 1] <- exp(-apply(opt.th, 1, Of.EB.cil,
      betad = betad, G0 = G0, ws = ws, th.prior = th.prior,
      rho.min = rho.min, rho.max = rho.max, V = V))
  } else if (method == 'EP') {
    opt.th[, ncol(opt.th) + 1] <- exp(-apply(opt.th, 1, Of.EP.cil,
      betad = betad, pj1 = pj1, th.prior = th.prior,
      rho.min = rho.min, rho.max = rho.max, V = V))
  } else {
    stop('supported theta search methods are "EB" and "EP".')
  }

  # Optimal values
  th.opt <- opt.th[which.max(opt.th[, ncol(opt.th)]), -ncol(opt.th)]
  th.opt <- as.numeric(th.opt)
  names(th.opt) <- paste('th', 1:length(th.opt) - 1, sep = '')

  # Output
  return(list(th.opt = th.opt))
}#; grid.opt.mt <- compiler::cmpfun(grid.opt.mt)

################################################################################
Of.EB.cil <- function(x, betad, G0, ws, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
  p1 <- c(rep(1/2, nchar(G0[1]) - length(p1)), p1)
  sc <- sum(sapply(G0, function(m) {
    gm <- as.numeric(unlist(strsplit(m, '')))
    pg <- prod((p1 ** gm) * (1 - p1) ** (1 - gm))
    wg <- ws[which(G0 == m)]
    return(wg * pg)
  }))
  pt1 <- log(sum(sc))  # log p(y | theta)

  # log p(theta)
  pt2 <- 0

  return(-(pt1 + pt2))
}

################################################################################
Gf.EB.cil <- function(x, betad, G0, ws, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # Compute (unnormalised) model inclusion probabilities
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
  sD <- nchar(G0[1]) - length(p1)
  p1 <- c(rep(1/2, sD), p1)
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
  gs[1] <- sum((mpprobs - p1)[-(1:sD)])
  for (k in 2:length(gs)) {
    aux0 <- (ncol(betad) / (length(x) - 1) - 1)
    idxs <- c(k - 1)
    if (aux0 > 0) { idxs <- c(idxs, (length(x) - 1) + aux0 * (k - 2) + 1:aux0) }
    aux1 <- mpprobs - p1
    aux2 <- rowSums(abs(as.matrix(betad[, idxs])))
    gs[k] <- sum(aux1[-(1:sD)] * aux2)
  }

  # grad log p(theta)
  t0 <- t1 <- 0

  # - grad log p(theta | y) = -(grad log p(y | theta) + grad log p(theta))
  return(-(gs + t0))
}

################################################################################
Of.EP.cil <- function(x, betad, pj1, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
  p1 <- c(rep(1/2, length(pj1) - length(p1)), p1)
  fj <- pj1 * p1 + (1 - pj1) * (1 - p1)
  pt1 <- sum(log(fj))

  # log p(theta)
  pt2 <- 0

  # - log p(theta | y) \propto -(log p(y | theta) + log p(theta))
  return(-(pt1 + pt2))
}

################################################################################
Gf.EP.cil <- function(x, betad, nt, pj1, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # grad log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
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
check.parameter.format <- function(rho.min, rho.max, th.range, tau, max.mod,
  lpen, eps, bvs.fit0, th.EP, D) {
################################################################################
  # Error messages accounted for
  err1 <- 'if informed, argument "rho.min" must be a scalar in the interval (0, 1/2).'
  err2 <- 'if informed, argument "rho.max" must be a scalar in the interval (1/2, 1).'
  err3 <- 'if informed, argument "th.range" must be numeric and have at least two distinct values.'
  err4 <- 'if informed, argument "tau" must be a positive scalar.'
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

  # Parameter: "rho.max"
  if (! is.null(rho.max) | length(rho.max) > 1) {
    if (! is.numeric(rho.max) | length(rho.max) > 1) {
      stop(err2)
    } else if (! inside(rho.max, c(0, 1/2), include.ab = FALSE)) {
      stop(err2)
    }
  }

  # Parameter: "th.range"
  if (! is.null(th.range)) {
    if (! is.numeric(th.range) | length(rho.max) == 1) {
      stop(err3)
    } else if (! length(unique(th.range)) == length(th.range)) {
      stop(err3)
    }
  }

  # Parameter: "tau"
  if (! is.null(tau) | length(tau) > 1) {
    if (! is.numeric(tau) | length(tau) > 1) {
      stop(err4)
    } else if (sign(tau) != 1) {
      stop(err4)
    }
  }

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
    if (class(bvs.fit0) != 'msfit') {
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
model.pprobs.cil <- function(y, D, X, I = NULL, mod1 = 'ginv', th.search = 'EB',
  th.prior = 'unif', beta.prior = 'nlp', rho.min = NULL, rho.max = NULL,
  th.range = NULL, tau = 0.348, max.mod = Inf, lpen = 'lambda.1se', eps = 1e-10,
  R = 1e5, bvs.fit0 = NULL, th.EP = NULL) {
################################################################################
  # Make sure parameter inputs are in the correct format
  check.parameter.format(rho.min = rho.min, rho.max = rho.max, tau = tau,
    th.range = th.range, max.mod = max.mod, lpen = lpen, eps = eps,
    bvs.fit0 = bvs.fit0, th.EP = th.EP, D = D)

  # Renaming (DEPENDENCIES: mombf, pracma, glmnet, hdm)
  ms <- modelSelection#mombf::modelSelection
  igp <- igprior(0.01, 0.01)#mombf::igprior(0.01, 0.01)
  bbp <- modelbbprior(1, 1)#mombf::modelbbprior(1, 1)
  nlp <- momprior(tau = tau)#mombf::momprior(tau = tau)
  zup <- zellnerprior(tau = nrow(X))#mombf::zellnerprior(tau = nrow(X))
  ufp <- modelunifprior()#mombf::modelunifprior()
  mlik <- nlpMarginal#mombf::nlpMarginal
  ginv <- pinv#pracma::pinv

  # These parameters must be integers 
  max.mod <- round(max.mod, 0)
  R <- round(R, 0)

  # Pre-computations
  p <- ncol(X)
  Z <- cbind(D, I, X)
  Di <- cbind(D, I)

  # Other fixed parameters
  if (is.null(rho.min)) {
    rho.min <- 1 / (ncol(Z)^2 + 1)
  }
  if (TRUE) {
  #if (is.null(rho.max)) {
    rho.max <- 1 - rho.min
  }

  # Coefficient prior
  if (beta.prior == 'nlp') {
    cp <- nlp
  } else if (beta.prior == 'zup') {
    cp <- zup
  } else {
    stop('arg. "beta.prior" currently supports: "nlp" (MOM) and "zup" (UIP)')
  }

  # MODULE 1: Estimate coefficients on exposure model
  betad <- matrix(NA, nrow = ncol(X), ncol = ncol(Di))
  for (i in 1:ncol(Di)) {
    # Define matrices
    A <- X
    b <- Di[, i]
    exc <- 2 * p  # phantom column: not exclude anything
    if (length(unique(Di[, i])) <= 3 &  # This controls sum-to-zero constraint
        all(unique(Di[, i]) %in% c(-1, 0, 1))) {
      b[which(b == -1)] <- 0
    }

    # Exposure family (support for Binomial and Gaussian ONLY)
    type.exp <- ifelse(length(unique(b)) == 2, 'binomial', 'gaussian')
    intcpt <- !any(apply(X, 2, var) == 0)

    if (length(table(b)) > 1 & min(table(b)) > 1) {
      if (mod1 == 'ginv') {
        betad[, i] <- as.matrix(ginv(t(A) %*% A) %*% t(A) %*% b)[-exc]
      } else if (mod1 %in% c('bvs', 'bma', 'bms')) {
        m1.fit <- ms(b, A, priorCoef = cp, priorVar = igp, priorDelta = bbp,
          niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)
        post.th <- colMeans(rnlp(msfit = m1.fit, priorCoef = cp, niter = 1e4))
        betad[, i] <- unname(post.th[2:(length(post.th) - 1)])[-exc]
      } else if (mod1 == 'lasso') {
        m1.fit <- cv.glmnet(A, b, intercept = intcpt)
        #m1.fit <- glmnet::cv.glmnet(A, b, intercept = intcpt)
        betad[, i] <- unname(coef(m1.fit, s = m1.fit[[lpen]])[-1, 1])[-exc]
      # } else if (mod1 == 'lasso_dml') {
      #   m1.fit <- hdm::rlasso(b ~ A, post = TRUE, intercept = intcpt)
      #   betad[, i] <- unname(m1.fit[['coefficients']])[-exc]
      } else if (mod1 == 'lasso_bic') {
        m1.fit <- lasso.bic(b, A, intercept = TRUE, family = type.exp)
        if (intcpt == FALSE) { m1.fit[['coef']][2] <- m1.fit[['coef']][1] }
        betad[, i] <- unname(m1.fit[['coef']][-1][-exc])
      } else if (mod1 == 'ridge') {
        m1.fit <- cv.glmnet(A, b, alpha = 0, intercept = intcpt)
        #m1.fit <- glmnet::cv.glmnet(A, b, alpha = 0, intercept = intcpt)
        betad[, i] <- unname(
          coef(m1.fit, s = m1.fit[['lambda.min']])[-1, ])[-exc]
      } else {
        stop('currently unsupported "mod1" method.')
      }
    } else {
      betad[, i] <- rep(0, nrow(betad))  # If treatment has no variability
    }
  }
  colnames(betad) <- colnames(Di)
  rownames(betad) <- colnames(X)

  # THETA SEARCH ###############################################################
  # Initial fit: th0 = 0; th1 = 0 (UNIFORM MODEL PRIOR)
  if (is.null(bvs.fit0)) {
    bvs.fit0 <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = ufp,
      niter = R, verbose = FALSE, center = FALSE, scale = FALSE)
  }
  # Posterior model probabilities
  pprobs <- postProb(bvs.fit0)[, c(1, 3)]

  # Marginal inclusion probabilities
  pj1 <- bvs.fit0[['margpp']]
  pj1[which(is.nan(pj1))] <- 0  # Potential problems in mombf with undefined
  pj1[which(pj1 == 1)] <- 1 - eps  # In case there are extreme values
  pj1[which(pj1 == 0)] <- eps

  # Set of visited models
  if (is.na(max.mod)) {  # Prevent a dictionary too big
    max.mod <- which(cumsum(pprobs[, ncol(pprobs)]) >= 0.9)[1]
  }
  G0 <- unname(sapply(
    as.character(pprobs[, 1]), binary.id, p = p, nt = ncol(Di), conc = TRUE))
  G0 <- G0[1:min(max.mod, length(G0))]  # Prevent a dictionary too big

  # Determine the values of hat{theta} #########################################
  # Set limits for optimisation (in principle: unbounded)
  lws <- rep(-Inf, ncol(D) + 1)
  ups <- rep(+Inf, ncol(D) + 1)
  s <- 1
  Rm <- diag(2)

  # Range for grid search
  if (is.null(th.range)) {
    if (ncol(D) == 1) {
      th.range <- seq(-40, 40, 1)
    } else if (ncol(D) == 2) {
      th.range <- seq(-40, 40, 2)
    } else if (ncol(D) <= 4) {
      th.range <- seq(-10, 10, 2)
    } else if (ncol(D) <= 7) {
      th.range <- seq(-7.5, 7.5, 2.5)
    } else if (ncol(D) <= 14) {
      th.range <- seq(-2, 2, 2)
      #th.range <- seq(-5, 5, 2.5)
    } else {
      th.range <- c(-2, 2)
    }
  }

  if (is.null(th.EP)) {
    # EP approximation (to set starting point)
    ths <- vector(ncol(D) + 1, mode = 'list')
    for (i in 1:length(ths)) { ths[[i]] <- th.range }
    EP.is <- grid.opt.mt(G0, betad, pj1, th.grid = expand.grid(ths),
      method = 'EP', th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
      V = s * Rm)

    # Load gradient functions and optimise
    st <- unlist(EP.is[[1]])
    opt.EP <- nlminb(st, objective = Of.EP.cil, gradient = Gf.EP.cil, 
      lower = lws, upper = ups, pj1 = pj1, betad = betad, 
      th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, V = s * Rm)

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
    auxprpr <- pinc.mt(betad, th = th.EP, rho.min = rho.min, rho.max = rho.max)
    auxprpr <- c(rep(1/2, ncol(betad)), auxprpr)
    auxprpr[which(auxprpr == 1)] <- 1 - eps
    auxprpr[which(auxprpr == 0)] <- eps
    auxmp <- modelbinomprior(p = auxprpr)#mombf::modelbinomprior(p = auxprpr)
    update.fit <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = auxmp,
      niter = round(R / 2), verbose = FALSE, center = FALSE, scale = FALSE)

    # Append model set
    upprobs <- postProb(update.fit)[, c(1, 3)]
    G1 <- unname(sapply(
      as.character(upprobs[, 1]), binary.id, p = p, nt = ncol(D), conc = TRUE))
    G0 <- unique(c(G0, G1[1:min(max.mod, length(G1))]))

    # Pre-compute marginal likelihoods
    ws <- unname(unlist(sapply(G0, function(m) {
      js <- which(unlist(strsplit(m, '')) == 1)
      return(mlik(js, y, Z, priorCoef = cp, logscale = FALSE))
    })))

    # Numerical problems with MLs that are too small
    if (all(ws < -650)) {  # We start having problems at exp(-750)
      ct <- abs(max(ws) + 650)  # Constant
      ws <- exp(ws + ct)  # Multiply by constant (NO EFFECT ON ARG MAX)
    } else {
      ws <- exp(ws)
    }

    # Optimise starting at the EP optimum
    opt.EB <- try(nlminb(th.EP, objective = Of.EB.cil, gradient = Gf.EB.cil,
      lower = lws, upper = ups, betad = betad, G0 = G0,
      ws = ws, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
      V = s * Rm), silent = TRUE)
    if (class(opt.EB) == 'try-error') {
      opt.EB <- try(nlminb(th.EP, objective = Of.EB.cil,
        gradient = Gf.EB.cil, lower = lws, upper = ups, betad = betad,
        G0 = G0, ws = ws, th.prior = th.prior, rho.min = rho.min,
        rho.max = rho.max, V = s * Rm))
    }

    # Optimal values
    if (opt.EB[['convergence']] %in% 0:1) {
      th.hat <- opt.EB[['par']]
    } else {
      stop('theta search under "EB" did not converge.')
    }
  } else if (th.search == 'EP') {
    # Optimal values
    th.hat <- th.EP
  } else {
    stop('theta search method unsupported -- choose amongst "EB" and "EP".')
  }

  # MODEL SEARCH (under fixed theta hat) #######################################
  # Prior probabilities under theta
  pg1cth <- c(rep(1/2, ncol(betad)),
              pinc.mt(betad, th = th.hat, rho.min = rho.min, rho.max = rho.max))
  pg1cth[which(pg1cth == 1)] <- 1 - eps
  pg1cth[which(pg1cth == 0)] <- eps

  # Define new asymmetric binomial prior
  nbp <- modelbinomprior(p = pg1cth)#mombf::modelbinomprior(p = pg1cth)

  # New model search
  new.fit <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = nbp,
    niter = R, verbose = FALSE, center = FALSE, scale = FALSE)

  # Posterior model probabilities
  pprobs <- postProb(new.fit)#mombf::postProb(new.fit)[, c(1, 3)]
  mpprobs <- new.fit[['margpp']]

  # Finish
  return(list(pprob.y = pprobs, mpprob.y = mpprobs, th.hat = th.hat,
              msfit = new.fit, init.msfit = bvs.fit0, th.EP = th.EP,
              mod1.coef = betad))
}#; model.pprobs.cil <- compiler::cmpfun(model.pprobs.cil)

################################################################################
bma.cil <- function(msfit, nt = 1, ret.mcmc = TRUE) {
################################################################################
  # Posterior draws
  draws <- rnlp(msfit = msfit, niter = 1e4)
  #draws <- mombf::rnlp(msfit = msfit, niter = 1e4)

  # BMA Treatment Effect Estimates
  idxs <- 1 + 1:nt
  if (nt == 1) {
    teff.est <- mean(draws[, idxs])
  } else {
    teff.est <- colMeans(draws[, idxs])
  }

  # BMA distribution
  if (nt == 1) {
    bma.distr <- quantile(draws[, idxs], seq(0, 1, 0.005))
  } else {
    bma.distr <- apply(draws[, idxs], 2, quantile, seq(0, 1, 0.005))
  }

  # Output
  out <- list(bma.e.cil = teff.est, bma.d.cil = bma.distr, mcmc.cil = NA)
  if (ret.mcmc == TRUE) { out[['mcmc.cil']] <- draws[, 1 + 1:nt] }

  # End
  return(invisible(out))
}

################################################################################
check.input.format <- function(y, D, X, I, R) {
################################################################################
  # Format errors that can be encountered
  err1 <- 'argument "y" must be of class "matrix" and contain one column.'
  err2 <- 'argument "D" must be of class "matrix" with at least one column and the same number of rows as "y".'
  err3 <- 'argument "X" must be of class "matrix" with at least one column and the same number of rows as "y".'
  err4 <- 'if argument "I" is informed, it must be of class "matrix" with at least one column and the same number of rows as "y".'
  err5 <- 'parameter "R" must be a single positive integer.'

  # Object "y"
  if (! 'matrix' %in% class(y)) {
    stop(err1)
  } else if (ncol(y) != 1) {
    stop(err1)
  }

  # Object "D"
  if (! 'matrix' %in% class(D)) {
    stop(err2)
  } else if (nrow(D) != nrow(y) | ncol(D) == 0) {
    stop(err2)
  }

  # Object "X"
  if (! 'matrix' %in% class(X)) {
    stop(err3)
  } else if (nrow(X) != nrow(y) | ncol(X) == 0) {
    stop(err3)
  }

  # Object "I"
  if (! is.null(I)) {
    if (! 'matrix' %in% class(I)) {
      stop(err4)
    } else if (nrow(I) != nrow(y) | ncol(I) == 0) {
      stop(err4)
    }
  }

  # Parameter "R"
  if (! is.numeric(R) | length(R) != 1) {
    stop(err5)
  } else if (R < 1) {
    stop(err5)
  }

  # Numeric objects
  if (! is.numeric(y)) { stop('"y" must be numeric.') }
  if (! is.numeric(D)) { stop('"D" must be numeric.') }
  if (! is.numeric(X)) { stop('"X" must be numeric.') }
  if (! is.null(I) & ! is.numeric(I)) {
    stop('if informed, "I" must be numeric.')
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
cil <- function(y, D, X, I = NULL, R = 1e4, th.search = 'EB',
  mod1 = 'lasso_bic', th.prior = 'unif', beta.prior = 'nlp', rho.min = NULL,
  rho.max = NULL, th.range = NULL, tau = 0.348, max.mod = Inf,
  lpen = 'lambda.1se', eps = 1e-10, bvs.fit0 = NULL, th.EP = NULL) {
################################################################################
  # Assert inputs are in the correct format
  check.input.format(y, D, X, I, R)

  # Posterior model probabilities
  pprobs <- model.pprobs.cil(y, D, X, I, R = R, th.search = th.search,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    rho.min = rho.min, rho.max = rho.max, th.range = th.range, tau = tau,
    max.mod = max.mod, lpen = lpen, eps = eps, bvs.fit0 = bvs.fit0,
    th.EP = th.EP)

  # BMA computation
  teff <- bma.cil(pprobs[['msfit']], nt = ncol(cbind(D, I)))

  # Output object
  out <- list(cil.teff = teff[['bma.e.cil']],  # BMA Point estimates
              cil.teff.postdist = teff[['bma.d.cil']],  # BMA Distribution
              cil.bma.mcmc = teff[['mcmc.cil']],  # BMA MCMC runs
              model.postprobs = pprobs[['pprob.y']],  # Model post. probs.
              marg.postprobs = pprobs[['mpprob.y']],  # Marginal post. probs.
              theta.hat = pprobs[['th.hat']],  # Estimated theta hat
              treat.coefs = pprobs[['mod1.coef']],  # Coefs. treatment models
              msfit = pprobs[['msfit']],  # Model selection fit of the model
              theta.EP = pprobs[['th.EP']],  # Estimated theta for EP method
              init.msfit = pprobs[['init.msfit']])  # Initial model sel. fit

  # Exit
  return(invisible(new('cilfit', out)))
  #return(invisible(out))
}

################################################################################
# Adapt relevant methods for objects of class "cilfit"
################################################################################
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
coef.cilfit <- function(object, ...) {
  qls <- c(0.025, 0.975)
  pes <- as.matrix(object[['cil.teff']], ncol = 1)
  if (length(object[['cil.teff']]) == 1) {
    cis <- quantile(object[['cil.teff.postdist']], probs = qls)
  } else {
    cis <- apply(object[['cil.teff.postdist']], 2, quantile, probs = qls)
  }
  mpp <- as.matrix(object[['marg.postprobs']][1:length(pes)], ncol = 1)
  out <- cbind(pes, t(as.matrix(cis, nrow = 1, ncol = 2)), mpp)
  colnames(out) <- c('estimate', '2.5%', '97.5%', 'margpp')
  rownames(out) <- paste('treat', 1:nrow(out), sep = '')
  return(out)
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
