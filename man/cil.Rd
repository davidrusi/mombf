\name{cil}
\alias{cil}
%\alias{modelsearchBlockDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Treatment effect estimation for linear models via Confounder Importance
  Learning using non-local priors.}
\description{
  Treatment effect estimation for linear models in the presence of
  multiple treatments and a potentially high-dimensional number of controls,
  i.e. \eqn{p \gg n} can be handled.

  Confounder Importance Learning (CIL) proposes an estimation framework where
  the importance of the relationship between treatments and controls is
  factored in into the establishment of prior inclusion probabilities for each
  of these controls on the response model. This is combined with the use of
  non-local priors to obtain BMA estimates and posterior model probabilities.

  \code{cil} is built on \code{modelSelection} and produces objects of type
  \code{cilfit}. Use \code{coef} and \code{postProb} to obtain treatment effect
  point estimates and posterior model probabilities, respectively, on this
  object class.
}
\usage{
cil(y, D, X, I = NULL, family = 'normal', familyD = 'normal',
  R = 1e4, Rinit = 500, th.search = 'EB', mod1 = 'lasso_bic',
  th.prior = 'unif', priorCoef = momprior(taustd=1),
  rho.min = NULL, rho.max = 0.95,
  th.range = NULL, max.mod = 2^20, lpen = 'lambda.1se',
  eps = 1e-10, bvs.fit0 = NULL, th.EP = NULL, center = TRUE, scale =
  TRUE, includevars, verbose = TRUE)
}
\arguments{
\item{y}{one-column matrix containing the observed responses. The response must
  be continuous (currently the only type supported)}
\item{D}{treatment matrix with numeric columns, continuous or discrete. Any finite
  number of treatments are supported. If only one treatment is provided, supply
  this object in the same format used for \code{y}}
\item{X}{matrix of controls with numeric columns, continuous or discrete. If only
  one treatment is provided, supply this object in the same format used for \code{y}}
\item{I}{matrix with the desired interaction terms between \code{D} and
  \code{X}. If not informed, i.e. supplied as the default \code{NULL}, this term
  will not be included into the response model}
\item{family}{Distribution of the outcome, e.g. 'normal', 'binomial' or
  'poisson'. See \code{help}(modelSelection) for a full list of options}
\item{familyD}{Distribution of the treatment(s). Only 'normal' or
  'binomial' currently allowed}
\item{R}{Number of MCMC iterations to be
  run by \code{modelSelection} on each stage of CIL (see argument \code{niter}
  therein)}
\item{Rinit}{MCMC iterations to estimate marginal posterior
  inclusion probabilities under a uniform model prior, needed for EP}
\item{th.search}{method to estimate theta values in the marginal prior inclusion
  probabilities of the CIL model. Options are: \code{EB} (Empirical Bayes, based
  on maximum marginal likelihood) and \code{EP} (Expectation propagation
  approximation)}
\item{mod1}{method to estimate the feature parameters corresponding to the
  influence of the controls on the treatments. Supported values for this
  argument are 'ginv' (generalised pseudo-inverse), \code{lasso} (see
  argument \code{lpen}), \code{lasso_bic} (default), and \code{ridge})}
\item{th.prior}{prior associated to the thetas for the Empirical Bayes
  estimation. Currently only \code{unif} (Uniform prior) is supported,
  effectively making the EB approach the maximisation of the marginal likelihood}
\item{priorCoef}{Prior on the response model parameters, see \code{modelSelection}}
\item{rho.min}{Lower bound on the covariate's prior inclusion probability. 
  By default, it is set to \eqn{1/p}, where p is the number of covariates}
\item{rho.max}{Upper bound on the covariate's prior inclusion probability}
\item{th.range}{sequence of values to be considered in the grid when searching
  for points to initialise the search for the optimal theta parameters. If left
  uninformed, the function will determine a computationally suitable grid
  depending on the number of parameters to be estimated}
\item{max.mod}{Maximum number of models considered when computing the marginal
  likelihood required by empirical Bayes.
  If set to \code{Inf} all visited models by the enumeration/MCMC
  are considered, but it might be computationally desirable to restrict this
  number when the dimension of \code{D} and/or \code{X} is large}
\item{lpen}{penalty type supplied to \code{glmnet} if \code{mod1} is set to
  \code{lasso}. Default is \code{lambda.1se} (see documentation corresponding to
  \code{glmnet} for options on how to set this parameter)}
\item{eps}{small scalar used to avoid round-offs to absolute zeroes or ones in
  marginal prior inclusion probabilities.}
\item{bvs.fit0}{object returned by \code{modelSelection} under \eqn{\theta = 0},
  used as a model exploration tool to compute EB approximation on the thetas.
  This argument is only supposed to be used in case of a second computation the
  model on the same data where \code{th.search} has ben changed to \code{EB},
  in order to avoid repeating the computation of the initial
  \code{modelSelection} fit. To use this argument, supply the object residing
  in the slot \code{init.msfit} of a \code{cilfit}-class object.}
\item{th.EP}{Optimal theta values under the EP approximation, obtained in a
  previous CIL run. This argument is only supposed to be used in case of
  a second computation the model on the same data where \code{th.search}
  has ben changed to \code{EB}, in order to save the cost of the EP search
  to initialise the optimisation algorithm. To use this argument, supply the
  object residing int the slot \code{th.hat} of a \code{cilfit}-class
  object.}
\item{center}{If \code{TRUE}, \code{y} and \code{x} are centered to have
  zero mean. Dummy variables corresponding to factors are NOT centered}
\item{scale}{If \code{TRUE}, \code{y} and columns in \code{x} are
  scaled to have variance=1. Dummy variables corresponding to factors are NOT scaled}
\item{includevars}{Logical vector of length ncol(x) indicating variables
that should always be included in the model, i.e. variable selection is
not performed for these variables}
\item{verbose}{Set \code{verbose==TRUE} to print iteration progress}
}
\details{
  We estimate treatment effects for the features present in the treatment
  matrix \code{D}. Features in \code{X}, which may or may not be causal
  factors of the treatments of interest, only act as controls and, therefore,
  are not used as inferential subjects.

  Confounder importance learning learns from data the amount of
  confounding in the data, based on a so-called confounding
  coefficient, and uses this to set covariate prior inclusion
  probabilities.
  Roughly, the coefficient measures to what extend covariates that are
  truly related to the outcome are also truly related (high confounding)
  or not related (no confounding) to the treatment(s).

  In high confounding, covariates in \code{X} that appear to be related
  to \code{D} are assigned high inclusion probability in the regression
  for \code{y}. In low confounding, they're assigned low prior inclusion
  probability. See function \code{plotprior} to visualize these prior
  probabilities. 
  Prior probabilities are regulated through a hyper-parameter theta that
  is set according to the method supplied to \code{th.search}.
  While the \code{EB} option obtains a more precise estimate,
  particularly when there are many covariates in \code{X}, the
  \code{EP} is much faster computationally and typically gives very
  similar results.

  See references for further details on implementation and computation.
}
\value{
  Object of class \code{cilfit}, which extends a list with elements
\item{cil.teff}{BMA estimates, 0.95 intervals and posterior inclusion
  probabilities for treatment effects in \code{D}}
\item{coef}{BMA inference for treatment effects and all other covariates}
\item{model.postprobs}{\code{matrix} returning the posterior model probabilities
  computed in the CIL model}
\item{margpp}{\code{numeric} vector containing the estimated marginal
  posterior inclusion probabilities of the featured treatments and controls}
\item{margprior}{Marginal prior inclusion probabilities, as estimated by
  CIL}
\item{margpp.unif}{Marginal posterior inclusion probabilities that would
  be obtained under a uniform model prior}
\item{theta.hat}{Values used for the hyper-parameter theta, estimated according
  to the argument \code{th.search} specified}
\item{treat.coefs}{Estimated weights of the effect of the control variables
  on each of the treatments, as estimated with the method specified in argument
  \code{mod1}}
\item{msfit}{Object returned by \code{modelSelection} (of class \code{msfit})
  of the final model estimated by CIL.}
\item{theta.EP}{Estimated values of theta using the EP algorithm. It coincides
  with \code{theta.hat} if the argument \code{th.search} is set to \code{EB}}
\item{init.msfit}{Initial \code{msfit} object used to estimate the inital model
  where all elements in theta are set to zero (used in the optimisation process
  of this hyper-parameter)}
}
\references{
  Torrens i Dinares M., Papaspiliopoulos O., Rossell D. Confounder
  importance learning for treatment effect inference.
  https://arxiv.org/abs/2110.00314, 2021, 1--48.
}
\author{
  Miquel Torrens
}
\seealso{
  \code{\link{postProb}} to obtain posterior model probabilities.
  
  \code{coef} for inference on the treatment parameters.

  \code{plotprior} to plot the estimated prior probabilities, as a
  function of 
  %\code{postProb} to obtain posterior model probabilities.
  %\code{coef.cilfit} for Bayesian model averaging estimates for the
  %treatment parameters.
}
\examples{
# Simulate data
set.seed(1)
X <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
beta_y <- matrix(c(rep(1, 6), rep(0, 44)), ncol = 1)
beta_d <- matrix(c(rep(1, 6), rep(0, 44)), ncol = 1)
alpha <- 1
d <- X \%*\% beta_d + rnorm(100)
y <- d * alpha + X \%*\% beta_y + rnorm(100)

# Confounder Importance Learning
fit1 <- cil(y = y, D = d, X = X, th.search = 'EP')

# BMA for treatment effects
coef(fit1)

# BMA for all covariates
head(fit1$coef)

# Estimated prior inclusion prob
# vs. treatment regression coefficients
plotprior(fit1)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ cil }

