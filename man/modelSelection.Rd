\name{modelSelection}
\alias{modelSelection}
\alias{modelsearchBlockDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian variable selection for linear models via non-local priors. }
\description{
  Bayesian model selection for linear, asymmetric linear,
  median and quantile regression under
  non-local or Zellner priors. p>>n can be handled.

  modelSelection enumerates all models when feasible
  and uses a Gibbs scheme otherwise.
  See \code{coef} and \code{coefByModel} for estimates and posterior
  intervals of regression coefficients, and \code{rnlp} for posterior samples.

  modelsearchBlockDiag seeks the highest posterior
  probability model using an iterative block search.
}
\usage{
modelSelection(y, x, data, smoothterms, nknots=9,
groups=1:ncol(x), constraints, center=TRUE, scale=TRUE,
enumerate, includevars=rep(FALSE,ncol(x)), models,
maxvars, niter=5000, thinning=1,
burnin=round(niter/10), family='normal', priorCoef,
priorGroup, priorDelta=modelbbprior(1,1),
priorConstraints,
priorVar=igprior(.01,.01),
priorSkew=momprior(tau=0.348), phi, deltaini=rep(FALSE,ncol(x)),
initSearch='greedy', method='auto', adj.overdisp='intercept',
hess='asymp', optimMethod, optim_maxit, initpar='none', B=10^5,
XtXprecomp= ifelse(ncol(x)<10^4,TRUE,FALSE), verbose=TRUE)

modelsearchBlockDiag(y, x, priorCoef=momprior(tau=0.348),
priorDelta=modelbbprior(1,1), priorVar=igprior(0.01,0.01),
blocksize=10, maxiter=10, maxvars=100, maxlogmargdrop=20,
maxenum=10, verbose=TRUE)
}
\arguments{
\item{y}{Either a formula with the regression equation or a vector with
  observed responses. The response can be either continuous or of class
  \code{Surv} (survival outcome). If \code{y} is a formula then \code{x},
  \code{groups} and \code{constraints} are automatically created}
\item{x}{Design matrix with linear covariates for which we want to
  assess if they have a linear effect on the response. Ignored if
  \code{y} is a formula}
\item{data}{If \code{y} is a formula then \code{data} should be a data
  frame containing the variables in the model}
\item{smoothterms}{Formula for non-linear covariates (cubic splines),
  modelSelection assesses if the variable has no effect, linear or
  non-linear effect. \code{smoothterms} can also be a design matrix or
  data.frame containing linear terms, for each column modelSelection
  creates a spline basis and tests no/linear/non-linear effects}
\item{nknots}{Number of spline knots. For cubic splines the non-linear
  basis adds knots-4 coefficients for each linear term, we recommend
  setting \code{nknots} to a small/moderate value}
\item{groups}{If variables in \code{x} such be added/dropped in groups,
  \code{groups} indicates the group that each variable corresponds to
  (by default each variable goes in a separate group)}
\item{constraints}{Constraints on the model space. List with length
  equal to the number of groups;
  if group[[i]]=c(j,k) then group i can only be in the model if groups j and k are also in the model}
\item{center}{If \code{TRUE}, \code{y} and \code{x} are centered to have
  zero mean. Dummy variables corresponding to factors are NOT centered}
\item{scale}{If \code{TRUE}, \code{y} and columns in \code{x} are
  scaled to have variance=1. Dummy variables corresponding to factors are NOT scaled}
\item{enumerate}{Default is \code{TRUE} if there's less than 15 variable
  groups. If \code{TRUE} all models with up to \code{maxvars} are
  enumerated, else Gibbs sampling is used to explore the model space}
\item{includevars}{Logical vector of length ncol(x) indicating variables
that should always be included in the model, i.e. variable selection is
not performed for these variables}
\item{models}{Optional logical matrix indicating the models to be
  enumerated with rows equal to the number of desired models and columns
to the number of variables in \code{x}. }
\item{maxvars}{When \code{enumerate==TRUE} only models with up to maxvars variables
  enumerated (defaults to all variables). In \code{modelsearchBlockDiag} a sequence of models
  is defined from 1 up to \code{maxvars}}
\item{niter}{Number of Gibbs sampling iterations}
\item{thinning}{MCMC thinning factor, i.e. only one out of each \code{thinning} iterations are reported. Defaults to thinning=1, i.e. no thinning}
\item{burnin}{Number of burn-in MCMC iterations. Defaults to
  \code{.1*niter}. Set to 0 for no burn-in}
\item{family}{Family of parametric distribution. Use
  'normal' for Normal errors, 'binomial' for logistic regression,
  'poisson' for Poisson regression.
  'twopiecenormal' for two-piece Normal,
  'laplace' for Laplace errors and 'twopiecelaplace' for double
  exponential.
  For 'auto' the errors are assumed continuous and their distribution
  is inferred from the data among
  'normal', 'laplace', 'twopiecenormal' and 'twopiecelaplace'.
  'laplace' corresponds to median regression and 'twopiecelaplace'
  to quantile regression. See argument \code{priorSkew}}
\item{priorCoef}{Prior on coefficients, created
  by \code{momprior}, \code{imomprior}, \code{emomprior} or
  \code{zellnerprior}.
  Prior dispersion is on coefficients/sqrt(scale) for Normal and
  two-piece Normal, and on coefficients/sqrt(2*scale) for Laplace
  and two-piece Laplace.}
\item{priorGroup}{Prior on grouped coefficients (e.g. categorical
  predictors with >2 categories, splines). Created by
  \code{groupmomprior}, \code{groupemomprior},
  \code{groupimomprior} or \code{groupzellnerprior}}
\item{priorDelta}{Prior on model space. Use \code{modelbbprior()}
  for Beta-Binomial prior, \code{modelbinomprior(p)} for Binomial
  prior with prior inclusion probability \code{p},
  \code{modelcomplexprior} for Complexity prior,
  or \code{modelunifprior()} for Uniform prior}
\item{priorConstraints}{Prior distribution on the number of terms
  subject to hierarchical constrains that are included in the model}
\item{priorVar}{Inverse gamma prior on scale parameter.
  For Normal outcomes variance=scale, for Laplace outcomes
  variance=2*scale}
\item{priorSkew}{Either a fixed value for tanh(alpha) where alpha is
  the asymmetry parameter or a prior on tanh(alpha).
  For \code{family=='twopiecelaplace'} setting alpha=a is equivalent
  to performing quantile regression for the quantile (1+a)/2.
  Ignored if \code{family} is 'normal' or 'laplace'.}
\item{phi}{The error variance in Gaussian models, typically this is
  unknown and is left missing}
\item{deltaini}{Logical vector of length \code{ncol(x)} indicating which
  coefficients should be initialized to be non-zero.
  Defaults to all variables being excluded from the model}
\item{initSearch}{Algorithm to refine
  \code{deltaini}. \code{initSearch=='greedy'} uses a greedy Gibbs
  sampling search. \code{initSearch=='SCAD'} sets \code{deltaini} to the
  non-zero elements in a SCAD fit with cross-validated regularization
  parameter. \code{initSearch=='none'} leaves \code{deltaini} unmodified}
\item{method}{Method to compute marginal likelihood.
  \code{method=='Laplace'} for Laplace approx, \code{method=='ALA'}
  for approximate Laplace approximation.
  \code{method=='MC'} for Importance Sampling, \code{method=='Hybrid'}
  for Hybrid Laplace-IS (only available for piMOM prior). See Details.}
  \code{method=='auto'} attempts to use exact calculations when
  possible, otherwise ALA if available, otherwise Laplace approx.
\item{adj.overdisp}{Only used when method=='ALA'. Over-dispersion
  adjustment in models with fixed dispersion parameter, as in logistic
  and Poisson regression. adj.overdisp='none' for no adjustment (not
  recommended, particularly for Poisson
  models). adj.overdisp='intercept' to estimate over-dispersion from the
  intercept-only model, and adj.overdisp='residuals' from the Pearson
  residuals of each model}
\item{hess}{Method to estimat the hessian in the Laplace approximation to the integrated
  likelihood under Laplace or asymmetric Laplace errors. When
  hess=='asymp' the asymptotic hessian is used, hess=='asympDiagAdj' a
  diagonal adjustment is applied (see Rossell and Rubio for details).}
\item{optimMethod}{Algorithm to maximize objective function when
  method=='Laplace'. Leave unspecified or set optimMethod=='auto' for an
  automatic choice. optimMethod=='LMA' uses modified
  Newton-Raphson algorithm, 'CDA' coordinate descent algorithm}
\item{optim_maxit}{Maximum number of iterations when method=='Laplace'}
\item{initpar}{Initial regression parameter values when finding the
  posterior mode to approximate the integrated likelihood. 'none', 'MLE',
  'L1', or a numeric vector with initial
  values. 'auto': if p<n/2 MLE is used, else L1 (regularization parameter set
  via BIC)}
%\item{initpar}{Initial regression parameter values when finding the
%  posterior mode to approximate the integrated likelihood. 'none', 'MLE',
%  'MLE-aisgd', 'L1', 'L2-aisgd', or a numeric vector with initial
%  values. If p<n/2 MLE is used, else L1 (regularization parameter set
%  via BIC). 'auto': if n>10,000 or p>200, then MLE-aisgd or L2-aisgd are
%  used. aisgd stands for averaged intrinsic stochastic gradient descent
%  (see function sgd in package sgd)}
\item{B}{Number of samples to use in Importance Sampling scheme. Ignored
  if \code{method=='Laplace'}}
\item{XtXprecomp}{Set to \code{TRUE} to pre-compute the Gram matrix x'x
  upfront (saves time), to \code{FALSE} to compute and store elements
  only as needed (saves memory)}
\item{verbose}{Set \code{verbose==TRUE} to print iteration progress}

\item{blocksize}{Maximum number of variables in a block. Careful, the
  cost of the algorithm is of order \code{2^blocksize}}
\item{maxiter}{Maximum number of iterations, each iteration includes a
  screening pass to add and subtract variables}
\item{maxlogmargdrop}{Stop the sequence of models when the drop in log
  p(y|model) is greater than \code{maxlogmargdrop}. This option avoids
  spending unnecessary time exploring overly large models}
\item{maxenum}{If the posterior mode found has less than \code{maxenum}
  variables then do a full enumeration of all its submodels}
}
\details{
  Let delta be the vector indicating inclusion/exclusion of each
  column of x in the model. The Gibbs algorithm sequentially samples from the
  posterior of each element in delta conditional on all the remaining
  elements in delta and the data.
  To do this it is necessary to evaluate the marginal likelihood for any
  given model. These have closed-form expression for the MOM prior, but
  for models with >15 variables these are expensive to compute and
  Laplace approximations are used instead (for the residual variance a
  log change of variables is used, which improves the approximation).
  For other priors closed forms
  are not available, so by default Laplace approximations are used.
  For the iMOM prior we also implement a Hybrid Laplace-IS
  which uses a Laplace approximation to evaluate the integral wrt beta
  and integrates wrt phi (residual variance) numerically.

  It should be noted that Laplace approximations tend to under-estimate
  the marginal densities when the MLE for some parameter is very close
  to 0. That is, it tends to be conservative in the sense of excluding
  more variables from the model than an exact calculation would.

  Finally, method=='plugin' provides a BIC-type approximation that is
  faster than exact or Laplace methods, at the expense of some
  accuracy. In non-sparse situations where models with many variables
  have large posterior probability method=='plugin' can be substantially
  faster.

  For more details on the methods used to compute marginal densities see
  Johnson & Rossell (2012).

  \code{modelsearchBlockDiag} uses the block search method described in
  Papaspiliopoulos & Rossell. Briefly, spectral clustering is run on
  X'X to cluster variables into blocks of \code{blocksize} and
  subsequently the Coolblock algorithm is used to define a sequence
  of models of increasing size. The exact integrated likelihood
  is evaluated for all models in this path, the best model chosen,
  and the scheme iteratively repeated to add and drop variables
  until convergence.
}
\value{
Object of class \code{msfit}, which extends a list with elements
\item{postSample}{\code{matrix} with posterior samples for the model
  indicator. \code{postSample[i,j]==1}
  indicates that variable j was included in the model in the MCMC
  iteration i}
\item{postOther}{\code{postOther}
  returns posterior samples for parameters other than the model
  indicator, i.e. basically hyper-parameters.
  If hyper-parameters were fixed in the model specification, \code{postOther} will be empty.}
\item{margpp}{Marginal posterior probability for inclusion of each
  covariate. This is computed by averaging marginal post prob for
  inclusion in each Gibbs iteration, which is much more accurate than
  simply taking \code{colMeans(postSample)}}.
\item{postMode}{Model with highest posterior probability amongst all those visited}
\item{postModeProb}{Unnormalized posterior prob of posterior mode (log scale)}
\item{postProb}{Unnormalized posterior prob of each visited model (log
  scale)}
\item{priors}{List with priors specified when calling \code{modelSelection}}
}
\references{
  Johnson V.E., Rossell D. Non-Local Prior Densities for Default
  Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B,
  2010, 72, 143-170.

  Johnson V.E., Rossell D. Bayesian model selection in high-dimensional
  settings. Journal of the American Statistical Association, 2012, 107,
  649-660.

  Papaspiliopoulos O., Rossell, D. Scalable Bayesian variable selection
  and model averaging under block orthogonal design. 2016

  Rossell D., Rubio F.J. Tractable Bayesian variable selection: beyond
  normality. 2016
}
\author{ David Rossell }
\seealso{
  \code{\link{msfit-class}} for details on the output.
  \code{\link{postProb}} to obtain posterior model probabilities.
  \code{coef.msfit} for Bayesian model averaging estimates and
  intervals. \code{predict.msfit} for BMA estimates and intervals
  for user-supplied covariate values.
  \code{plot.msfit} for an MCMC diagnostic plot showing estimated
  marginal posterior inclusion probabilities vs. iteration number.
 \code{\link{rnlp}} to obtain posterior samples for the coefficients.
 \code{\link{nlpMarginal}} to compute marginal densities for a given model.


}
\examples{
#Simulate data
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x \%*\% theta + rnorm(100)

#Specify prior parameters
priorCoef <- momprior(tau=0.348)
priorDelta <- modelunifprior()

#Alternative model space prior: 0.5 prior prob for including any covariate
priorDelta <- modelbinomprior(p=0.5)

#Alternative: Beta-Binomial prior for model space
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)

#Model selection
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE,
priorCoef=priorCoef, priorDelta=priorDelta)
postProb(fit1) #posterior model probabilities

fit1$margpp #posterior marginal inclusion prob

coef(fit1) #BMA estimates, 95% intervals, marginal post prob
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

