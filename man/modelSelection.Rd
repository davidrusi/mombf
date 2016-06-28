\name{modelSelection}
\alias{modelSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian variable selection for linear models via non-local priors. }
\description{
  Bayesian model selection for linear models using
  non-local priors. The
  algorithm uses a Gibbs scheme and can handle p>>n cases.
  See \code{rnlp} to obtain posterior samples for the coefficients.
}
\usage{
modelSelection(y, x, center=TRUE, scale=TRUE, niter=10^4, thinning=1,
burnin=round(niter/10), family='normal', priorCoef=momprior(tau=0.348),
priorDelta=modelbbprior(alpha.p=1,beta.p=1),
priorVar=igprior(alpha=.01,lambda=.01),
priorSkew=momprior(tau=0.348), phi, deltaini=rep(FALSE,ncol(x)),
initSearch='greedy', method='auto', optimMethod='LMA', B=10^5, verbose=TRUE)
}
\arguments{
\item{y}{Vector with observed responses}
\item{x}{Design matrix with all potential predictors}
\item{center}{If \code{center==TRUE}, \code{y} and \code{x} are centered to have zero mean, therefore eliminating the need to include an intercept term in x}
\item{scale}{If \code{scale==TRUE}, \code{y} and columns in \code{x} are scaled to have standard deviation 1}
\item{niter}{Number of Gibbs sampling iterations}
\item{thinning}{MCMC thinning factor, i.e. only one out of each \code{thinning} iterations are reported. Defaults to thinning=1, i.e. no thinning}
\item{burnin}{Number of burn-in MCMC iterations. Defaults to
  \code{.1*niter}. Set to 0 for no burn-in}
\item{family}{Residual distribution. Possible values are
  'normal','twopiecenormal','laplace', 'twopiecelaplace', or 'auto'. For
  the latter the residual distribution is inferred from the data.}
\item{priorCoef}{Prior distribution for the coefficients. Object created
  with \code{momprior}, \code{imomprior}, \code{emomprior} or \code{zellnerprior}.}
\item{priorDelta}{Prior on model space. Use \code{modelbbprior()} for Beta-Binomial prior, \code{modelbinomprior(p)} for Binomial prior with prior inclusion probability \code{p}, or \code{modelunifprior()} for Uniform prior}
\item{priorVar}{Inverse gamma prior on residual variance, created with
  \code{igprior()}}
\item{priorSkew}{prior on residual skewness parameter, assumed to be of
  the same family as priorCoef. Ignored if \code{family} is 'normal' or 'laplace'.}
\item{phi}{Residual variance. Typically this is unknown and therefore left missing. If specified argument priorVar is ignored}
\item{deltaini}{Logical vector of length \code{ncol(x)} indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model}
\item{initSearch}{Algorithm to refine \code{deltaini}. \code{initSearch=='greedy'} uses a greedy Gibbs sampling search. \code{initSearch=='SCAD'} sets \code{deltaini} to the non-zero elements in a SCAD fit with cross-validated regularization parameter. \code{initSearch=='none'} leaves \code{deltaini} unmodified}
\item{method}{Method to compute marginal likelihood. The default is to
  used closed-form expressions whenever possible (only available for
  pMOM priors with up to 15 covariates) and Laplace approximations otherwise.
  \code{method=='Laplace'} for Laplace approx, \code{method=='plugin'}
  for BIC-type approximation,
  \code{method=='MC'} for Importance Sampling, \code{method=='Hybrid'}
  for Hybrid Laplace-IS (only available for piMOM prior). See Details.}
\item{hess}{Method to estimat the hessian in the Laplace approximation to the integrated
  likelihood under Laplace or asymmetric Laplace errors. When
  hess=='asymp' the asymptotic hessian is used, hess=='asympDiagAdj' a
  diagonal adjustment is applied (see Rossell and Rubio for details).}
\item{optimMethod}{Algorithm to maximize objective function when method=='Laplace' or method=='MC'. Only used for family=='twopiecenormal'. optimMethod=='LMA' uses modified Newton-Raphson algorithm, 'CDA' coordinate descent algorithm}
\item{B}{Number of samples to use in Importance Sampling scheme. Ignored if \code{method=='Laplace'}}
\item{verbose}{Set \code{verbose==TRUE} to print iteration progress}
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
\item{coef}{Estimated coefficients (via posterior mode) for highest
  posterior probability model}
}
\references{
  Johnson V.E., Rossell D. Non-Local Prior Densities for Default
  Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B,
  2010, 72, 143-170.
  
  Johnson V.E., Rossell D. Bayesian model selection in high-dimensional
  settings. Journal of the American Statistical Association, 2012, 107,
  649-660.

  Rossell D., Rubio F.J. Tractable Bayesian variable selection: beyond
  normal errors. 2016
}
\author{ David Rossell }
\seealso{
 See \code{\link{msfit-class}} for details on the output, postProb to
 extract posterior model probabilities.
 See \code{\link{rnlp}} to obtain posterior samples for the coefficients.
 To compute marginal densities for a given model see
  \code{\link{pmomMarginalK}}, \code{\link{pmomMarginalU}},
  \code{\link{pimomMarginalK}}, \code{\link{pimomMarginalU}}.
}
\examples{
#Simulate data
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x \%*\% theta + rnorm(100)

#Specify prior parameters
priorCoef <- imomprior(tau=1)
priorDelta <- modelunifprior()
priorVar <- igprior(alpha=.01,lambda=.01)

#Alternative prior for model space: 0.5 prior prob for including any covariate
priorDelta <- modelbinomprior(p=0.5)

#Alternative prior for model space: Beta hyper-prior for prob of inclusion
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)

#Model selection
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE, niter=10^2,
priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar, phi=1,
method='Laplace')
fit1$postMode
fit1$margpp
postProb(fit1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

