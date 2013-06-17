\name{rnlp}
\alias{rnlp}
\alias{rnlp-methods}
\alias{rnlp,numeric,matrix,missing-method}
\alias{rnlp,numeric,matrix,msfit-method}
\title{
Posterior sampling for Non-Local Priors
}
\description{
Gibbs sampler to draw posterior samples for linear model
regression coefficients and residual variance when the coefficient prior
is a product non-local prior.

Both sampling conditional on a model and Bayesian model averaging are
implemented (see Details).
}
\usage{
rnlp(y, x, msfit, priorCoef, priorVar, niter = 10^3, burnin = round(niter/10), 
thinning = 1)
}
\arguments{
  \item{y}{Vector with observed responses}
  \item{x}{Design matrix with all potential predictors}
  \item{msfit}{Object of class \code{msfit}, e.g. as returned by
    \code{modelSelection}. If left missing sampling under the full model
    \code{y ~ x} is performed, otherwise posterior model samples in
    \code{msfit} are used.}
  \item{priorCoef}{Prior distribution for the coefficients. Must be object of class \code{msPriorSpec} with slot \code{priorType} set to 'coefficients'. Possible values for slot \code{priorDistr} are 'pMOM', 'piMOM' and 'peMOM'}
  \item{priorVar}{Prior on residual variance. Must be object of class \code{msPriorSpec} with slot \code{priorType} set to 'nuisancePars'. Slot \code{priorDistr} must be equal to 'invgamma'}
  \item{niter}{Number of MCMC iterations}
  \item{burnin}{Number of burn-in MCMC iterations. Defaults to \code{.1*niter}. Set to 0 for no burn-in}
  \item{thinning}{MCMC thinning factor, i.e. only one out of each \code{thinning} iterations are reported. Defaults to no thinning}
}
\details{
  The algorithm is implemented for product MOM (pMOM), product iMOM
  (piMOM) and product eMOM (peMOM) priors.
  The algorithm combines an
  orthogonalization that provides low serial correlation
  with a latent truncation representation that allows fast sampling.

  When argument \code{msfit} is left missing, posterior sampling is for
  the full model regressing \code{y} on all covariates in \code{x}.
  
  When \code{msfit} is specified each model is drawn with
  probability given by \code{postProb(msfit)}. In this case, a Bayesian
  Model Averaging estimate of the regression coefficients can be
  obtained by applying \code{colMeans} to the \code{rnlp} ouput matrix.
}
\value{
Matrix with posterior samples
}
\references{
D. Telesca and D. Rossell. Technical report. 2013
}
\author{
David Rossell
}
\seealso{
\code{\link{modelSelection}} to perform model selection and compute
posterior model probabilities.
For more details on prior specification see \code{\link{msPriorSpec-class}}.
}
\examples{
#Generate data
set.seed(2)
n <- 10^3; tau <- 0.133; x <- rmvnorm(n,sigma=matrix(c(2,1,1,2),nrow=2))
thtrue <- c(0.5,1); phitrue <- 1
y <- thtrue[1]*x[,1] + thtrue[2]*x[,2] + rnorm(n,sd=sqrt(phitrue))

#Specify prior parameters
priorCoef <- imomprior(tau=1)
priorVar <- igprior(alpha=.01,lambda=.01)

th <- rnlp(y=y, x=x, niter=100, priorCoef=priorCoef, priorVar=priorVar)
colMeans(th)
acf(th[,1])[1]
}
\keyword{ models }
\keyword{ distribution }
