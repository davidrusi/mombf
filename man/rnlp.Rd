\name{rnlp}
\alias{rnlp}
\alias{rnlp-methods}
\alias{rnlp,missing,missing,missing,missing,msfit,missing,missing-method}
\alias{rnlp,ANY,matrix,missing,missing,missing,missing,missing-method}
\alias{rnlp,ANY,matrix,missing,missing,msfit,missing,missing-method}
\alias{rnlp,missing,missing,numeric,matrix,missing,missing,missing-method}
\alias{rnlp,ANY,matrix,missing,missing,missing,character,character-method}
\title{
Posterior sampling for regression parameters
}
\description{
Gibbs sampler for linear, generalized linear and survival models
under product non-local priors, Zellner's prior and a Normal
approximation to the posterior.
Both sampling conditional on a model and Bayesian model averaging are
implemented (see Details).

If x and y not specified samples from non-local priors/posteriors with
density proportional to d(theta) N(theta; m, V) are produced,
where d(theta) is the non-local penalty term.
}
\usage{
rnlp(y, x, m, V, msfit, outcometype, family, priorCoef, priorGroup,
priorVar, priorprec, isgroup,
niter=10^3, burnin=round(niter/10), thinning=1, pp='norm')
}
\arguments{
  \item{y}{Vector with observed responses. When \code{class(y)=='Surv'}
    sampling is based on the Cox partial likelihood, else a linear model
  is assumed.}
  \item{x}{Design matrix with all potential predictors}
  \item{m}{Mean for the Normal kernel}
  \item{V}{Covariance for the Normal kernel}
  \item{msfit}{Object of class \code{msfit} returned by
    \code{modelSelection}. If specified Bayesian model averaging
    posterior samples are returned, according to posterior model
    probabilities in \code{msfit}, and then arguments \code{y},
    \code{x}, \code{m}, \code{V} etc. If \code{msfit} is missing then
    posterior samples under the full model \code{y ~ x} are returned}
  \item{outcometype}{Type of outcome. Possible values are "Continuous",
    "glm" or "Survival"}
  \item{family}{Assumed family for the family. Some possible values are
    "normal", "binomial logit" and "Cox"}
  \item{priorCoef}{Prior distribution for the coefficients. Ignored if
    \code{msfit} is supplied. Must be object of class
    \code{msPriorSpec}, e.g. created by \code{momprior},
    \code{emomprior}, \code{imomprior}, \code{zellnerprior}}
  \item{priorGroup}{Prior on grouped coefficients (e.g. categorical
  predictors with >2 categories, splines), as passed to \code{modelSelection}}
  \item{priorVar}{Prior on residual variance. Ignored if \code{msfit}
    supplied. Must be object of class \code{msPriorSpec}, e.g. created
    with \code{igprior}}
  \item{priorprec}{Prior precision matrix (inverse
    covariance). Currently, this is only needed for the ICAR+ prior}
  \item{isgroup}{Logical vector where \code{TRUE} indicates that the
    variable is part of a group, e.g. one of several dummy indicators
    for a discrete covariate}
  \item{niter}{Number of MCMC iterations}
  \item{burnin}{Number of burn-in MCMC iterations. Defaults to \code{.1*niter}. Set to 0 for no burn-in}
  \item{thinning}{MCMC thinning factor, i.e. only one out of each \code{thinning} iterations are reported. Defaults to no thinning}
  \item{pp}{When \code{msfit} is provided this is the method to compute posterior model probabilities,
  which determine the sampled models. Can be 'norm' or 'exact', see \code{postProb} for details. }
}
\details{
  The algorithm is implemented for product MOM (pMOM), product iMOM
  (piMOM) and product eMOM (peMOM) priors.
  The algorithm combines an
  orthogonalization that provides low serial correlation
  with a latent truncation representation that allows fast sampling.

  When \code{y} and \code{x} are specified sampling is for the linear
  regression posterior.
  When argument \code{msfit} is left missing, posterior sampling is for
  the full model regressing \code{y} on all covariates in \code{x}.
  When \code{msfit} is specified each model is drawn with
  probability given by \code{postProb(msfit)}. In this case, a Bayesian
  Model Averaging estimate of the regression coefficients can be
  obtained by applying \code{colMeans} to the \code{rnlp} ouput matrix.

  When \code{y} and \code{x} are left missing, sampling is from a
  density proportional to d(theta) N(theta; m,V), where d(theta) is the
  non-local penalty (e.g. d(theta)=prod(theta^(2r)) for the pMOM prior).
}
\value{
Matrix with posterior samples
}
\references{
  D. Rossell and D. Telesca. Non-local priors for high-dimensional
  estimation, 2014. http://arxiv.org/pdf/1402.5107v2.pdf
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
#Simulate data
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x \%*\% theta + rnorm(100)
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE)

th <- rnlp(msfit=fit1, niter=100)
colMeans(th)
}
\keyword{ models }
\keyword{ distribution }
