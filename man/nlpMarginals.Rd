\name{nlpmarginals}
\alias{nlpmarginals}
\alias{nlpMarginal}
\alias{pimomMarginalK}
\alias{pimomMarginalU}
\alias{pmomMarginalK}
\alias{pmomMarginalU}
\title{ Marginal density of the observed data for linear regression with
  Normal, two-piece Normal, Laplace or two-piece Laplace residuals
  under non-local and Zellner priors }
\description{
  The marginal density of the data, i.e. the likelihood
  integrated with respect to the given prior distribution on the
  regression coefficients of the variables included in the model and
  an inverse gamma prior on the residual variance.

  \code{nlpMarginal} is the general function, the remaining ones
  correspond to particular cases and are kept for backwards
  compatibility with old code, and will be deprecated in the future.
}
\usage{
nlpMarginal(sel, y, x, data, smoothterms, nknots=9, groups=1:ncol(x),
family="normal", priorCoef, priorGroup,
priorVar=igprior(alpha=0.01,lambda=0.01), priorSkew=momprior(tau=0.348),
phi, method='auto', adj.overdisp='intercept', hess='asymp', optimMethod,
optim_maxit, initpar='none', B=10^5, logscale=TRUE, XtX, ytX)

pimomMarginalK(sel, y, x, phi, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX)
pimomMarginalU(sel, y, x, alpha=0.001, lambda=0.001, tau=1,
method='Laplace', B=10^5, logscale=TRUE, XtX, ytX)
pmomMarginalK(sel, y, x, phi, tau, r=1, method='auto', B=10^5,
logscale=TRUE, XtX, ytX)
pmomMarginalU(sel, y, x, alpha=0.001, lambda=0.001, tau=1,
r=1, method='auto', B=10^5, logscale=TRUE, XtX, ytX)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sel}{Vector with indexes of columns in x to be included in the model.
  Ignored if \code{y} is a formula}
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
\item{family}{Residual distribution. Possible values are
  'normal','twopiecenormal','laplace', 'twopiecelaplace'}
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
\item{priorVar}{Inverse gamma prior on scale parameter, created by
  \code{igprior(). For Normal variance=scale, for Laplace variance=2*scale.}}
\item{priorSkew}{Either a number fixing tanh(alpha) where alpha is the
  asymmetry parameter or a prior on residual skewness parameter,
  assumed to be of
  the same family as priorCoef. Ignored if \code{family} is 'normal' or
  'laplace'.}
\item{method}{Method to approximate the integral. See
  \code{help(modelSelection)}.}
\item{adj.overdisp}{Only used for method=='ALA'. Over-dispersion adjustment for models with fixed
  dispersion parameter such as logistic and Poisson regression}
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
  posterior mode to approximate the integrated likelihood. See
  help(modelSelection)}
\item{B}{Number of Monte Carlo samples to use (ignored unless
  \code{method=='MC'})}
\item{logscale}{If \code{logscale==TRUE} the log marginal density is returned.}
\item{XtX}{Optionally, specify the matrix X'X.
  Useful when the function must be called a large number of times.}
\item{ytX}{Optionally, specify the vector y'X.
  Useful when the function must be called a large number of times.}
\item{phi}{Disperson parameter. See \code{help(modelSelection)}}
\item{alpha}{Prior for phi is inverse gamma \code{alpha/2},
  \code{lambda/2}}
\item{lambda}{Prior for phi is inverse gamma \code{alpha/2},
  \code{lambda/2}}
\item{tau}{Prior dispersion parameter for MOM and iMOM priors (see details)}
\item{r}{Prior power parameter for MOM prior is \code{2*r}}
}
\details{
  The marginal density of the data is equal to the integral of
  N(y;x[,sel]*theta,phi*I) * pi(theta|phi,tau) * IG(phi;alpha/2,lambda/2)
  with respect to theta, where pi(theta|phi,tau) is a non-local prior
  and IG denotes the density of an inverse gamma.

  \code{pmomMarginalK} and \code{pimomMarginalK} assume that the
  residual variance is known and therefore the inverse-gamma term in the
  integrand can be ommitted.

  The product MOM and iMOM densities can be evaluated using the
  functions \code{dmom} and \code{dimom}.
}
\value{
  Marginal density of the observed data under the specified prior.
}
\references{
  Johnson V.E., Rossell D. Non-Local Prior Densities for Default Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B, 2010, 72, 143-170.
  See http://rosselldavid.googlepages.com for technical reports.
}
\author{ David Rossell }
\seealso{ \code{\link{modelSelection}} to perform model selection based
  on product non-local priors.
  \code{\link{momunknown}}, \code{\link{imomunknown}}, \code{\link{momknown}}, \code{\link{imomknown}}
  to compute Bayes factors for additive MOM and iMOM priors }
\examples{
x <- matrix(rnorm(100*2),ncol=2)
y <- x \%*\% matrix(c(.5,1),ncol=1) + rnorm(nrow(x))
pmomMarginalK(sel=1, y=y, x=x, phi=1, tau=1, method='Laplace')
pmomMarginalK(sel=1:2, y=y, x=x, phi=1, tau=1, method='Laplace')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

