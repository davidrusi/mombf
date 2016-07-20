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
nlpMarginal(sel, y, x, family="normal", priorCoef=momprior(tau=0.348),
priorVar=igprior(alpha=0.01,lambda=0.01), priorSkew=momprior(tau=0.348),
method='auto', hess='asymp', optimMethod='CDA', B=10^5, logscale=TRUE, XtX, ytX)

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
\item{sel}{Vector with indexes of columns in x to be included in the model}
\item{y}{Vector with observed responses}
\item{x}{Design matrix with covariates. Only the columns specified in
  \code{sel} are included in the model, the rest are disregarded}
\item{family}{Residual distribution. Possible values are
  'normal','twopiecenormal','laplace', 'twopiecelaplace'}
\item{priorCoef}{Prior on coefficients, created
  by \code{momprior}, \code{imomprior}, \code{emomprior} or
  \code{zellnerprior}.
  Prior dispersion is on coefficients/sqrt(scale) for Normal and
  two-piece Normal, and on coefficients/sqrt(2*scale) for Laplace
  and two-piece Laplace.}
\item{priorVar}{Inverse gamma prior on scale parameter, created by
  \code{igprior(). For Normal variance=scale, for Laplace variance=2*scale.}}
\item{priorSkew}{prior on residual skewness parameter, assumed to be of
  the same family as priorCoef. Ignored if \code{family} is 'normal' or
  'laplace'.}
\item{method}{Method to approximate the integral. method=='auto' uses
  closed-form expressions whenever possible and Laplace approximations
  otherwise. method=='Laplace' for Laplace approx.
  method=='MC' for Monte Carlo}
\item{hess}{Method to estimat the hessian in the Laplace approximation to the integrated
  likelihood under Laplace or asymmetric Laplace errors. When
  hess=='asymp' the asymptotic hessian is used, hess=='asympDiagAdj' a
  diagonal adjustment is applied (see Rossell and Rubio for details).}
\item{optimMethod}{Algorithm to maximize objective function when
  method=='Laplace' or method=='MC'. optimMethod=='LMA' uses modified
  Newton-Raphson algorithm, 'CDA' coordinate descent algorithm}
\item{B}{Number of Monte Carlo samples to use (ignored unless
  \code{method=='MC'})}
\item{logscale}{If \code{logscale==TRUE} the log marginal density is returned.}
\item{XtX}{Optionally, specify the matrix X'X.
  Useful when the function must be called a large number of times.}
\item{ytX}{Optionally, specify the vector y'X.
  Useful when the function must be called a large number of times.}
\item{phi}{Residual variance, assumed to be known by
  \code{pimomMarginalK} and \code{pmomMarginalK}}
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
  to compute Bayes factors for additive MOM and iMOM priors.
  \code{\link{mode2g}} for prior elicitation.  }
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

