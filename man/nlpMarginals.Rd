\name{nlpmarginals}
\alias{nlpmarginals}
\alias{pimomMarginalK}
\alias{pimomMarginalU}
\alias{pmomMarginalK}
\alias{pmomMarginalU}
\title{ Marginal density of the observed data for linear models under a
  given model, using non-local prior distributions for the regression coefficients.}
\description{
  The marginal density of the observed data, i.e. the likelihood
  integrated with respect to a non-local prior distribution on the
  regression coefficients of the variables included in the model and
  (optionally) an
  inverse gamma prior on the residual variance.
  \code{pimomMarginalK} and \code{pimomMarginalU} use a product iMOM prior on
  the regression coefficients for the known and unknown variance cases
  (respectively).
  \code{pmomMarginalK} and \code{pmomMarginalU} use a product MOM prior.
}
\usage{
pimomMarginalK(sel, y, x, phi, tau=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX) 
pimomMarginalU(sel, y, x, alpha=0.001, lambda=0.001, tau=1,
method='Laplace', B=10^5, logscale=TRUE, XtX, ytX)
pmomMarginalK(sel, y, x, phi, tau, r=1, method='Laplace', B=10^5,
logscale=TRUE, XtX, ytX)
pmomMarginalU(sel, y, x, alpha=0.001, lambda=0.001, tau=1,
r=1, method='Laplace', B=10^5, logscale=TRUE, XtX, ytX)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sel}{Vector with indexes of columns in x to be included in the model}
\item{y}{Vector with observed responses}
\item{x}{Design matrix with covariates. Only the columns specified in
  \code{sel} are included in the model, the rest are disregarded}
\item{phi}{Residual variance, assumed to be known by
  \code{pimomMarginalK} and \code{pmomMarginalK}}
\item{alpha}{Prior for phi is inverse gamma \code{alpha/2},
  \code{lambda/2}}
\item{lambda}{Prior for phi is inverse gamma \code{alpha/2},
  \code{lambda/2}}
\item{tau}{Prior dispersion parameter for MOM and iMOM priors (see details)}
\item{r}{Prior power parameter for MOM prior is \code{2*r}}
\item{method}{Method to approximate the integral. 'Laplace' for Laplace
approx. 'MC' for Monte Carlo. 'Plug-in' for plug-in estimate. 'Hybrid'
for hybrid Laplace-MC (only available for \code{pimomMarginalU}). See
help for \code{modelSelection} for details.}
\item{B}{Number of Monte Carlo samples to use (ignored unless
  \code{method=='MC'})}
\item{logscale}{If \code{logscale==TRUE} the log marginal density is returned.}
\item{XtX}{Optionally, specify the matrix X'X.
  Useful when the function must be called a large number of times.}
\item{ytX}{Optionally, specify the vector y'X.
  Useful when the function must be called a large number of times.}
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

