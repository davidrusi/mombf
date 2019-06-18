\name{dmom}
\alias{dmom}
\alias{dimom}
\alias{pmom}
\alias{pimom}
\alias{qmom}
\alias{qimom}
\alias{demom}
\alias{demom-methods}
\alias{demom,data.frame-method}
\alias{demom,matrix-method}
\alias{demom,vector-method}
\alias{pemom}
\alias{dmomigmarg}
\alias{pmomigmarg}
\alias{demomigmarg}
\alias{pemomigmarg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Non-local prior density, cdf and quantile functions. }
\description{
  \code{dmom}, \code{dimom} and \code{demom} return the density for the
  moment, inverse moment and exponential moment priors.
  \code{pmom}, \code{pimom} and \code{pemom} return the distribution function for the univariate
  moment, inverse moment and exponential moment priors (respectively).
  \code{qmom} and \code{qimom} return the quantiles for the univariate
  moment and inverse moment priors.
  \code{dmomigmarg} returns the marginal density implied by a
  MOM(x;tau*phi)*Invgamma(phi;a/2,b/2), \code{pmomigmarg} its cdf.
  Analogously \code{demomigmarg} and \code{demomigmarg} for
  eMOM(x;tau*phi)*Invgamma(phi;a/2,b/2)
}
\usage{
dmom(x, tau, a.tau, b.tau, phi=1, r=1, V1, baseDensity='normal', nu=3,
logscale=FALSE, penalty='product')
dimom(x, tau=1, phi=1, V1, logscale=FALSE, penalty='product')
demom(x, tau, a.tau, b.tau, phi=1, logscale=FALSE)

pmom(q, V1 = 1, tau = 1)
pimom(q, V1 = 1, tau = 1, nu = 1)
pemom(q, tau, a.tau, b.tau)

qmom(p, V1 = 1, tau = 1)
qimom(p, V1 = 1, tau = 1, nu = 1)

dmomigmarg(x,tau,a,b,logscale=FALSE)
pmomigmarg(x,tau,a,b)

demomigmarg(x,tau,a,b,logscale=FALSE)
pemomigmarg(x,tau,a,b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ In the univariate setting, \code{x} is a vector with the
    values at which to evaluate the density. In the multivariate setting
    it is a matrix with an observation in each row.}
  \item{q}{ Vector of quantiles.}
  \item{p}{ Vector of probabilities.}
  \item{V1}{ Scale matrix (ignored if \code{penalty=='product'}). Defaults to 1 in univariate setting and
    the identity matrix in the multivariate setting. }
  \item{tau}{ Prior dispersion parameter is \code{tau*phi}. See
    details. }
  \item{a.tau}{If \code{tau} is left missing, an Inverse Gamma(a.tau/2,b.tau/2)
  is placed on \code{tau}. In this case \code{dmom} and \code{demom} return the density
  marginalized with respect to \code{tau}.}
  \item{b.tau}{See \code{a.tau}.}
  \item{phi}{ Prior dispersion parameter is \code{tau*phi}. See
    details. }
  \item{r}{Prior power parameter for MOM prior is \code{2*r}}
  \item{baseDensity}{For \code{baseDensity=='normal'} a Normal MOM prior
    is used, for \code{baseDensity=='laplace'} a Laplace MOM prior,
    for \code{baseDensity=='t'} a T MOM prior with \code{nu}
  degrees of freedom is used.}
  \item{nu}{ Prior parameter indicating the degrees of freedom for the
    quadratic T MOM and iMOM prior densities. The
    tails of the inverse moment prior are proportional to the tails of a
    multivariate T with \code{nu} degrees of freedom.}
  \item{penalty}{\code{penalty=='product'} indicates that product MOM/iMOM should
    be used. \code{penalty=='quadratic'} indicates quadratic iMOM. See Details.}
  \item{logscale}{ For \code{logscale==TRUE}, \code{dimom} returns the
    natural log of the prior density.}
  \item{a}{The marginal prior on phi is IG(a/2,b/2)}
  \item{b}{The marginal prior on phi is IG(a/2,b/2)}
}
\details{
 For \code{type=='quadratic'} the density is as follows.
 Define the quadratic form q(theta)= (theta-theta0)' *
   solve(V1) * (theta-theta0) / (tau*phi).
   The normal moment prior density is proportional to
   q(theta)*dmvnorm(theta,theta0,tau*phi*V1).
   The T moment prior is proportional to
   q(theta)*dmvt(theta,theta0,tau*phi*V1,df=nu).
   The inverse moment prior density is proportional to
   \code{q(theta)^(-(nu+d)/2) * exp(-1/q(theta))}.

   pmom, pimom and qimom use closed-form expressions, while qmom uses
   nlminb to find quantiles numerically.
   Only the univariate version is implemented. In this case the product
   MOM is equivalent to the quadratic MOM. The same happens for the
   iMOM.

   \code{dmomigmarg} returns the marginal density

   p(x)= int MOM(x;0,tau*phi) IG(phi;a/2,b/2) dphi


}
\value{
  Prior density, cumulative distribution function or quantile.
}
\references{
    Johnson V.E., Rossell D. Non-Local Prior Densities for Default
  Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B,
  2010, 72, 143-170.

  Johnson V.E., Rossell D. Bayesian model selection in high-dimensional
  settings. Journal of the American Statistical Assocation, 2012, 107,
  649-660

See http://rosselldavid.googlepages.com for technical
  reports. }
\author{ David Rossell }
\examples{
#evaluate and plot the moment and inverse moment priors
library(mombf)
tau <- 1
thseq <- seq(-3,3,length=1000)
plot(thseq,dmom(thseq,tau=tau),type='l',ylab='Prior density')
lines(thseq,dimom(thseq,tau=tau),lty=2,col=2)
}
\keyword{ distribution }
\keyword{ htest }
