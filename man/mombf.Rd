\name{mombf}
\alias{mombf}
\alias{mombf.lm}
\alias{imombf}
\alias{imombf.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Moment and inverse moment Bayes factors for linear models. }
\description{
  \code{mombf} computes moment Bayes factors to test whether a subset of
  regression coefficients are equal to some user-specified value.
  \code{imombf} computes inverse moment Bayes factors.
}
\usage{
mombf(lm1, coef, g, prior.mode, baseDensity='normal', nu=3, theta0,
logbf=FALSE, B=10^5)
imombf(lm1, coef, g, prior.mode, nu = 1, theta0 , method='adapt',
nquant=100, B = 10^5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lm1}{ Linear model fit, as returned by \code{lm1}. }
  \item{coef}{ Vector with indexes of coefficients to be
    tested. e.g. \code{coef==c(2,3)}
    and \code{theta0==c(0,0)} tests \code{coef(lm1)[2]=coef(lm1)[3]=0}. }
  \item{g}{ Vector with prior parameter values. See \code{dmom} and
    \code{dimom} for details. }
  \item{prior.mode}{ If specified, \code{g} is set such that the prior
    mode is \code{prior.mode} }
  \item{baseDensity}{Density upon which the Mom prior is
    based. \code{baseDensity=='normal'} results in the normal Mom prior,
  \code{baseDensity=='t'} in the t Mom prior with \code{nu} degrees of freedom.}
  \item{nu}{For \code{mombf}, \code{nu} specifies the degrees of freedom
    of the t Mom prior. It is ignored unless
    \code{baseDensity=='t'}. \code{nu} defaults to 3.
    For \code{imombf}, \code{nu} specifies the degrees of freedom for
    the inverse moment prior (see
    \code{dimom} for details). Defaults to \code{nu=1}, which Cauchy-like
    tails.}
  \item{theta0}{ Null value for the regression coefficients. Defaults to
    0. }
  \item{logbf}{ If \code{logbf==TRUE} the natural logarithm of the Bayes
    factor is returned.}
  \item{method}{ Numerical integration method to compute the bivariate
    integral (only used by \code{imombf}).
    For \code{method=='adapt'}, the inner integral is evaluated (via \code{integrate}) at a series of
    \code{nquant} quantiles of the residual variance posterior distribution, and then
    averaged as described in Johnson (1992).
    Set \code{method=='MC'} to use Monte Carlo integration.}
  \item{nquant}{ Number of quantiles at which to evaluate the integral
    for known \code{sigma}. Only used if \code{method=='adapt'}.}
  \item{B}{Number of Monte Carlo samples to estimate the T Mom and the inverse moment
    Bayes factor. Only used in \code{mombf} if \code{baseDensity=='t'}. Only used in \code{imombf} if \code{method=='MC'}.}
}
\details{
  These functions actually call \code{momunknown} and
  \code{imomunknown}, but they have a simpler interface.
  See \code{dmom} and \code{dimom} for details on the moment and inverse
  moment priors.
}
\value{
\code{mombf} returns the moment Bayes factor to compare the model where
\code{theta!=theta0}
with the null model where \code{theta==theta0}. Large values favor the
alternative model; small values favor the null.
\code{imombf} returns
inverse moment Bayes factors.
}
\references{ See http://rosselldavid.googlepages.com for technical
  reports.
For details on the quantile integration, see Johnson, V.E. A Technique for Estimating Marginal Posterior Densities in Hierarchical Models
Using Mixtures of Conditional Densities. Journal of the American Statistical Association, Vol. 87, No. 419. (Sep., 1992), pp. 852-860.
}
\author{ David Rossell }
\seealso{\code{\link{nlpMarginal}} for a better interface to
  integrated likelihoods and \code{\link{modelSelection}} to also
  search over the model space}
\examples{
##compute Bayes factor for Hald's data
data(hald)
lm1 <- lm(hald[,1] ~ hald[,2] + hald[,3] + hald[,4] + hald[,5])

# Set g so that interval (-0.2,0.2) has 5\% prior probability
# (in standardized effect size scale)
priorp <- .05; q <- .2
gmom <- priorp2g(priorp=priorp,q=q,prior='normalMom')
gimom <- priorp2g(priorp=priorp,q=q,prior='iMom')

mombf(lm1,coef=2,g=gmom) #moment BF
imombf(lm1,coef=2,g=gimom,B=10^5) #inverse moment BF

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

