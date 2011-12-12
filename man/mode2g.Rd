\name{mode2g}
\alias{mode2g}
\alias{g2mode}
\alias{priorp2g}
\title{ Moment and inverse moment prior elicitation. }
\description{
  \code{mode2g} finds the \code{g} value corresponding to a given prior mode.
  \code{g2mode} finds the prior mode corresponding to a given \code{g} value.
  \code{priorp2g} finds the \code{g} value giving \code{priorp} prior
  probability to the interval (\code{-q},\code{q}).
  \cr
  All routines operate in the standardized effect sizes scale.
}
\usage{
mode2g(prior.mode, prior=c("iMom", "normalMom", "tMom"), nu=1, dim=1)
g2mode(g, prior=c("iMom", "normalMom", "tMom"), nu=1, dim=1)
priorp2g(priorp, q, nu=1, prior=c("iMom", "normalMom", "tMom"))
}
\arguments{
  \item{prior.mode}{Prior mode for the quadratic form
    (theta-theta0)' * solve(Sigma) * (theta-theta0)/(n*g*sigma^2),
     where sigma is the dispersion parameter and Sigma is given by the
     design matrix.  }
  \item{prior}{\code{prior=='normalMom'} does computations for the normal
    moment prior, \code{prior=='tMom'} for the T moment prior,
    \code{prior=='iMom'} does computations for the inverse moment prior.
    Currently \code{prior=='tMom'} is not implemented in \code{priorp2g}.}
  \item{nu}{ Prior degrees of freedom for the T moment prior or the iMom
    prior (ignored if \code{prior=='normalMom'}).}
  \item{dim}{ Dimensionality of the parameter, i.e. \code{dim==1} for
    univariate, \code{dim==2} for bivariate and so on. }
  \item{g}{Prior parameter. See \code{dimom} for details.}
  \item{priorp}{ \code{priorp2g} returns g giving \code{priorp} prior
    probability to the interval \code{(-q,q)}. }
  \item{q}{ \code{priorp2g} returns g giving \code{priorp} prior
    probability to the interval \code{(-q,q)}.}
}
\details{
  See \code{dmom} and \code{dimom} for details on the meaning of the
  prior parameters.
}
\value{
  \code{mode2g} returns the value of the prior parameter \code{g} matching
  the given mode.

  \code{g2mode} returns the prior mode for a given prior parameter \code{g}.

  \code{priorp2g} returns g giving \code{priorp} prior probability to the
  interval \code{(-q,q)}.
}
\references{
 See http://rosselldavid.googlepages.com for technical reports.
}
\author{
  David Rossell \email{rosselldavid@gmail.com}
}
\seealso{
  \code{\link{dmom}},
  \code{\link{dimom}},
  \code{\link{mombf}},
  \code{\link{imombf}}
}
\examples{
#find g value giving a prior mode for (theta/(sigma*n*Sigma))^2 at 0.2^2
data(hald)
lm1 <- lm(hald[, 1] ~ hald[, 2] + hald[, 3] + hald[, 4] + hald[, 5])
prior.mode <- .2
gmom <- mode2g(prior.mode^2, prior='normalMom')
gtmom <- mode2g(prior.mode^2, prior='tMom', nu=3)
gimom <- mode2g(prior.mode^2, prior='iMom')
gmom
gimom

#find g value giving 0.05 probability to interval (-.2,.2)
priorp <- .05; q <- .2
gmom <- priorp2g(priorp=priorp, q=q, prior='normalMom')
gimom <- priorp2g(priorp=priorp, q=q, prior='iMom')
gmom
gimom
}
\keyword{distribution}
\keyword{htest}

